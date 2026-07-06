# app.py
import streamlit as st
import pandas as pd
from Bio import SeqIO
import subprocess
import os
import tempfile
import plotly.graph_objects as go
import plotly.express as px

# Importy z modułów zewnętrznych
from report import generate_pdf_report
from host_database import load_host_database
from protein_classifier import kategoryzuj_bialko, CATEGORY_COLORS
from src.evidence_engine import EvidenceEngine

# GUI
st.set_page_config(page_title="Fagometr", page_icon="🧬", layout="wide")

PROFESSIONAL_COLORS = [
    "#1f77b4", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
    "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#ff7f0e"
]

st.markdown("""
<style>
    /* Zaokrąglone przyciski i płynne animacje */
    .stButton > button {
        border-radius: 8px;
        transition: 0.2s ease-in-out;
    }
    .stButton > button:hover {
        box-shadow: 0 4px 8px rgba(0,0,0,0.1);
        transform: translateY(-2px);
    }
    /* Delikatnie zaokrąglony pasek boczny */
    [data-testid="stSidebar"] {
        border-right: 1px solid #e0e0e0;
        border-radius: 0 15px 15px 0;
    }
    /* Estetyczne ramki dla tabel i statusów */
    [data-testid="stTable"], [data-testid="stDataFrame"] {
        border-radius: 10px;
        overflow: hidden;
        box-shadow: 0 1px 4px rgba(0,0,0,0.05);
    }
    /* Akcenty zakładek */
    .stTabs [data-baseweb="tab-list"] {
        gap: 10px;
    }
    .stTabs [data-baseweb="tab"] {
        border-radius: 8px 8px 0 0;
    }
</style>
""", unsafe_allow_html=True)

# Zmienne globalne
FILE_NAMES = {
    "fasta": "input.fasta",
    "proteins": "proteins.faa",
    "gff": "genes.gff",
    "annot": "annot.tsv"
}


# FUNKCJE POMOCNICZE
def has_orfs():
    return "orf_df" in st.session_state

def has_annot():
    return "annot_df" in st.session_state

if "logs" not in st.session_state:
    st.session_state["logs"] = []

def log_step(msg):
    st.session_state["logs"].append(msg)
    st.toast(msg)

@st.cache_resource(ttl=3600)
def get_engine():
    db_path = "data/processed/phage_host_database.parquet"
    if os.path.exists(db_path):
        return EvidenceEngine(db_path=db_path)
    return None

@st.cache_resource(ttl=3600)
def safe_load_host_database():
    try:
        return load_host_database()
    except Exception as e:
        st.error(f"⚠️ Błąd podczas ładowania bazy gospodarzy: {e}")
        return None, []

# 4. timeout zabezpieczający przed zawieszeniem programów
def run_command(cmd, task_name, timeout=3600):
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=timeout)
        return True
    except subprocess.TimeoutExpired:
        st.error(f"🔴 Przekroczono limit czasu ({timeout}s) dla {task_name}. Proces zawiesił się i został przerwany.")
        return False
    except subprocess.CalledProcessError as e:
        st.error(f"🔴 Błąd wykonania {task_name}.")
        with st.expander("Szczegóły błędu (stderr)"):
            st.code(e.stderr)
        return False
    except FileNotFoundError:
        st.error(f"🔴 Nie znaleziono programu: {cmd[0]}. Upewnij się, że jest zainstalowany i dostępny w PATH.")
        return False

# GUI

HOST_DB_DF, KNOWN_GENERA = safe_load_host_database()

with st.sidebar:
    st.title("⚙️ Opcje i Status")
    if HOST_DB_DF is not None:
        st.success(f"✅ Baza aktywna ({len(KNOWN_GENERA)} hostów)")
    else:
        st.warning("⚠️ Tryb ograniczony (Brak bazy)")

    uploaded_file = st.file_uploader("Wgraj genom (FASTA)", type=["fasta", "fa"])

st.title("🧬 Fagometr: Profesjonalna Analiza Genomowa")

if uploaded_file:
    # Wykrywanie zmiany pliku i reset stanu
    file_id = uploaded_file.file_id
    if "current_file_id" not in st.session_state or st.session_state["current_file_id"] != file_id:
        st.session_state["current_file_id"] = file_id

        # Usuwamy stare wyniki
        for key in ["orf_df", "annot_df", "host_results"]:
            st.session_state.pop(key, None)
        st.session_state["logs"] = []

        # Samooczyszczający się katalog za pomocą TemporaryDirectory
        if "temp_dir_obj" in st.session_state:
            st.session_state["temp_dir_obj"].cleanup()

        st.session_state["temp_dir_obj"] = tempfile.TemporaryDirectory(prefix="fagometr_")
        st.session_state["work_dir"] = st.session_state["temp_dir_obj"].name

    work_dir = st.session_state["work_dir"]
    fasta_path = os.path.join(work_dir, FILE_NAMES["fasta"])

    # Zapis nowego pliku
    with open(fasta_path, "wb") as f:
        f.write(uploaded_file.getbuffer())

    tab1, tab2, tab3, tab4 = st.tabs([
        "1. Wykrywanie ORF",
        "2. Adnotacja DIAMOND",
        "3. Predykcja Gospodarza",
        "4. Raport PDF"
    ])

    # ZAKŁADKA 1: ORF
    with tab1:
        if st.button("Uruchom Prodigal"):
            proteins_path = os.path.join(work_dir, FILE_NAMES["proteins"])
            gff_path = os.path.join(work_dir, FILE_NAMES["gff"])

            cmd = ["prodigal.exe", "-i", fasta_path, "-a", proteins_path, "-o", gff_path, "-q"]

            with st.status("Prodigal analizuje genom...", expanded=True) as status:
                st.write("Uruchamianie narzędzia...")
                if run_command(cmd, "Prodigal", timeout=600):  # Krótszy timeout dla Prodigala
                    if os.path.exists(proteins_path) and os.path.getsize(proteins_path) > 0:
                        records = list(SeqIO.parse(proteins_path, "fasta"))
                        data = []
                        for r in records:
                            parts = r.description.split(" # ")
                            if len(parts) >= 4:
                                data.append({
                                    "ID": r.id,
                                    "Start": int(parts[1]),
                                    "Stop": int(parts[2]),
                                    "Nić": "+" if parts[3] == "1" else "-",
                                    "Długość (aa)": len(r.seq)
                                })

                        st.session_state["orf_df"] = pd.DataFrame(data)
                        log_step(f"Zakończono: Znaleziono {len(data)} ORF.")
                        status.update(label="Analiza zakończona sukcesem!", state="complete", expanded=False)
                    else:
                        status.update(label="Błąd: Brak wyników", state="error")
                        st.error("Prodigal nie wygenerował pliku wyjściowego.")
                else:
                    status.update(label="Błąd analizy", state="error")

        if has_orfs():
            st.subheader("Interaktywna Mapa Genomu")
            fig = go.Figure()
            max_stop = st.session_state["orf_df"]["Stop"].max()
            fig.add_trace(go.Scatter(
                x=[0, max_stop], y=[0, 0], mode="lines",
                line=dict(color="black", width=2), hoverinfo="skip", showlegend=False
            ))

            arrow_len = max_stop * 0.015
            has_ann = has_annot()

            for _, row in st.session_state["orf_df"].iterrows():
                start, stop, strand, g_id = row["Start"], row["Stop"], row["Nić"], row["ID"]
                y_base = 0.4 if strand == "+" else -0.4
                color = "#1f77b4" if strand == "+" else "#d62728"
                cat_text = ""

                if has_ann:
                    annot_match = st.session_state["annot_df"][st.session_state["annot_df"]["query"] == g_id]
                    if not annot_match.empty:
                        kategoria = annot_match.iloc[0]["Kategoria"]
                        color = CATEGORY_COLORS.get(kategoria, CATEGORY_COLORS.get("Pozostałe", "#808080"))
                        cat_text = f"<br>Kategoria: {kategoria}"

                al = min(arrow_len, (stop - start) * 0.4)
                if strand == '+':
                    x_poly = [start, stop - al, stop, stop - al, start, start]
                else:
                    x_poly = [stop, start + al, start, start + al, stop, stop]

                y_poly = [y_base - 0.25, y_base - 0.25, y_base, y_base + 0.25, y_base + 0.25, y_base - 0.25]
                hover_text = f"<b>{g_id}</b><br>Start: {start}<br>Stop: {stop}<br>Nić: {strand}<br>Długość: {row['Długość (aa)']} aa{cat_text}"

                fig.add_trace(go.Scatter(
                    x=x_poly, y=y_poly, fill="toself", fillcolor=color,
                    line=dict(color="black", width=1), name=g_id,
                    text=hover_text, hoverinfo="text", showlegend=False
                ))

            fig.update_layout(
                yaxis=dict(showticklabels=False, range=[-1, 1], zeroline=False),
                xaxis=dict(title="Pozycja (nt)", showgrid=True),
                height=350, plot_bgcolor="white",
                margin=dict(l=10, r=10, t=30, b=30)
            )
            st.plotly_chart(fig, use_container_width=True)
            st.dataframe(st.session_state["orf_df"], use_container_width=True)

    # ZAKŁADKA 2: DIAMOND
    with tab2:
        st.info("Wymaga przeprowadzenia detekcji ORF (Zakładka 1)." if not has_orfs() else "Gotowy do adnotacji.")

        if st.button("Uruchom Adnotację DIAMOND", disabled=not has_orfs()):
            proteins_path = os.path.join(work_dir, FILE_NAMES["proteins"])
            annot_path = os.path.join(work_dir, FILE_NAMES["annot"])

            cmd = [
                "diamond.exe", "blastp", "-d", "phage_db.dmnd", "-q", proteins_path,
                "-o", annot_path, "--outfmt", "6", "qseqid", "sseqid", "pident",
                "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
                "evalue", "bitscore", "qlen", "slen", "stitle", "-k", "50"
            ]

            with st.status("DIAMOND szuka homologów...", expanded=True) as status:
                st.write("Wykonywanie dopasowania...")
                # Timeout dla DIAMOND-a to ważna sprawa, bo bywa kapryśny
                if run_command(cmd, "DIAMOND", timeout=3600):
                    if os.path.exists(annot_path) and os.path.getsize(annot_path) > 0:
                        df = pd.read_csv(annot_path, sep="\t", names=[
                            "query", "subject", "pident", "length", "mismatch", "gapopen",
                            "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "stitle"
                        ])

                        df["stitle"] = df["stitle"].fillna("")
                        df["Kategoria"] = df["stitle"].apply(kategoryzuj_bialko)

                        st.session_state["annot_df"] = df
                        log_step(f"Zakończono adnotację: Zidentyfikowano {len(df)} trafień.")
                        status.update(label="Adnotacja zakończona", state="complete", expanded=False)
                    else:
                        status.update(label="Brak trafień", state="error")
                        st.warning("DIAMOND nie znalazł żadnych homologów.")
                else:
                    status.update(label="Błąd analizy", state="error")

        if has_annot():
            st.markdown("### 🔬 Kategorie funkcjonalne")
            counts = st.session_state["annot_df"]["Kategoria"].value_counts()

            fig_pie = px.pie(
                values=counts.values,
                names=counts.index,
                hole=0.4,
                color_discrete_sequence=PROFESSIONAL_COLORS
            )
            fig_pie.update_traces(
                textposition='inside',
                textinfo='percent+label',
                marker=dict(line=dict(color='#ffffff', width=2))
            )
            fig_pie.update_layout(
                showlegend=False,
                margin=dict(t=150, b=150, l=150, r=150),
                height=700,
                autosize=True
            )
            st.plotly_chart(fig_pie, use_container_width=True)
            df_display = st.session_state["annot_df"].copy()
            df_display["evalue"] = df_display["evalue"].apply(lambda x: "{:.2e}".format(x))
            st.dataframe(df_display, width='stretch')

    # ZAKŁADKA 3: HOST
    with tab3:
        st.info("Wymaga przeprowadzenia adnotacji (Zakładka 2)." if not has_annot() else "Gotowy do predykcji.")

        if st.button("🚀 Uruchom Predykcję Gospodarza", type="primary", disabled=not has_annot()):
            with st.spinner("Przetwarzanie dowodów biologicznych..."):
                engine = get_engine()
                if engine:
                    try:
                        results = engine.predict_host(st.session_state["annot_df"], st.session_state["orf_df"])
                        st.session_state["host_results"] = results
                        log_step("Predykcja gospodarza zakończona sukcesem.")
                    except Exception as e:
                        st.error(f"Wystąpił błąd podczas predykcji: {e}")
                else:
                    st.error("Błąd: Nie znaleziono pliku bazy predykcyjnej (phage_host_database.parquet).")

        if "host_results" in st.session_state:
            results = st.session_state["host_results"]
            if isinstance(results, list) and len(results) > 0:
                st.subheader("📊 Rozkład prawdopodobieństwa")
                df_plot = pd.DataFrame(results)
                total_score = df_plot['Score'].sum()
                if total_score > 0:
                    df_plot['Procent'] = (df_plot['Score'] / total_score) * 100
                else:
                    df_plot['Procent'] = 0
                fig_host = px.bar(
                    df_plot.head(10),
                    x='Procent',
                    y='Host',
                    orientation='h',
                    color_discrete_sequence=[PROFESSIONAL_COLORS[0]],
                    title="Top Kandydaci na Gospodarza (w %)"
                )
                fig_host.update_traces(texttemplate='%{x:.1f}%', textposition='outside')

                fig_host.update_layout(yaxis={'categoryorder': 'total ascending'})
                fig_host.update_xaxes(title_text='Prawdopodobieństwo (%)')
                st.plotly_chart(fig_host, use_container_width=True)

                st.subheader("📋 Lista kandydatów")
                res_data = [
                    {
                        "Gospodarz (Host)": r["Host"],
                        "Zgodność": f"{(r['Score'] / total_score) * 100:.1f}%" if total_score > 0 else "0%",
                        "E-value": f"{r.get('E-value', 0):.2e}"
                    }
                    for r in results
                ]
                st.table(pd.DataFrame(res_data))
            else:
                st.warning("Brak znaczących wyników predykcji. Prawdopodobnie brak mocnych homologów.")

    # --- ZAKŁADKA 4: RAPORT ---
    with tab4:
        st.markdown("### Generowanie podsumowania PDF")

        if not has_annot() or "host_results" not in st.session_state:
            st.warning("Analiza nie jest kompletna. Raport może nie zawierać wszystkich danych.")

        if st.button("📄 Generuj PDF", disabled=not has_orfs()):
            with st.spinner("Trwa składanie raportu..."):
                try:
                    pdf_path = generate_pdf_report(
                        uploaded_file.name,
                        st.session_state.get("orf_df"),
                        st.session_state.get("annot_df"),
                        st.session_state.get("host_results")
                    )

                    with open(pdf_path, "rb") as f:
                        st.download_button(
                            label="⬇️ Pobierz Raport",
                            data=f,
                            file_name=f"Raport_Fagometr_{uploaded_file.name}.pdf",
                            mime="application/pdf"
                        )
                    log_step("Raport wygenerowany pomyślnie.")
                except Exception as e:
                    st.error(f"🔴 Nie udało się wygenerować raportu: {e}")

# SEKCJA LOGÓW
st.markdown("---")
with st.expander("📜 Logi systemowe operacji"):
    if not st.session_state["logs"]:
        st.write("Brak operacji w tej sesji.")
    else:
        for log in reversed(st.session_state["logs"]):
            st.text(f"• {log}")

# Stopka
st.markdown("---")
st.markdown("""
<div style="text-align: center; color: #666; padding: 10px 0; font-family: sans-serif;">
    <p style="margin: 0; font-size: 1.1em;"><b>Fagometr</b> | wersja 2.0</p>
    <p style="margin: 0; font-size: 0.9em;">&copy; 2026</p>
    <p style="margin: 0; font-size: 0.9em;">Autor: Michał Chojnacki</p>
</div>
""", unsafe_allow_html=True)