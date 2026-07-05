# app.py
import streamlit as st
import pandas as pd
from Bio import SeqIO
import subprocess
import os
import plotly.graph_objects as go
import plotly.express as px
from report import generate_pdf_report

# Importy z modułów
from host_database import load_host_database
from protein_classifier import kategoryzuj_bialko, CATEGORY_COLORS
from src.evidence_engine import EvidenceEngine

# KONFIGURACJA STRONY
st.set_page_config(layout="wide")

# INICJALIZACJA SILNIKA
@st.cache_resource
def get_engine():
    if os.path.exists("data/processed/phage_host_database.parquet"):
        return EvidenceEngine(db_path="data/processed/phage_host_database.parquet")
    return None
TAX_DF = pd.read_csv("taxonomy_final.csv")

# System logów
if "logs" not in st.session_state: st.session_state["logs"] = []
def add_log(msg): st.session_state["logs"].append(msg)

# Inicjalizacja Danych
HOST_DB_DF, KNOWN_GENERA = load_host_database()

with st.sidebar:
    st.title("⚙️ Opcje i Status")
    if HOST_DB_DF is not None:
        st.success(f"✅ Baza aktywna ({len(KNOWN_GENERA)} hostów)")
    else:
        st.warning("⚠️ Tryb ograniczony")

    uploaded_file = st.file_uploader("Wgraj genom (FASTA)", type=["fasta", "fa"])

st.title("🧬 BioAnalyzer PRO: Profesjonalna Analiza Genomowa")

if uploaded_file:
    with open("input.fasta", "wb") as f:
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
            with st.spinner("Prodigal analizuje genom..."):
                add_log("Uruchomiono Prodigal...")
                subprocess.run(["prodigal.exe", "-i", "input.fasta", "-a", "proteins.faa", "-o", "genes.gff", "-q"],
                               check=True)
                records = list(SeqIO.parse("proteins.faa", "fasta"))
                data = [{"ID": r.id, "Start": int(r.description.split(" # ")[1]),
                         "Stop": int(r.description.split(" # ")[2]),
                         "Nić": "+" if r.description.split(" # ")[3] == "1" else "-", "Długość (aa)": len(r.seq)} for r
                        in records]
                st.session_state["orf_df"] = pd.DataFrame(data)
                add_log("ORF wykryte pomyślnie.")

        if "orf_df" in st.session_state:
            st.subheader("Interaktywna Mapa Genomu")
            fig = go.Figure()
            max_stop = st.session_state["orf_df"]["Stop"].max()
            fig.add_trace(
                go.Scatter(x=[0, max_stop], y=[0, 0], mode="lines", line=dict(color="black", width=2), hoverinfo="skip",
                           showlegend=False))

            arrow_len = max_stop * 0.015
            has_annot = "annot_df" in st.session_state

            for _, row in st.session_state["orf_df"].iterrows():
                start, stop, strand, g_id = row["Start"], row["Stop"], row["Nić"], row["ID"]
                y_base = 0.4 if strand == "+" else -0.4
                color = "#1f77b4" if strand == "+" else "#d62728"
                cat_text = ""

                if has_annot:
                    annot_match = st.session_state["annot_df"][st.session_state["annot_df"]["query"] == g_id]
                    if not annot_match.empty:
                        kategoria = annot_match.iloc[0]["Kategoria"]
                        color = CATEGORY_COLORS.get(kategoria, CATEGORY_COLORS["Pozostałe"])
                        cat_text = f"<br>Kategoria: {kategoria}"

                al = min(arrow_len, (stop - start) * 0.4)
                x_poly = [start, stop - al, stop, stop - al, start, start] if strand == '+' else [stop, start + al,
                                                                                                  start, start + al,
                                                                                                  stop, stop]
                y_poly = [y_base - 0.25, y_base - 0.25, y_base, y_base + 0.25, y_base + 0.25, y_base - 0.25]
                hover_text = f"<b>{g_id}</b><br>Start: {start}<br>Stop: {stop}<br>Nić: {strand}<br>Długość: {row['Długość (aa)']} aa{cat_text}"
                fig.add_trace(
                    go.Scatter(x=x_poly, y=y_poly, fill="toself", fillcolor=color, line=dict(color="black", width=1),
                               name=g_id, text=hover_text, hoverinfo="text", showlegend=False))

            fig.update_layout(yaxis=dict(showticklabels=False, range=[-1, 1], zeroline=False),
                              xaxis=dict(title="Pozycja (nt)", showgrid=True), height=350, plot_bgcolor="white",
                              margin=dict(l=10, r=10, t=30, b=30))
            st.plotly_chart(fig, use_container_width=True)
            st.dataframe(st.session_state["orf_df"], use_container_width=True)

    # ZAKŁADKA 2: DIAMOND
    with tab2:
            if st.button("Uruchom Adnotację DIAMOND"):
                with st.spinner("DIAMOND szuka homologów..."):
                    add_log("Uruchomiono DIAMOND...")

                    cmd = [
                        "diamond.exe", "blastp", "-d", "phage_db.dmnd", "-q", "proteins.faa",
                        "-o", "annot.tsv", "--outfmt", "6", "qseqid", "sseqid", "pident",
                        "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
                        "evalue", "bitscore", "qlen", "slen", "stitle", "-k", "50"
                    ]
                    subprocess.run(cmd, check=True)

                    df = pd.read_csv("annot.tsv", sep="\t", names=[
                        "query", "subject", "pident", "length", "mismatch", "gapopen",
                        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "stitle"
                    ])

                    df["stitle"] = df["stitle"].fillna("")
                    df["Kategoria"] = df["stitle"].apply(kategoryzuj_bialko)

                    st.session_state["annot_df"] = df
                    add_log("Adnotacja zakończona.")

            if "annot_df" in st.session_state:
                col1, col2 = st.columns([1, 2])
                with col1:
                    st.markdown("### 🥧 Podział funkcjonalny")
                    counts = st.session_state["annot_df"]["Kategoria"].value_counts()
                    fig_pie = px.pie(values=counts.values, names=counts.index, hole=0.3)
                    st.plotly_chart(fig_pie, use_container_width=True)
                with col2:
                    st.dataframe(st.session_state["annot_df"], use_container_width=True)

    # ZAKŁADKA 3: HOST
    with tab3:
        if st.button("🚀 Uruchom Predykcję Gospodarza", type="primary"):
            if "annot_df" in st.session_state:
                add_log("Predykcja gospodarza w toku...")
                engine = get_engine()
                if engine:
                    results = engine.predict_host(st.session_state["annot_df"], st.session_state["orf_df"])
                else:
                    st.error("Nie znaleziono silnika predykcji.")
                    results = []
                st.session_state["host_results"] = results
                add_log(f"Predykcja zakończona.")
            else:
                st.error("Uruchom najpierw zakładki 1 i 2.")

        if "host_results" in st.session_state:
            results = st.session_state["host_results"]
            if isinstance(results, list) and len(results) > 0:
                # 1. Wykres (nad tabelą, pełna szerokość)
                st.subheader("📊 Rozkład prawdopodobieństwa")
                df_plot = pd.DataFrame(results)
                fig_host = px.pie(df_plot, values='Score', names='Host', hole=0.4)
                fig_host.update_layout(margin=dict(t=30, b=30, l=30, r=30))
                st.plotly_chart(fig_host, use_container_width=True)

                # 2. Tabela (pod wykresem, bez Evidence)
                st.subheader("📋 Lista kandydatów")
                res_data = [
                    {
                        "Host": r["Host"],
                        "Score": f"{r['Score']:.2f}",
                        "E-value": f"{r['E-value']:.2e}"
                    }
                    for r in results
                ]
                st.table(pd.DataFrame(res_data))
            else:
                st.info("Brak wyników predykcji.")

    # ZAKŁADKA 4: RAPORT
    with tab4:
        if st.button("Generuj PDF"):
            pdf_path = generate_pdf_report(uploaded_file.name, st.session_state["orf_df"],
                                           st.session_state.get("annot_df"), st.session_state.get("host_results"))
            with open(pdf_path, "rb") as f: st.download_button("📄 Pobierz Raport", f, "Raport.pdf", "application/pdf")
            add_log("Raport wygenerowany.")

# SEKCJA LOGÓW NA DOLE
st.markdown("---")
with st.expander("📜 Logi operacji"):
    for log in reversed(st.session_state["logs"]):
        st.text(f"• {log}")