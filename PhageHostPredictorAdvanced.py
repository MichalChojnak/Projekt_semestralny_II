import streamlit as st
import pandas as pd
from Bio import SeqIO
import subprocess
import os
import plotly.graph_objects as go
from fpdf import FPDF

st.set_page_config(layout="wide", page_title="BioAnalyzer PRO v6.1")


# --- MODUŁ BAZY WIEDZY ---
@st.cache_data
def load_host_database():
    db_path = "host_database.csv"
    if os.path.exists(db_path):
        try:
            df = pd.read_csv(db_path)
            # Pobieramy unikalne nazwy rodzajów bakterii (np. Escherichia, Pseudomonas)
            genera = set(df['host_genus'].dropna().str.capitalize().unique())
            # Oczyszczamy z przypadkowych słów, które mogły się pobrać z NCBI
            banned_words = {"Phage", "Bacteriophage", "Virus", "Unidentified", "Unknown"}
            genera = {g for g in genera if len(g) > 2 and g not in banned_words}
            return df, genera
        except Exception as e:
            return None, set()
    return None, set()


# Ładowanie bazy przy starcie
HOST_DB_DF, KNOWN_GENERA = load_host_database()

# --- Moduł 2: Protein Importance ---
PROTEIN_WEIGHTS = {
    "Receptor Binding Protein (RBP)": 12,
    "Tail Fiber": 10,
    "Tail Spike": 10,
    "Depolimerazy": 9,
    "Adsorpcja": 8,
    "Baseplate": 6,
    "Lizyny i Holiny": 2,
    "Kapsyd": 0,
    "Białka portalu": 0,
    "Polimerazy": 0,
    "Terminazy": 0,
    "Integrazy": 0,
    "Pozostałe": 1
}

CATEGORY_COLORS = {
    "Receptor Binding Protein (RBP)": "#ff7f0e",
    "Tail Fiber": "#ffbb78",
    "Tail Spike": "#ff9896",
    "Depolimerazy": "#8c564b",
    "Adsorpcja": "#e377c2",
    "Baseplate": "#1f77b4",
    "Lizyny i Holiny": "#d62728",
    "Kapsyd": "#2ca02c",
    "Białka portalu": "#17becf",
    "Polimerazy": "#9467bd",
    "Terminazy": "#c5b0d5",
    "Integrazy": "#f7b6d2",
    "Pozostałe": "#7f7f7f",
    "Brak adnotacji (ORFan)": "#c7c7c7"
}


def kategoryzuj_bialko(opis):
    opis = str(opis).lower()
    if any(x in opis for x in ["receptor binding", "rbp"]): return "Receptor Binding Protein (RBP)"
    if "tail fiber" in opis: return "Tail Fiber"
    if "tail spike" in opis: return "Tail Spike"
    if "depolymerase" in opis: return "Depolimerazy"
    if "adsorption" in opis: return "Adsorpcja"
    if "baseplate" in opis: return "Baseplate"
    if any(x in opis for x in ["lysin", "endolysin", "holin"]): return "Lizyny i Holiny"
    if any(x in opis for x in ["capsid", "coat", "head"]): return "Kapsyd"
    if "portal" in opis: return "Białka portalu"
    if "polymerase" in opis: return "Polimerazy"
    if "terminase" in opis: return "Terminazy"
    if "integrase" in opis: return "Integrazy"
    return "Pozostałe"


def usun_pl(tekst):
    zamienniki = {'ą': 'a', 'ć': 'c', 'ę': 'e', 'ł': 'l', 'ń': 'n', 'ó': 'o', 'ś': 's', 'ź': 'z', 'ż': 'z',
                  'Ą': 'A', 'Ć': 'C', 'Ę': 'E', 'Ł': 'L', 'Ń': 'N', 'Ó': 'O', 'Ś': 'S', 'Ź': 'Z', 'Ż': 'Z'}
    for pl, asc in zamienniki.items():
        tekst = tekst.replace(pl, asc)
    return tekst


# --- Interfejs (Pasek boczny) ---
with st.sidebar:
    st.title("⚙️ Opcje i Status")

    # Wyświetlanie statusu nowej bazy
    st.markdown("### Status Bazy Wiedzy")
    if HOST_DB_DF is not None:
        st.success(
            f"✅ Baza aktywna\n- **{len(HOST_DB_DF):,}** białek fagowych\n- **{len(KNOWN_GENERA)}** unikalnych rodzajów gospodarzy")
    else:
        st.error("❌ Brak pliku `host_database.csv`. Engine będzie działał w trybie ograniczonym.")

    st.markdown("---")
    uploaded_file = st.file_uploader("Wgraj genom (FASTA)", type=["fasta", "fa"])

st.title("🧬 BioAnalyzer PRO: Profesjonalna Analiza Genomowa")

if uploaded_file:
    with open("input.fasta", "wb") as f:
        f.write(uploaded_file.getbuffer())

    tab1, tab2, tab3, tab4 = st.tabs(
        ["1. Wykrywanie ORF i Mapa", "2. Adnotacja i Funkcje", "3. Predykcja Hosta (Evidence V2)", "4. Raport PDF"])

    # ZAKŁADKA 1: ORF i WIZUALIZACJA
    with tab1:
        if st.button("Uruchom wykrywanie ORF (ATG, GTG, TTG)"):
            with st.spinner("Prodigal analizuje..."):
                exe_path = os.path.join(os.getcwd(), "prodigal.exe")
                subprocess.run(
                    [exe_path, "-i", "input.fasta", "-a", "proteins.faa", "-o", "genes.gff", "-q", "-g", "11"],
                    check=True)
                records = list(SeqIO.parse("proteins.faa", "fasta"))

                # ZMODYFIKOWANA LINIJKA - dodane pole "Sekwencja (aa)"
                data = [{"ID": r.id, "Start": int(r.description.split(" # ")[1]),
                         "Stop": int(r.description.split(" # ")[2]),
                         "Nić": "+" if r.description.split(" # ")[3] == "1" else "-", "Długość (aa)": len(r.seq),
                         "Sekwencja (aa)": str(r.seq)} for r in records]

                st.session_state["orf_df"] = pd.DataFrame(data)

        if "orf_df" in st.session_state:
            has_annot = "annot_df" in st.session_state
            st.subheader("Interaktywna Mapa Genomu")

            fig = go.Figure()
            max_stop = st.session_state["orf_df"]["Stop"].max()
            fig.add_trace(
                go.Scatter(x=[0, max_stop], y=[0, 0], mode="lines", line=dict(color="black", width=2), hoverinfo="skip",
                           showlegend=False))

            if has_annot:
                for cat, color in CATEGORY_COLORS.items():
                    fig.add_trace(
                        go.Scatter(x=[None], y=[None], mode='markers', marker=dict(size=10, color=color), name=cat))

            arrow_len = max_stop * 0.015

            for _, row in st.session_state["orf_df"].iterrows():
                start, stop, strand, g_id = row["Start"], row["Stop"], row["Nić"], row["ID"]
                y_base = 0.4 if strand == "+" else -0.4
                color = "rgba(31, 119, 180, 0.8)" if strand == "+" else "rgba(214, 39, 40, 0.8)"
                cat_text = ""

                if has_annot:
                    annot_match = st.session_state["annot_df"][st.session_state["annot_df"]["query"] == g_id]
                    if not annot_match.empty:
                        kategoria = annot_match.iloc[0]["Kategoria"]
                        color = CATEGORY_COLORS.get(kategoria, CATEGORY_COLORS["Pozostałe"])
                        cat_text = f"<br>Kategoria: {kategoria}"
                    else:
                        color = CATEGORY_COLORS["Brak adnotacji (ORFan)"]
                        cat_text = "<br>Kategoria: Brak adnotacji (ORFan)"

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
                              xaxis=dict(title="Pozycja (nt)", showgrid=True), height=400, plot_bgcolor="white",
                              legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1))
            st.plotly_chart(fig, use_container_width=True)

            st.subheader("Tabela zidentyfikowanych ORF")
            st.dataframe(st.session_state["orf_df"], use_container_width=True)

    # ZAKŁADKA 2: ADNOTACJA
    with tab2:
        if st.button("Uruchom Adnotację"):
            with st.spinner("DIAMOND przeszukuje Swiss-Prot..."):
                diamond_path = os.path.join(os.getcwd(), "diamond.exe")
                cmd = [diamond_path, "blastp", "-d", "sprot_db.dmnd", "-q", "proteins.faa", "-o", "annot.tsv", "-k",
                       "1", "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
                       "qend", "sstart", "send", "evalue", "bitscore", "stitle"]
                subprocess.run(cmd, check=True)
                df = pd.read_csv("annot.tsv", sep="\t",
                                 names=["query", "subject", "pident", "length", "mismatch", "gap", "qstart", "qend",
                                        "sstart", "send", "evalue", "bitscore", "stitle"])
                df["evalue_fmt"] = df["evalue"].map('{:.2e}'.format)
                df["Kategoria"] = df["subject"].apply(kategoryzuj_bialko)
                st.session_state["annot_df"] = df

        if "annot_df" in st.session_state:
            st.dataframe(st.session_state["annot_df"][
                             ["query", "subject", "stitle", "Kategoria", "pident", "bitscore", "evalue_fmt"]],
                         use_container_width=True)
        if "annot_df" in st.session_state:
            # --- DODAJ TO: Wykres statystyk kategorii ---
            st.subheader("Statystyki Adnotacji")
            cat_counts = st.session_state["annot_df"]["Kategoria"].value_counts()

            fig_pie = go.Figure(data=[go.Pie(labels=cat_counts.index, values=cat_counts.values, hole=0.3)])
            fig_pie.update_layout(height=400, margin=dict(t=0, b=0, l=0, r=0))
            st.plotly_chart(fig_pie, use_container_width=True)
            # -------------------------------------------
    # ZAKŁADKA 3: PREDYKCJA HOSTA (EVIDENCE ENGINE V2)
    with tab3:
        st.subheader("Wielowymiarowy Host Evidence Engine (V2)")
        if "annot_df" in st.session_state:
            host_evidence = {}
            evidence_details = []

            for _, row in st.session_state["annot_df"].iterrows():
                title, pident, bitscore, kategoria = str(row["stitle"]).lower(), row["pident"], row["bitscore"], row[
                    "Kategoria"]

                # ZMIANA: Wyzerowanie wagi "Pozostałe" zgodnie z Twoim życzeniem
                waga = PROTEIN_WEIGHTS.get(kategoria, 0) if kategoria != "Pozostałe" else 0

                # Obliczanie Confidence Score dla pojedynczego białka
                if waga > 0:
                    confidence = (bitscore / 500) * (pident / 100) * waga

                    found_host = None

                    # SZUKANIE W NOWEJ BAZIE WIEDZY (z obsługą "phage")
                    if KNOWN_GENERA:
                        for genus in KNOWN_GENERA:
                            # Teraz sprawdza czy nazwa bakterii jest w opisie LUB czy jest przed "phage"
                            if genus.lower() in title or f"{genus.lower()} phage" in title:
                                found_host = f"{genus} spp."
                                break

                    if found_host:
                        host_evidence[found_host] = host_evidence.get(found_host, 0) + confidence
                        evidence_details.append(
                            {"Host": found_host, "Białko": kategoria, "Confidence Score": round(confidence, 2),
                             "Opis": row["stitle"]})

            if host_evidence:
                st.session_state["host_scores"] = host_evidence
                total_score = sum(host_evidence.values())

                sorted_hosts = sorted(host_evidence.items(), key=lambda item: item[1], reverse=True)
                top_host, top_score = sorted_hosts[0]
                top_prob = (top_score / total_score) * 100

                st.markdown("═══════════════════════════════════")
                st.markdown(f"### Najbardziej prawdopodobny host")
                st.markdown(f"## 🦠 **{top_host}** ({top_prob:.1f}%)")
                st.progress(int(top_prob))

                st.markdown("**Kluczowe dowody białkowe z Bazy NCBI:**")
                top_evidence = [e for e in evidence_details if e["Host"] == top_host]
                for ev in sorted(top_evidence, key=lambda x: x["Confidence Score"], reverse=True)[:5]:
                    st.markdown(f"✔ Homolog: **{ev['Białko']}** (Score: {ev['Confidence Score']})")
                st.markdown("═══════════════════════════════════")

                st.subheader("Pełny Ranking Kandydatów")
                ranking_data = []
                for host, score in sorted_hosts:
                    prob = (score / total_score) * 100
                    conf_level = "Wysoki" if prob > 70 else ("Średni" if prob > 20 else "Niski")
                    ranking_data.append(
                        {"Host": host, "Prawdopodobieństwo (%)": f"{prob:.1f}%", "Punkty Dowodowe": round(score, 2),
                         "Zaufanie": conf_level})

                st.table(pd.DataFrame(ranking_data))

            else:
                st.session_state["host_scores"] = {}
                st.warning(
                    "Silnik znalazł białka strukturalne, ale nie zdołał powiązać ich z żadnym gospodarzem ze 109-tysięcznej bazy. Fag może infekować rzadką, nieopisaną bakterię.")
        else:
            st.info("Najpierw wykonaj adnotację w Zakładce 2.")

    # ZAKŁADKA 4: GENEROWANIE RAPORTU
    with tab4:
        st.subheader("Eksport Wyników Analizy")
        st.write("Wygeneruj podsumowanie wszystkich operacji w formie pliku PDF.")

        if st.button("Generuj PDF"):
            if "orf_df" not in st.session_state:
                st.error("Najpierw wygeneruj ORF (Zakładka 1)!")
            else:
                pdf = FPDF()
                pdf.add_page()
                pdf.set_font("Arial", 'B', 16)
                pdf.cell(0, 10, "Raport Analizy Genomowej - BioAnalyzer PRO", ln=True, align="C")
                pdf.ln(10)

                pdf.set_font("Arial", '', 12)
                pdf.cell(0, 10, usun_pl(f"Wgrany plik: {uploaded_file.name}"), ln=True)
                pdf.cell(0, 10, usun_pl(f"Liczba zidentyfikowanych ORF: {len(st.session_state['orf_df'])}"), ln=True)

                if HOST_DB_DF is not None:
                    pdf.cell(0, 10, usun_pl(f"Silnik Predykcji: Baza NCBI ({len(HOST_DB_DF):,} bialek)"), ln=True)

                if "annot_df" in st.session_state:
                    pdf.ln(5)
                    pdf.set_font("Arial", 'B', 12)
                    pdf.cell(0, 10, usun_pl("Adnotacja funkcjonalna - Kategorie:"), ln=True)
                    pdf.set_font("Arial", '', 12)

                    counts = st.session_state["annot_df"]["Kategoria"].value_counts()
                    for cat, count in counts.items():
                        pdf.cell(0, 8, usun_pl(f"- {cat}: {count} bialek"), ln=True)

                if "host_scores" in st.session_state and st.session_state["host_scores"]:
                    pdf.ln(5)
                    pdf.set_font("Arial", 'B', 12)
                    pdf.cell(0, 10, usun_pl("Prawdopodobny gospodarz (Wyniki z Evidence V2):"), ln=True)
                    pdf.set_font("Arial", '', 12)

                    scores = st.session_state["host_scores"]
                    total = sum(scores.values())
                    for host, pts in scores.items():
                        prob = (pts / total) * 100
                        pdf.cell(0, 8, usun_pl(f"- {host}: {prob:.1f}% dopasowan"), ln=True)

                pdf.output("raport_bioanalyzer.pdf")
                with open("raport_bioanalyzer.pdf", "rb") as pdf_file:
                    st.download_button(
                        label="Pobierz Raport PDF",
                        data=pdf_file,
                        file_name="Raport_BioAnalyzer.pdf",
                        mime="application/pdf"
                    )
                st.success("Raport został wygenerowany! Kliknij przycisk wyżej, aby go pobrać.")

else:
    st.info("Wgraj plik FASTA w panelu bocznym, aby rozpocząć analizę.")