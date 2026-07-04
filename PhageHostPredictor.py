import streamlit as st
import pandas as pd
from Bio import SeqIO
import subprocess
import os
import plotly.express as px
from fpdf import FPDF

st.set_page_config(layout="wide", page_title="BioAnalyzer PRO v3.0")


# --- Funkcje Bioinformatyczne ---
def kategoryzuj_bialko(opis):
    opis = str(opis).lower()
    if any(x in opis for x in ["lysin", "endolysin", "holin"]): return "Lizyny i Holiny"
    if any(x in opis for x in ["tail fiber", "tail spike", "receptor binding", "rbp"]): return "Adhezyny (RBP)"
    if any(x in opis for x in ["capsid", "coat", "head"]): return "Kapsyd"
    if "tail" in opis: return "Białka ogonka"
    if "polymerase" in opis: return "Polimerazy"
    if "depolymerase" in opis: return "Depolimerazy"
    if "portal" in opis: return "Białka portalu"
    if "integrase" in opis: return "Integrazy"
    return "Pozostałe"


TAXA_MAP = {
    "coli": "Escherichia coli",
    "aureus": "Staphylococcus aureus",
    "subtilis": "Bacillus subtilis",
    "pyogenes": "Streptococcus pyogenes",
    "typhimurium": "Salmonella typhimurium",
    "pneumoniae": "Klebsiella pneumoniae",
    "aeruginosa": "Pseudomonas aeruginosa"
}


# --- Funkcja Pomocnicza (Usuwanie polskich znaków dla FPDF) ---
def usun_pl(tekst):
    zamienniki = {'ą': 'a', 'ć': 'c', 'ę': 'e', 'ł': 'l', 'ń': 'n', 'ó': 'o', 'ś': 's', 'ź': 'z', 'ż': 'z',
                  'Ą': 'A', 'Ć': 'C', 'Ę': 'E', 'Ł': 'L', 'Ń': 'N', 'Ó': 'O', 'Ś': 'S', 'Ź': 'Z', 'Ż': 'Z'}
    for pl, asc in zamienniki.items():
        tekst = tekst.replace(pl, asc)
    return tekst


# --- Interfejs ---
st.title("🧬 BioAnalyzer PRO: Profesjonalna Analiza Genomowa")
uploaded_file = st.sidebar.file_uploader("Wgraj genom (FASTA)", type=["fasta", "fa"])

if uploaded_file:
    with open("input.fasta", "wb") as f:
        f.write(uploaded_file.getbuffer())

    tab1, tab2, tab3, tab4 = st.tabs(
        ["1. Wykrywanie ORF i Mapa", "2. Adnotacja i Funkcje", "3. Predykcja Hosta", "4. Raport PDF"])

    # ZAKŁADKA 1: ORF i WIZUALIZACJA
    with tab1:
        if st.button("Uruchom wykrywanie ORF (ATG, GTG, TTG)"):
            with st.spinner("Prodigal analizuje..."):
                exe_path = os.path.join(os.getcwd(), "prodigal.exe")
                # Flaga -g 11 wymusza uwzględnianie kodonów start ATG, GTG, TTG
                subprocess.run(
                    [exe_path, "-i", "input.fasta", "-a", "proteins.faa", "-o", "genes.gff", "-q", "-g", "11"],
                    check=True)

                records = list(SeqIO.parse("proteins.faa", "fasta"))
                data = []
                for r in records:
                    parts = r.description.split(" # ")
                    data.append({
                        "ID": r.id,
                        "Start": int(parts[1]),
                        "Stop": int(parts[2]),
                        "Nić": "+" if parts[3] == "1" else "-",
                        "Długość (aa)": len(r.seq)
                    })
                st.session_state["orf_df"] = pd.DataFrame(data)

        if "orf_df" in st.session_state:
            df_plot = st.session_state["orf_df"].copy()
            st.subheader("Interaktywna Mapa Genomu")
            st.info(
                "💡 Użyj kółka myszy, aby przybliżyć (zoom), oraz najedź kursorem na bloki, aby zobaczyć szczegóły ORF.")

            # Przygotowanie danych do Plotly
            df_plot["Długość_nt"] = df_plot["Stop"] - df_plot["Start"]
            df_plot["Contig"] = "Genom"

            # Generowanie wykresu horyzontalnego
            fig = px.bar(df_plot, base="Start", x="Długość_nt", y="Contig", color="Nić",
                         orientation='h', hover_name="ID",
                         hover_data={"Contig": False, "Start": True, "Stop": True, "Długość_nt": False,
                                     "Długość (aa)": True},
                         color_discrete_map={"+": "royalblue", "-": "indianred"})

            st.plotly_chart(fig, use_container_width=True)
            st.dataframe(st.session_state["orf_df"], use_container_width=True)

    # ZAKŁADKA 2: ADNOTACJA
    with tab2:
        if st.button("Uruchom Adnotację"):
            with st.spinner("DIAMOND przeszukuje Swiss-Prot..."):
                diamond_path = os.path.join(os.getcwd(), "diamond.exe")
                cmd = [diamond_path, "blastp", "-d", "sprot_db.dmnd", "-q", "proteins.faa", "-o", "annot.tsv",
                       "-k", "1", "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch",
                       "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle"]
                subprocess.run(cmd, check=True)

                df = pd.read_csv("annot.tsv", sep="\t",
                                 names=["query", "subject", "pident", "length", "mismatch", "gap", "qstart", "qend",
                                        "sstart", "send", "evalue", "bitscore", "stitle"])
                df["evalue"] = df["evalue"].map('{:.2e}'.format)
                df["Kategoria"] = df["subject"].apply(kategoryzuj_bialko)
                st.session_state["annot_df"] = df

        if "annot_df" in st.session_state:
            col1, col2 = st.columns([1, 2])
            with col1:
                st.subheader("Podsumowanie funkcjonalne")
                st.bar_chart(st.session_state["annot_df"]["Kategoria"].value_counts())
            with col2:
                st.subheader("Szczegóły adnotacji")
                st.dataframe(st.session_state["annot_df"][["query", "subject", "stitle", "Kategoria", "evalue"]],
                             use_container_width=True)

            annotated_ids = st.session_state["annot_df"]["query"].unique()
            missed = st.session_state["orf_df"][~st.session_state["orf_df"]["ID"].isin(annotated_ids)]
            st.warning(f"Liczba białek bez dopasowania (ORFans): {len(missed)}")
            with st.expander("Zobacz listę nieprzypisanych ORF"):
                st.write(missed)

    # ZAKŁADKA 3: PREDYKCJA HOSTA
    with tab3:
        st.subheader("Predykcja Hosta")
        if "annot_df" in st.session_state:
            scores = {}
            for title in st.session_state["annot_df"]["stitle"]:
                for key, full_name in TAXA_MAP.items():
                    if key in str(title).lower():
                        scores[full_name] = scores.get(full_name, 0) + 1

            if scores:
                st.session_state["host_scores"] = scores
                total = sum(scores.values())
                cols = st.columns(len(scores))
                for i, (host, pts) in enumerate(scores.items()):
                    cols[i].metric(label=host, value=f"{(pts / total) * 100:.1f}%")
            else:
                st.session_state["host_scores"] = {}
                st.error("Brak rozpoznawalnych nazw gospodarzy w wynikach adnotacji.")
        else:
            st.info("Najpierw wykonaj adnotację w zakładce nr 2.")

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
                    pdf.cell(0, 10, usun_pl("Prawdopodobny gospodarz (Wyniki):"), ln=True)
                    pdf.set_font("Arial", '', 12)

                    scores = st.session_state["host_scores"]
                    total = sum(scores.values())
                    for host, pts in scores.items():
                        prob = (pts / total) * 100
                        pdf.cell(0, 8, usun_pl(f"- {host}: {prob:.1f}% dopasowan"), ln=True)

                # Zapis i przygotowanie pliku do pobrania
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