import streamlit as st
import pandas as pd
from Bio import SeqIO
import subprocess
import os

# Konfiguracja strony
st.set_page_config(layout="wide", page_title="BioAnalyzer PRO v2.0")


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


# Baza mapowania taksonomicznego dla predykcji hosta
TAXA_MAP = {
    "coli": "Escherichia coli",
    "aureus": "Staphylococcus aureus",
    "subtilis": "Bacillus subtilis",
    "pyogenes": "Streptococcus pyogenes",
    "typhimurium": "Salmonella typhimurium",
    "pneumoniae": "Klebsiella pneumoniae",
    "aeruginosa": "Pseudomonas aeruginosa"
}

# --- Interfejs ---
st.title("🧬 BioAnalyzer PRO: Profesjonalna Analiza Genomowa")
uploaded_file = st.sidebar.file_uploader("Wgraj genom (FASTA)", type=["fasta", "fa"])

if uploaded_file:
    with open("input.fasta", "wb") as f:
        f.write(uploaded_file.getbuffer())

    tab1, tab2, tab3 = st.tabs(["1. Wykrywanie ORF", "2. Adnotacja i Funkcje", "3. Predykcja Hosta"])

    # ZAKŁADKA 1
    with tab1:
        if st.button("Uruchom wykrywanie ORF"):
            with st.spinner("Prodigal analizuje..."):
                exe_path = os.path.join(os.getcwd(), "prodigal.exe")
                subprocess.run([exe_path, "-i", "input.fasta", "-a", "proteins.faa", "-o", "genes.gff", "-q"],
                               check=True)
                records = list(SeqIO.parse("proteins.faa", "fasta"))
                data = [{"ID": r.id, "Start": r.description.split(" # ")[1], "Stop": r.description.split(" # ")[2],
                         "Długość (aa)": len(r.seq)} for r in records]
                st.session_state["orf_df"] = pd.DataFrame(data)
        if "orf_df" in st.session_state: st.dataframe(st.session_state["orf_df"], use_container_width=True)

    # ZAKŁADKA 2
    with tab2:
        if st.button("Uruchom Adnotację"):
            with st.spinner("DIAMOND przeszukuje Swiss-Prot..."):
                diamond_path = os.path.join(os.getcwd(), "diamond.exe")
                # Komenda DIAMOND z pobraniem 'stitle' (pełny opis gatunku)
                cmd = [diamond_path, "blastp", "-d", "sprot_db.dmnd", "-q", "proteins.faa", "-o", "annot.tsv",
                       "-k", "1", "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch",
                       "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle"]
                subprocess.run(cmd, check=True)

                # Wczytanie wyników
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

            # Analiza braków
            annotated_ids = st.session_state["annot_df"]["query"].unique()
            missed = st.session_state["orf_df"][~st.session_state["orf_df"]["ID"].isin(annotated_ids)]
            st.warning(f"Liczba białek bez dopasowania (ORFans): {len(missed)}")
            with st.expander("Zobacz listę nieprzypisanych ORF"):
                st.write(missed)

    # ZAKŁADKA 3
    with tab3:
        st.subheader("Predykcja Hosta")
        if "annot_df" in st.session_state:
            scores = {}
            # Analiza gatunków z kolumny stitle
            for title in st.session_state["annot_df"]["stitle"]:
                for key, full_name in TAXA_MAP.items():
                    if key in str(title).lower():
                        scores[full_name] = scores.get(full_name, 0) + 1

            if scores:
                total = sum(scores.values())
                cols = st.columns(len(scores))
                for i, (host, pts) in enumerate(scores.items()):
                    cols[i].metric(label=host, value=f"{(pts / total) * 100:.1f}%")
            else:
                st.error("Brak rozpoznawalnych nazw gospodarzy w wynikach adnotacji.")
        else:
            st.info("Najpierw wykonaj adnotację w zakładce nr 2.")

else:
    st.info("Wgraj plik FASTA w panelu bocznym.")