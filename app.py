# app.py
import streamlit as st
import pandas as pd
from Bio import SeqIO
import subprocess
import os
import plotly.graph_objects as go

# Importy z naszych nowych modułów
from host_database import load_host_database
from protein_classifier import kategoryzuj_bialko, CATEGORY_COLORS
from evidence_engine import predict_host, calculate_genome_features

st.set_page_config(layout="wide", page_title="BioAnalyzer PRO v7.0 Modular")

# --- Inicjalizacja Danych ---
HOST_DB_DF, KNOWN_GENERA = load_host_database()

with st.sidebar:
    st.title("⚙️ Opcje i Status")
    if HOST_DB_DF is not None:
        st.success(f"✅ Baza aktywna ({len(KNOWN_GENERA)} hostów)")
    else:
        st.warning("⚠️ Tryb ograniczony (Brak host_database.csv)")
    uploaded_file = st.file_uploader("Wgraj genom (FASTA)", type=["fasta", "fa"])

st.title("🧬 BioAnalyzer PRO: Profesjonalna Analiza Genomowa")

if uploaded_file:
    with open("input.fasta", "wb") as f:
        f.write(uploaded_file.getbuffer())

    tab1, tab2, tab3, tab4 = st.tabs([
        "1. Wykrywanie ORF",
        "2. Adnotacja DIAMOND",
        "3. Predykcja Hosta & XAI",
        "4. Raport PDF"
    ])

    # --- ZAKŁADKA 1 ---
    with tab1:
        if st.button("Uruchom Prodigal"):
            with st.spinner("Szukanie ORF..."):
                subprocess.run(["prodigal.exe", "-i", "input.fasta", "-a", "proteins.faa", "-o", "genes.gff", "-q"],
                               check=True)
                records = list(SeqIO.parse("proteins.faa", "fasta"))
                data = [{
                    "ID": r.id,
                    "Start": int(r.description.split(" # ")[1]),
                    "Stop": int(r.description.split(" # ")[2]),
                    "Nić": "+" if r.description.split(" # ")[3] == "1" else "-",
                    "Długość (aa)": len(r.seq)
                } for r in records]
                st.session_state["orf_df"] = pd.DataFrame(data)

        if "orf_df" in st.session_state:
            st.dataframe(st.session_state["orf_df"], use_container_width=True)

    # --- ZAKŁADKA 2 ---
    with tab2:
        if st.button("Uruchom Adnotację DIAMOND"):
            with st.spinner("Szukanie homologów..."):
                cmd = ["diamond.exe", "blastp", "-d", "sprot_db.dmnd", "-q", "proteins.faa", "-o", "annot.tsv",
                       "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                       "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle"]
                subprocess.run(cmd, check=True)

                df = pd.read_csv("annot.tsv", sep="\t", names=[
                    "query", "subject", "pident", "length", "mismatch", "gap",
                    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle"
                ])
                # Punkt 11: Zmiana subject na stitle!
                df["Kategoria"] = df["stitle"].apply(kategoryzuj_bialko)
                st.session_state["annot_df"] = df

        if "annot_df" in st.session_state:
            st.dataframe(st.session_state["annot_df"][["query", "subject", "Kategoria", "evalue", "bitscore"]],
                         use_container_width=True)

    # --- ZAKŁADKA 3 (Punkt 3, 14, 15) ---
    with tab3:
        if st.button("Przewiduj Gospodarza"):
            if "annot_df" in st.session_state and "orf_df" in st.session_state:
                # Genome features
                features = calculate_genome_features(st.session_state["orf_df"])
                st.session_state["genome_features"] = features

                # Wywołanie zewnętrznego silnika (Punkt 3)
                results = predict_host(
                    st.session_state["annot_df"],
                    st.session_state["orf_df"],
                    HOST_DB_DF,
                    KNOWN_GENERA
                )

                # Punkt 10: Zapisujemy wyniki do sesji dla PDF
                st.session_state["host_results"] = results
            else:
                st.error("Uruchom najpierw zakładki 1 i 2.")

        # Wyświetlanie Explainable AI (Punkt 15)
        if "host_results" in st.session_state and st.session_state["host_results"]:
            st.subheader("Wyniki Predykcji (Explainable AI)")
            best_host, data = st.session_state["host_results"][0]

            col1, col2 = st.columns([1, 1])
            with col1:
                st.metric(label="Najbardziej prawdopodobny gospodarz", value=best_host)
                st.markdown("### Dlaczego?")
                for prot, count in data["proteins_found"].items():
                    st.markdown(f"✔ **{count}x** {prot}")

                st.markdown(f"✔ Znaleziono łącznie **{len(data['details'])}** niezależnych homologów")

            with col2:
                # Ekstrakcja z genomu (Punkt 14)
                st.markdown("### Cechy genomu")
                gf = st.session_state["genome_features"]
                st.write(f"- **Długość:** {gf['genome_length']:,} bp")
                st.write(f"- **Zagęszczenie genów:** {gf['coding_density']:.1f}%")
                st.write(f"- **Średnie ORF:** {gf['average_orf_len']:.1f} aa")

    # --- ZAKŁADKA 4 ---
    with tab4:
        st.write("Generowanie raportu PDF opartego na `st.session_state['host_results']`.")
        # Tutaj wywołasz swój moduł report.py, który poprawnie zaczyta nowe wyniki.