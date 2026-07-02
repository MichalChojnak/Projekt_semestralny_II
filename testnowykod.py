import streamlit as st
import pandas as pd
import os
import subprocess
import re
import io
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from dna_features_viewer import GraphicFeature, GraphicRecord

st.set_page_config(page_title="Local Windows Phage Annotator", layout="wide")
st.title("🧬 Lokalny Pipeline Fagowy — 100% Windows (No Cloud)")
st.write(
    "Wszystkie obliczenia i bazy danych są przetwarzane lokalnie na Twoim komputerze bez wysyłania danych na zewnętrzne serwery.")

# --- KONFIGURACJA ŚCIEŻEK (LOKALNY WINDOWS) ---
st.sidebar.header("1. Ścieżki systemowe Windows")
uploaded_file = st.sidebar.file_uploader("Wgraj genom faga (FASTA)", type=['fasta', 'fa', 'fna'])

# Ścieżki do lokalnych plików .exe na Windowsie
HMMER_BIN = st.sidebar.text_input("Ścieżka do hmmscan.exe", value="bin/hmmscan.exe")
HMM_DB = st.sidebar.text_input("Lokalna baza HMM (np. phrogs.hmm)", value="bazy/phrogs_annotated.hmm")
HMM_EVALUE = st.sidebar.number_input("Próg czułości HMM (E-value)", value=1e-4, format="%e")

WORKDIR = "phage_local_windows_db"
if not os.path.exists(WORKDIR):
    os.makedirs(WORKDIR)

FASTA_PATH = os.path.join(WORKDIR, "genome.fasta")
PROTEINS_FAA = os.path.join(WORKDIR, "predicted_proteins.faa")
HMMER_OUT = os.path.join(WORKDIR, "hmmer_local_results.txt")
FINAL_CSV = os.path.join(WORKDIR, "local_final_annotations.csv")


# --- LOKALNE SILNIKI BIOINFORMATYCZNE ---

def local_biopython_orf_finder(genome_fasta, min_aa_len=60):
    """Lokalny detektor genów (6 ramek odczytu) - działa w 100% na Windows"""
    record = list(SeqIO.parse(genome_fasta, "fasta"))[0]
    sequence = record.seq
    orfs = []
    orf_counter = 1

    for strand, seq in [(1, sequence), (-1, sequence.reverse_complement())]:
        for frame in range(3):
            trans = seq[frame:].translate(table=11)
            trans_len = len(trans)
            aa_start = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len
                if (aa_end - aa_start) >= min_aa_len:
                    if strand == 1:
                        start_pos = frame + aa_start * 3 + 1
                        end_pos = frame + (aa_end + 1) * 3
                    else:
                        gen_len = len(sequence)
                        start_pos = gen_len - (frame + (aa_end + 1) * 3) + 1
                        end_pos = gen_len - (frame + aa_start * 3)

                    orfs.append({
                        "ORF_ID": f"ORF_{orf_counter}",
                        "START": start_pos,
                        "STOP": end_pos,
                        "STRAND": "+" if strand == 1 else "-",
                        "SEQUENCE": str(trans[aa_start:aa_end])
                    })
                    orf_counter += 1
                aa_start = aa_end + 1

    df = pd.DataFrame(orfs)
    records = [SeqRecord(Seq(row["SEQUENCE"]), id=row["ORF_ID"]) for _, row in df.iterrows()]
    SeqIO.write(records, PROTEINS_FAA, "fasta")
    return df


def parse_hmmer_with_host(tblout_path):
    """
    Parser wyciągający funkcję oraz taksonomię gospodarza (Host) z lokalnego pliku HMMER.
    Szuka znaczników typu [Escherichia coli] lub 'host: Staphylococcus' w liniach opisu bazy danych.
    """
    hmmer_data = []
    if not os.path.exists(tblout_path):
        return pd.DataFrame()

    with open(tblout_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split(None, 22)
            if len(parts) >= 23:
                target_name = parts[0]  # Nazwa hitu z bazy (np. PHROG_00123)
                query_name = parts[2]  # Nasz ORF_X
                evalue = float(parts[4])
                description = parts[22].strip()  # Pełny opis zawierający funkcję i gospodarza

                # REGEKS: Szukamy nazwy organizmu w nawiasach kwadratowych [ ]
                host_match = re.search(r'\[(.*?)\]', description)
                host_detected = host_match.group(1) if host_match else "Nieznany (Szeroki zakres)"

                # Czyszczenie opisu z tagów taksonomicznych na potrzeby ładnego wyświetlania funkcji
                clean_desc = re.sub(r'\[(.*?)\]', '', description).strip()

                hmmer_data.append({
                    "ORF_ID": query_name,
                    "HIT_DOMENY": target_name,
                    "EVALUE": evalue,
                    "FUNKCJA": clean_desc,
                    "WYKRYTY_GOSPODARZ": host_detected
                })

    df = pd.DataFrame(hmmer_data)
    if not df.empty:
        # Zostawiamy tylko najlepsze dopasowanie dla każdego genu
        df = df.drop_duplicates(subset=["ORF_ID"], keep="first")
    return df


# --- INTERFEJS STREMLIT ---
if uploaded_file is not None:
    with open(FASTA_PATH, "wb") as f:
        f.write(uploaded_file.getbuffer())
    st.sidebar.success("Genom załadowany lokalnie.")

    tab1, tab2, tab3 = st.tabs(["1. Detekcja ORF", "2. Lokalny HMMER & Analiza Gospodarza", "3. Mapa Genomu"])

    with tab1:
        st.header("Krok 1: Wykrywanie genów strukturalnych faga")
        if st.button("Uruchom Lokalny Gene-Finder"):
            with st.spinner("Przeszukiwanie nici DNA..."):
                df_orfs = local_biopython_orf_finder(FASTA_PATH)
                st.success(f"Zakończono. Wykryto {len(df_orfs)} genów i zapisano sekwencje białkowe.")
                st.dataframe(df_orfs[["ORF_ID", "START", "STOP", "STRAND"]], use_container_width=True)
                df_orfs.to_csv(os.path.join(WORKDIR, "orfs.csv"), index=False)

    with tab2:
        st.header("Krok 2: Skanowanie Profilami HMM i Mapowanie Pochodzenia")
        st.write("Program uruchomi lokalny plik `hmmscan.exe` i przeanalizuje homologiczne domeny białkowe.")

        if not os.path.exists(PROTEINS_FAA):
            st.warning("Najpierw wykonaj Krok 1.")
        else:
            if st.button("Uruchom Lokalny HMMER"):
                with st.spinner("HMMER przeszukuje lokalną bazę danych (to może chwilę potrwać)..."):
                    # Wywołanie lokalnego narzędzia .exe na Windowsie przez subprocess
                    cmd = f'"{HMMER_BIN}" --tblout "{HMMER_OUT}" -E {HMM_EVALUE} "{H_DB}" "{PROTEINS_FAA}"'
                    # Uwaga: dla celów testowych upewnij się, że przekazujesz poprawne ścieżki w cudzysłowach
                    cmd = f'"{HMMER_BIN}" --tblout {HMMER_OUT} -E {HMM_EVALUE} {HMM_DB} {PROTEINS_FAA}'

                    res = subprocess.run(cmd, shell=True, capture_output=True, text=True)

                    if res.returncode != 0:
                        st.error("Błąd podczas uruchamiania lokalnego hmmscan.exe. Sprawdź logi:")
                        st.code(res.stderr)
                    else:
                        st.success("Skanowanie HMMER zakończone pomyślnie!")

                        # Integracja wyników
                        df_orfs = pd.read_csv(os.path.join(WORKDIR, "orfs.csv"))
                        df_hmm = parse_hmmer_with_host(HMMER_OUT)

                        if df_hmm.empty:
                            st.warning("HMMER nie znalazł żadnych dopasowań powyżej progu E-value.")
                            df_master = df_orfs.copy()
                            df_master["FUNKCJA"] = "Hypothetical protein"
                            df_master["WYKRYTY_GOSPODARZ"] = "Brak danych"
                        else:
                            df_master = pd.merge(df_orfs, df_hmm, on="ORF_ID", how="left")
                            df_master["FUNKCJA"] = df_master["FUNKCJA"].fillna("Hypothetical protein")
                            df_master["WYKRYTY_GOSPODARZ"] = df_master["WYKRYTY_GOSPODARZ"].fillna("Brak danych")

                        df_master.to_csv(FINAL_CSV, index=False)
                        st.dataframe(df_master[["ORF_ID", "START", "STOP", "FUNKCJA", "WYKRYTY_GOSPODARZ", "EVALUE"]],
                                     use_container_width=True)

                        # --- SEKCJA PREDYKCJI GOSPODARZA ---
                        st.subheader("🎯 Wynik Analizy Zakresu Gospodarza (Host Range)")
                        # Filtrujemy tylko geny odpowiedzialne za rozpoznawanie komórki (RBP / Tail fibers)
                        rbp_keywords = ["tail", "fiber", "receptor", "binding", "baseplate", "spike"]
                        rbp_hits = df_master[df_master["FUNKCJA"].str.lower().str.contains('|'.join(rbp_keywords))]

                        if not rbp_hits.empty:
                            st.info(
                                "Wykryto białka aparatu infekcyjnego faga. Oto organizmy źródłowe homologicznych domen:")
                            st.dataframe(rbp_hits[["ORF_ID", "FUNKCJA", "WYKRYTY_GOSPODARZ"]], use_container_width=True)

                            # Pobieramy najczęstszego gospodarza z wykrytych markerów RBP
                            real_hosts = rbp_hits[rbp_hits["WYKRYTY_GOSPODARZ"] != "Brak danych"]["WYKRYTY_GOSPODARZ"]
                            if not real_hosts.empty:
                                predicted_host = real_hosts.mode()[0]
                                st.metric(label="⚠️ Prawdopodobny gospodarz docelowy:", value=predicted_host)
                        else:
                            st.warning("Nie wykryto białek ogonka o znanej homologii taksonomicznej.")

    with tab3:
        st.header("Krok 3: Wizualizacja genomu")
        if os.path.exists(FINAL_CSV):
            df_master = pd.read_csv(FINAL_CSV)
            max_c = int(df_master["STOP"].max() + 500)
            zoom = st.slider("Wycinek genomu", 0, max_c, (0, min(max_c, 20000)))

            if st.button("Generuj mapę fizyczną"):
                features = []
                for _, row in df_master.iterrows():
                    if row["START"] >= zoom[0] and row["STOP"] <= zoom[1]:
                        s = 1 if row["STRAND"] == "+" else -1
                        func = str(row["FUNKCJA"]).lower()

                        color = "#ff4b4b" if "tail" in func or "fiber" in func else (
                            "#f0ad4e" if "hypothetical" in func else "#4b96ff")
                        features.append(
                            GraphicFeature(start=int(row["START"]), end=int(row["STOP"]), strand=s, color=color,
                                           label=row["ORF_ID"]))

                if features:
                    record = GraphicRecord(sequence_length=zoom[1], features=features)
                    fig, ax = plt.subplots(figsize=(12, 4))
                    record.plot(ax=ax, with_ruler=True)
                    st.pyplot(fig)
        else:
            st.info("Anotuj genom w Kroku 2, aby odblokować mapowanie.")