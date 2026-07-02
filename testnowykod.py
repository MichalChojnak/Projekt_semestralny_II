import streamlit as st
import subprocess
import pandas as pd
import os
import io
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from dna_features_viewer import GraphicFeature, GraphicRecord

# Konfiguracja okna Streamlit
st.set_page_config(page_title="Phage Hybrid Annotator", layout="wide")
st.title("🧬 Phage Annotator PRO — Pipeline Hybrydowy (MMseqs2 + HMMER)")
st.write(
    "Profesjonalne środowisko anotacji genomów fagowych łączące wyszukiwanie sekwencyjne oraz profilowe modele Markowa (HMM).")

# --- PANEL BOCZNY: ŚCIEŻKI I BAZY DANYCH ---
st.sidebar.header("1. Konfiguracja Środowiska")
uploaded_file = st.sidebar.file_uploader("Wgraj genom faga (FASTA)", type=['fasta', 'fas', 'fna', 'fa'])

st.sidebar.subheader("Bazy Referencyjne")
mmseqs_db = st.sidebar.text_input("Ścieżka do bazy MMseqs2", value="bazy/phrogs_mmseqs")
hmmer_db = st.sidebar.text_input("Ścieżka do bazy HMMER (.hmm)", value="bazy/phrogs_hmmer.hmm")

# Progi odcięcia (Cut-offs) dla e-value
st.sidebar.subheader("Progi czułości (E-value)")
mmseqs_evalue = st.sidebar.number_input("MMseqs2 Max E-value", value=1e-5, format="%e")
hmmer_evalue = st.sidebar.number_input("HMMER Max E-value", value=1e-3, format="%e")

# Katalog roboczy
WORKDIR = "phage_hybrid_workdir"
if not os.path.exists(WORKDIR):
    os.makedirs(WORKDIR)

FASTA_PATH = os.path.join(WORKDIR, "genome.fasta")
GENES_TSV = os.path.join(WORKDIR, "phanotate_genes.tsv")
PROTEINS_FAA = os.path.join(WORKDIR, "predicted_proteins.faa")
MMSEQS_OUT = os.path.join(WORKDIR, "mmseqs_results.tsv")
HMMER_OUT = os.path.join(WORKDIR, "hmmer_results.txt")
TMP_DIR = os.path.join(WORKDIR, "tmp")


# --- FUNKCJE BIOINFORMATYCZNE ---

def translate_phanotate_to_proteins(genome_fasta, phanotate_tsv, output_faa):
    """Translacja współrzędnych PHANOTATE na sekwencje aminokwasowe (Biopython)"""
    record = list(SeqIO.parse(genome_fasta, "fasta"))[0]
    genome_seq = record.seq
    df = pd.read_csv(phanotate_tsv, sep=r'\s+', comment='#', header=None)

    proteins = []
    for idx, row in df.iterrows():
        start = int(row[0]) - 1
        end = int(row[1])
        strand = str(row[2])

        gene_seq = genome_seq[start:end]
        if strand in ['-', '-1']:
            gene_seq = gene_seq.reverse_complement()

        protein_seq = gene_seq.translate(table=11, to_stop=True)
        orf_id = f"ORF_{idx + 1}"
        prot_record = SeqRecord(protein_seq, id=orf_id, description=f"coords:{row[0]}-{row[1]}({strand})")
        proteins.append(prot_record)

    SeqIO.write(proteins, output_faa, "fasta")
    return len(proteins)


def parse_hmmer_tblout(tblout_path):
    """Parser dla specyficznego formatu wyjściowego --tblout z programu HMMER"""
    hmmer_data = []
    if not os.path.exists(tblout_path):
        return pd.DataFrame()

    with open(tblout_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            # Format tblout rozdzielany jest spacjami, opis jest na samym końcu (od 22 kolumny)
            parts = line.split(None, 22)
            if len(parts) >= 23:
                target_name = parts[0]
                query_name = parts[2]
                evalue = float(parts[4])
                description = parts[22].strip()
                hmmer_data.append({
                    "ORF_ID": query_name,
                    "HMM_TARGET": target_name,
                    "HMM_EVALUE": evalue,
                    "HMM_DESC": description
                })
    return pd.DataFrame(hmmer_data)


# --- GŁÓWNY POTOK INTERFEJSU ---
if uploaded_file is not None:
    with open(FASTA_PATH, "wb") as f:
        f.write(uploaded_file.getbuffer())
    st.sidebar.success("Genom załadowany poprawnie.")

    tab1, tab2, tab3, tab4 = st.tabs([
        "1. Gene Calling (PHANOTATE)",
        "2. Hybrydowa Anotacja (MMseqs2 + HMMER)",
        "3. Zakres Gospodarza",
        "4. Interaktywna Mapa Genomu"
    ])

    # ==========================================
    # KROK 1: PHANOTATE
    # ==========================================
    with tab1:
        st.header("Identyfikacja Ram Odczytu (ORF)")
        if st.button("Uruchom PHANOTATE", key="btn_phanotate"):
            with st.spinner("PHANOTATE mapuje genom faga..."):
                cmd = f"phanotate {FASTA_PATH} > {GENES_TSV}"
                res = subprocess.run(cmd, shell=True, capture_output=True, text=True)

                if res.returncode == 0 and os.path.exists(GENES_TSV) and os.path.getsize(GENES_TSV) > 0:
                    num_prot = translate_phanotate_to_proteins(FASTA_PATH, GENES_TSV, PROTEINS_FAA)
                    st.success(f"Sukces! Wykryto {num_prot} genów i przetłumaczono je na sekwencje białkowe (.faa).")

                    df_genes = pd.read_csv(GENES_TSV, sep=r'\s+', comment='#', header=None)
                    df_genes.columns = ["START", "STOP", "STRAND", "FRAME", "SCORE"]
                    df_genes.insert(0, "ORF_ID", [f"ORF_{i + 1}" for i in range(len(df_genes))])
                    st.dataframe(df_genes, use_container_width=True)
                else:
                    st.error("Błąd podczas uruchamiania PHANOTATE. Sprawdź logi systemowe.")
                    st.code(res.stderr)

    # ==========================================
    # KROK 2: HYBRYDOWA ANOTACJA (MMseqs2 -> HMMER)
    # ==========================================
    with tab2:
        st.header("Anotacja dwuetapowa (Szybkie dopasowanie + Analiza profili HMM)")
        st.info(
            "Algorytm najpierw mapuje białka przez MMseqs2. Jeśli funkcja jest nieznana (Hypothetical), HMMER analizuje strukturę domeny pod kątem odległej homologii.")

        if not os.path.exists(PROTEINS_FAA):
            st.warning("Uruchom najpierw Krok 1, aby wygenerować sekwencje białkowe.")
        else:
            if st.button("Uruchom Analizę Hybrydową", key="btn_annotation"):
                # --- FAZA A: MMseqs2 ---
                with st.spinner("Faza A: Przeszukiwanie bazy sekwencji (MMseqs2)..."):
                    cmd_mmseqs = f"mmseqs easy-search {PROTEINS_FAA} {mmseqs_db} {MMSEQS_OUT} {TMP_DIR} --format-output 'query,target,pident,evalue,theader' -e {mmseqs_evalue} --v 0"
                    res_mm = subprocess.run(cmd_mmseqs, shell=True, capture_output=True, text=True)

                    if res_mm.returncode != 0:
                        st.error("Błąd krytyczny MMseqs2. Sprawdź ścieżkę bazy danych.")
                        st.code(res_mm.stderr)
                        st.stop()

                # --- FAZA B: HMMER ---
                with st.spinner("Faza B: Skanowanie profilowe 'ciemnej materii' (HMMER hmmscan)..."):
                    cmd_hmmer = f"hmmscan --tblout {HMMER_OUT} -E {hmmer_evalue} --cpu 2 {hmmer_db} {PROTEINS_FAA}"
                    res_hm = subprocess.run(cmd_hmmer, shell=True, capture_output=True, text=True)

                    if res_hm.returncode != 0:
                        st.error("Błąd krytyczny HMMER. Upewnij się, że baza została przygotowana komendą hmmpress.")
                        st.code(res_hm.stderr)
                        st.stop()

                st.success("Integracja danych zakończona pomyślnie!")

                # --- INTEGRACJA WYNIKÓW DO MASTER TABELI ---
                df_m_genes = pd.read_csv(GENES_TSV, sep=r'\s+', comment='#', header=None)
                df_m_genes.columns = ["START", "STOP", "STRAND", "FRAME", "SCORE"]
                df_m_genes.insert(0, "ORF_ID", [f"ORF_{i + 1}" for i in range(len(df_m_genes))])

                # Wczytanie MMseqs2
                if os.path.exists(MMSEQS_OUT) and os.path.getsize(MMSEQS_OUT) > 0:
                    df_mm = pd.read_csv(MMSEQS_OUT, sep='\t', header=None,
                                        names=["ORF_ID", "MM_TARGET", "MM_PIDENT", "MM_EVALUE", "MM_TITLE"])
                    # Zatrzymujemy tylko najlepszy hit dla każdego ORF
                    df_mm = df_mm.drop_duplicates(subset=["ORF_ID"], keep="first")
                else:
                    df_mm = pd.DataFrame(columns=["ORF_ID", "MM_TITLE", "MM_EVALUE"])

                # Wczytanie HMMER
                df_hm = parse_hmmer_tblout(HMMER_OUT)
                if not df_hm.empty:
                    df_hm = df_hm.drop_duplicates(subset=["ORF_ID"], keep="first")

                # Łączenie w jedną ramkę danych
                df_master = pd.merge(df_m_genes, df_mm, on="ORF_ID", how="left")
                df_master = pd.merge(df_master, df_hm, on="ORF_ID", how="left")

                # Logika hybrydowa przypisywania funkcji ostatecznej
                final_annotations = []
                for _, row in df_master.iterrows():
                    mm_title = str(row.get("MM_TITLE", "NaN"))
                    hm_desc = str(row.get("HMM_DESC", "NaN"))

                    # Definiujemy, co uznajemy za brak jasnej funkcji
                    is_mm_hypo = mm_title == "NaN" or "hypothetical" in mm_title.lower() or "unknown" in mm_title.lower()
                    is_hm_valid = hm_desc != "NaN" and "hypothetical" not in hm_desc.lower()

                    if is_mm_hypo and is_hm_valid:
                        # HMMER uratował gen - brak trafienia w MMseqs2, ale jest trafienie domenowe HMM
                        final_annotations.append(f"[HMMER] {hm_desc}")
                    elif mm_title != "NaN":
                        final_annotations.append(f"[MMseqs2] {mm_title}")
                    else:
                        final_annotations.append("Hypothetical protein")

                df_master["OSTATECZNA_FUNKCJA"] = final_annotations
                df_master.to_csv(os.path.join(WORKDIR, "master_annotations.csv"), index=False)

                st.dataframe(df_master[["ORF_ID", "START", "STOP", "OSTATECZNA_FUNKCJA", "MM_EVALUE", "HMM_EVALUE"]],
                             use_container_width=True)

    # ==========================================
    # KROK 3: PREDIKCJA GOSPODARZA
    # ==========================================
    with tab3:
        st.header("Analiza zakresu gospodarza (Host Range)")
        master_path = os.path.join(WORKDIR, "master_annotations.csv")

        if not os.path.exists(master_path):
            st.warning("Ukończ krok 2 (Anotacja Hybrydowa), aby wygenerować tabele funkcjonalną.")
        else:
            df_master = pd.read_csv(master_path)
            rbp_keywords = ["tail fiber", "receptor binding", "spike", "baseplate", "tail protein",
                            "tail fiber protein"]

            rbp_df = df_master[df_master["OSTATECZNA_FUNKCJA"].str.lower().str.contains('|'.join(rbp_keywords))]

            if not rbp_df.empty:
                st.success("Wykryto molekularne determinanty gospodarza (RBP):")
                st.dataframe(rbp_df[["ORF_ID", "START", "STOP", "OSTATECZNA_FUNKCJA"]], use_container_width=True)

                # Próba wyciągnięcia taksonomii bakterii
                hosts = []
                for f in rbp_df["OSTATECZNA_FUNKCJA"]:
                    if "[" in f and "]" in f:
                        hosts.append(f.split("[")[1].split("]")[0])

                if hosts:
                    st.metric(label="Najbardziej prawdopodobny gospodarz docelowy faga:", value=hosts[0])
                    st.caption("Predykcja oparta na homologii domenowej aparatów infekcyjnych ogonka.")
                else:
                    st.info(
                        "Wykryto strukturalne białka ogonka faga, lecz nagłówki baz danych nie zawierają bezpośredniej nazwy gatunkowej bakterii w nawiasach kwadratowych.")
            else:
                st.warning(
                    "Brak wyraźnych markerów białek wiążących receptory. Genom może należeć do faga bezogonowego lub zawierać nieopisane dotąd motywy strukturalne.")

    # ==========================================
    # KROK 4: WIZUALIZACJA GENOMU
    # ==========================================
    with tab4:
        st.header("Fizyczna mapa adnotacji genomu")
        master_path = os.path.join(WORKDIR, "master_annotations.csv")

        if os.path.exists(master_path):
            df_master = pd.read_csv(master_path)
            max_coord = int(df_master["STOP"].max() + 500)

            zoom_range = st.slider("Wybierz okno widoku genomu (bp)", 0, max_coord, (0, min(max_coord, 20000)))

            if st.button("Generuj wizualizację mapy"):
                features = []
                rbp_keywords = ["tail fiber", "receptor binding", "spike", "baseplate"]

                for _, row in df_master.iterrows():
                    if row["START"] >= zoom_range[0] and row["STOP"] <= zoom_range[1]:
                        func_lower = str(row["OSTATECZNA_FUNKCJA"]).lower()
                        strand_dir = 1 if str(row["STRAND"]) in ["+", "1", "+1"] else -1

                        # Dobór palety kolorów i etykiet
                        if any(kw in func_lower for kw in rbp_keywords):
                            color = "#ff4b4b"  # Czerwony dla maszynerii infekcyjnej (RBP)
                            label = f"RBP ({row['ORF_ID']})"
                        elif "hypothetical" in func_lower:
                            color = "#f0ad4e"  # Pomarańczowy dla hipotetycznych
                            label = row["ORF_ID"]
                        else:
                            color = "#4b96ff"  # Niebieski dla stałych funkcji (np. kapsyd, enzymy)
                            label = str(row["OSTATECZNA_FUNKCJA"]).replace("[MMseqs2]", "").replace("[HMMER]",
                                                                                                    "").strip()[
                                        :25] + "..."

                        features.append(GraphicFeature(start=int(row["START"]), end=int(row["STOP"]), strand=strand_dir,
                                                       color=color, label=label))

                if features:
                    record = GraphicRecord(sequence_length=zoom_range[1], features=features)
                    fig, ax = plt.subplots(1, 1, figsize=(15, 5))
                    record.plot(ax=ax, with_ruler=True, elevate_outline_annotations=True)
                    ax.set_xlim(zoom_range[0], zoom_range[1])

                    st.pyplot(fig)

                    # Eksport do PDF
                    pdf_buf = io.BytesIO()
                    fig.savefig(pdf_buf, format="pdf", bbox_inches='tight')
                    pdf_buf.seek(0)
                    st.download_button("Pobierz mapę (PDF)", pdf_buf, "mapa_hybrydowa.pdf", "application/pdf")
                else:
                    st.warning("Brak zidentyfikowanych genów w tym wycinku genomu.")
        else:
            st.info("Mapa będzie dostępna po wygenerowaniu adnotacji w Kroku 2.")
else:
    st.info("Wgraj genom fagowy w formacie FASTA w lewym panelu bocznym, aby rozpocząć proces analizy.")