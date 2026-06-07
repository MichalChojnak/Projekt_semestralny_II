import streamlit as st
from Bio import SeqIO
import pandas as pd
import io


# --- FUNKCJE POMOCNICZE I OBLICZENIA BIOLOGICZNE ---

def n_content(seq: str) -> float:
    return seq.count('N') / len(seq) if len(seq) > 0 else 0


def gc_content(seq: str) -> float:
    seq_upper = seq.upper()
    g = seq_upper.count('G')
    c = seq_upper.count('C')
    return ((g + c) / len(seq_upper) * 100) if len(seq_upper) > 0 else 0


def process_uploaded_file(uploaded_file) -> list[dict]:
    """Wczytuje plik z przeglądarki (w pamięci) i parsuje go przez Biopython."""
    records = []
    # Streamlit zwraca plik jako bajty, musimy go zdekodować do tekstu dla Biopythona
    stringio = io.StringIO(uploaded_file.getvalue().decode("utf-8"))

    # Próba zgadnięcia formatu (fasta lub fastq)
    fmt = "fasta" if uploaded_file.name.lower().endswith((".fasta", ".fa", ".fna")) else "fastq"

    for record in SeqIO.parse(stringio, fmt):
        clean_seq = str(record.seq).upper().replace("-", "")
        records.append({
            "id": record.id,
            "length": len(clean_seq),
            "gc_pct": gc_content(clean_seq),
            "n_pct": n_content(clean_seq),
            "sequence": clean_seq,
            "source": uploaded_file.name
        })
    return records


# --- INTERFEJS UŻYTKOWNIKA (GUI) ---

# 1. Ustawienia strony
st.set_page_config(page_title="Phage Host Predictor - QC", page_icon="🧬", layout="wide")
st.title("🧬 Moduł Walidacji Danych Sekwencyjnych")
st.markdown("Wgraj pliki FASTA/FASTQ, aby sprawdzić ich jakość przed uruchomieniem predykcji gospodarza.")

# 2. Pasek boczny z parametrami filtrowania
st.sidebar.header("Parametry Filtrowania")
min_len = st.sidebar.number_input("Minimalna długość sekwencji (bp)", min_value=1, value=5000, step=100)
max_n_pct = st.sidebar.slider("Maksymalna zawartość 'N' (%)", min_value=0.0, max_value=100.0, value=5.0) / 100.0

# 3. Pole do wgrywania plików
uploaded_files = st.file_uploader("Przeciągnij i upuść pliki sekwencji", accept_multiple_files=True,
                                  type=['fasta', 'fa', 'fna', 'fastq', 'fq'])

if uploaded_files:
    st.info(f"Wczytano {len(uploaded_files)} plik(ów). Rozpoczynam analizę...")

    all_records = []
    for f in uploaded_files:
        all_records.extend(process_uploaded_file(f))

    if not all_records:
        st.error("Nie udało się odczytać żadnych sekwencji z podanych plików.")
    else:
        # 4. Filtracja danych
        accepted = []
        rejected = []
        allowed_chars = set("ACGTN")

        for r in all_records:
            invalid_chars = set(r["sequence"]) - allowed_chars
            if invalid_chars:
                r["reject_reason"] = "Niedozwolone znaki"
                rejected.append(r)
            elif r["length"] < min_len:
                r["reject_reason"] = "Zbyt krótka"
                rejected.append(r)
            elif r["n_pct"] > max_n_pct:
                r["reject_reason"] = "Zbyt dużo przerw (N)"
                rejected.append(r)
            else:
                accepted.append(r)

        # 5. GENEROWANIE RAPORTU / DASHBOARDU
        st.header("📊 Raport Jakości (QC)")

        # Kolumny z głównymi metrykami (KPIs)
        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Wszystkie sekwencje", len(all_records))
        col2.metric("Zaakceptowane", len(accepted))
        col3.metric("Odrzucone", len(rejected))

        if accepted:
            df_accepted = pd.DataFrame(accepted)
            avg_len = df_accepted["length"].mean()
            avg_gc = df_accepted["gc_pct"].mean()

            col4.metric("Średnia długość (zaakcept.)", f"{avg_len:,.0f} bp")

            st.subheader("Szczegóły zaakceptowanych sekwencji")
            # Wyświetlamy ładną tabelę (bez samej sekwencji, żeby nie zamulić przeglądarki)
            st.dataframe(df_accepted[["id", "length", "gc_pct", "n_pct", "source"]].style.format(
                {"gc_pct": "{:.1f}%", "n_pct": "{:.2%}"}), use_container_width=True)

        if rejected:
            with st.expander("Pokaż odrzucone sekwencje i powody"):
                df_rejected = pd.DataFrame(rejected)
                st.dataframe(df_rejected[["id", "length", "reject_reason", "source"]], use_container_width=True)

        # 6. Opcja pobrania połączonego i oczyszczonego pliku FASTA
        if accepted:
            st.success("✅ Walidacja zakończona. Możesz pobrać ujednolicony plik do dalszej analizy.")

            # Generowanie zawartości pliku FASTA w pamięci
            fasta_output = io.StringIO()
            for r in accepted:
                fasta_output.write(f">{r['id']} source:{r['source']} len:{r['length']}\n")
                # Łamanie linii co 80 znaków
                for i in range(0, r["length"], 80):
                    fasta_output.write(r["sequence"][i:i + 80] + "\n")

            st.download_button(
                label="📥 Pobierz oczyszczony plik FASTA",
                data=fasta_output.getvalue(),
                file_name="zwalidowane_fagi.fasta",
                mime="text/plain"
            )