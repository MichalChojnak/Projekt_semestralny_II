from __future__ import annotations
import csv
from pathlib import Path
from Bio import SeqIO  # Wymaga instalacji: pip install biopython


# --- FUNKCJE POMOCNICZE ---

def n_content(seq: str) -> float:
    """Zwraca ułamek występowania nukleotydu 'N' w sekwencji."""
    if not seq:
        return 0.0
    return seq.count('N') / len(seq)


def parse_and_clean_file(file_path: Path) -> list[dict]:
    """
    Rozpoznaje format pliku, parsuje go za pomocą Biopython,
    usuwa przerwy (gapy) i konwertuje na listę słowników.
    """
    records = []
    # Automatyczne wykrywanie formatu na podstawie rozszerzenia
    ext = file_path.suffix.lower()
    if ext in [".fasta", ".fa", ".fna"]:
        fmt = "fasta"
    elif ext in [".fastq", ".fq"]:
        fmt = "fastq"
    elif ext in [".gbk", ".gb"]:
        fmt = "genbank"
    else:
        print(f"Pominięto plik {file_path.name}: nieznany format.")
        return records

    try:
        for record in SeqIO.parse(file_path, fmt):
            # Normalizacja: wielkie litery i usunięcie przerw (gapów)
            clean_seq = str(record.seq).upper().replace("-", "")
            records.append({
                "id": record.id,
                "description": record.description,
                "sequence": clean_seq,
                "source_file": file_path.name
            })
    except Exception as e:
        print(f"Błąd podczas czytania pliku {file_path.name}: {e}")

    return records


# --- GŁÓWNA LOGIKA FILTROWANIA (Z TWOICH ZAJĘĆ, ROZBUDOWANA) ---

def filter_sequences(records: list[dict], min_len: int, max_n_pct: float) -> tuple[list[dict], list[dict]]:
    """
    Filtruje sekwencje według parametrów jakości oraz poprawności alfabetu.
    Zwraca: (accepted, rejected)
    """
    accepted, rejected = [], []
    allowed_chars = set("ACGTN")

    for r in records:
        seq = r["sequence"]
        reject = False
        reason = ""

        # 1. Sprawdzenie niedozwolonych znaków (tylko A, T, G, C, N są dozwolone)
        invalid_chars = set(seq) - allowed_chars

        if invalid_chars:
            reject = True
            reason = f"INVALID_ALPHABET (znaleziono: {', '.join(invalid_chars)})"
        # 2. Sprawdzenie długości minimalnej
        elif len(seq) < min_len:
            reject = True
            reason = f"TOO_SHORT (len={len(seq)}, min={min_len})"
        # 3. Sprawdzenie zawartości N
        elif n_content(seq) > max_n_pct:
            reject = True
            reason = f"HIGH_N (n_pct={n_content(seq):.0%}, max={max_n_pct:.0%})"

        if reject:
            rejected.append({**r, "reason": reason})
        else:
            accepted.append(r)

    return accepted, rejected


# --- ZAPIS WYNIKÓW ---

def save_to_fasta(records: list[dict], out_path: Path) -> None:
    """Zapisuje przefiltrowane słowniki do zunifikowanego pliku FASTA."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        for r in records:
            f.write(f">{r['id']} {r['description']} | source:{r['source_file']}\n")
            # Łamanie długich sekwencji na linie po 80 znaków (standard FASTA)
            seq = r["sequence"]
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")


def save_param_compare(results_by_variant: dict, out_path: Path) -> None:
    """Zapisuje tabelę porównawczą wyników dla obu wariantów z zajęć."""
    if not results_by_variant:
        return

    rows = []
    for variant, data in results_by_variant.items():
        rows.append({
            "variant": variant,
            "min_len": data["params"]["min_len"],
            "max_n_pct": data["params"]["max_n_pct"],
            "total_input": data["total"],
            "accepted": data["accepted"],
            "rejected": data["rejected"],
            "accepted_pct": round(data["accepted"] / data["total"] * 100, 1) if data["total"] else 0,
        })
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


# --- URUCHOMIENIE PIPELINE'U ---

def run_pipeline(input_dir: str, output_dir: str, min_len: int = 5000, max_n_pct: float = 0.05):
    """Główna funkcja orkiestrująca działanie programu."""
    in_path = Path(input_dir)
    out_path = Path(output_dir)

    all_records = []

    print(f"Skanowanie katalogu: {in_path}...")
    for file_path in in_path.iterdir():
        if file_path.is_file():
            # Parsowanie i unifikacja do pamięci
            all_records.extend(parse_and_clean_file(file_path))

    print(f"Wczytano {len(all_records)} sekwencji. Rozpoczęcie filtracji...")

    # Filtrowanie
    accepted, rejected = filter_sequences(all_records, min_len, max_n_pct)

    # Zapis zaakceptowanych do jednego pliku FASTA
    final_fasta = out_path / "unified_valid_phages.fasta"
    save_to_fasta(accepted, final_fasta)
    print(f"Zapisano {len(accepted)} poprawnych sekwencji do {final_fasta}")

    # Generowanie raportu z Twojego fragmentu
    results_dict = {
        "variant_1": {
            "params": {"min_len": min_len, "max_n_pct": max_n_pct},
            "total": len(all_records),
            "accepted": len(accepted),
            "rejected": len(rejected)
        }
    }
    report_file = out_path / "qc_report.csv"
    save_param_compare(results_dict, report_file)
    print(f"Zapisano raport odrzutów do {report_file}")


# Jeśli uruchamiamy plik bezpośrednio:
if __name__ == "__main__":
    # Zakładając, że pliki wejściowe masz w folderze 'raw_data', a wyniki chcesz w 'processed_data'
    run_pipeline(input_dir="raw_data", output_dir="processed_data", min_len=5000, max_n_pct=0.05)