import pandas as pd
from Bio import SeqIO
import re


def parse_phrogs_fasta(file_path):
    """Parsuje combined_phrogs.fasta (format z !! separatorami)"""
    data = []
    for record in SeqIO.parse(file_path, "fasta"):
        parts = [p.strip() for p in record.description.split("!!")]
        data.append({
            "accession": parts[5] if len(parts) > 5 else None,
            "phrog_id": parts[0].replace("##", "").strip(),
            "function_phrog": parts[4] if len(parts) > 4 else None
        })
    return pd.DataFrame(data)


def parse_ncbi_fasta(file_path):
    """Parsuje ncbi_phages.fasta (format z nawiasami kwadratowymi [])"""
    data = []
    for record in SeqIO.parse(file_path, "fasta"):
        # record.id to np. YP_011633683.1
        # record.description to pełny nagłówek

        # Wyciąganie hosta z nawiasów []
        match = re.search(r'\[(.*?)\]', record.description)
        host = match.group(1) if match else "Unknown"

        data.append({
            "accession": record.id,
            "host": host,
            "description": record.description.split('|')[
                -1].strip() if '|' in record.description else record.description
        })
    return pd.DataFrame(data)


def run_pipeline():
    # 1. Parsowanie plików
    print("Parsowanie plików FASTA...")
    df_phrog = parse_phrogs_fasta("combined_phrogs.fasta")
    df_ncbi = parse_ncbi_fasta("ncbi_phages.fasta")

    # 2. Integracja (Etap 3 z Twojego planu)
    print("Łączenie danych...")
    # Merge po kolumnie 'accession'
    final_df = pd.merge(df_phrog, df_ncbi, on="accession", how="inner")

    # 3. Zapis wyniku
    output_file = "final_analysis_report.csv"
    final_df.to_csv(output_file, index=False)

    print(f"Sukces! Wynik zapisano do {output_file}")
    print(final_df.head())


if __name__ == "__main__":
    run_pipeline()