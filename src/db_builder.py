import pandas as pd
import os


def build_database():
    print("Rozpoczynam budowę bazy danych...")

    # --- KONFIGURACJA ŚCIEŻEK ---
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    path_ncbi = os.path.join(base_dir, "data", "raw", "sequences.csv")
    path_mapping = os.path.join(base_dir, "data", "raw", "mapping.tsv")
    path_phrog = os.path.join(base_dir, "data", "raw", "phrog_annot_v4.tsv")
    out_dir = os.path.join(base_dir, "data", "processed")
    out_path = os.path.join(out_dir, "phage_host_database.parquet")

    # --- SPRAWDZENIE ---
    files_to_check = [path_ncbi, path_mapping, path_phrog]
    for p in files_to_check:
        if not os.path.exists(p):
            print(f"BŁĄD: Nie znaleziono pliku: {p}")
            return

    # --- WCZYTYWANIE DANYCH ---
    print("Wczytuję pliki...")
    # low_memory=False wyeliminuje ostrzeżenie DtypeWarning
    df_ncbi = pd.read_csv(path_ncbi, low_memory=False)

    # Mapowanie z DIAMOND
    df_mapping = pd.read_csv(path_mapping, sep='\t', names=['Accession', 'phrog'])

    # Adnotacje PHROG
    df_phrog = pd.read_csv(path_phrog, sep='\t')

    # --- KONWERSJA TYPÓW (Naprawa błędu merge) ---
    # Konwertujemy kolumnę 'phrog' na tekst w obu tabelach, żeby typy się zgadzały
    df_mapping['phrog'] = df_mapping['phrog'].astype(str)
    df_phrog['phrog'] = df_phrog['phrog'].astype(str)

    # --- ŁĄCZENIE (MERGE) ---
    print("Łączę dane...")
    # 1. NCBI + Mapowanie
    full_df = df_ncbi.merge(df_mapping, on='Accession', how='left')

    # 2. Wynik + Adnotacje PHROG
    full_df = full_df.merge(df_phrog, on='phrog', how='left')

    # --- ZAPIS ---
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    full_df.to_parquet(out_path, engine="pyarrow")
    print("-" * 30)
    print(f"SUKCES! Baza zapisana w: {out_path}")
    print(f"Liczba wierszy w bazie: {len(full_df)}")


if __name__ == "__main__":
    build_database()