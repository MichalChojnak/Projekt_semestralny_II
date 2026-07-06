# host_database.py
import pandas as pd
import os


def load_host_database(db_path="host_database.csv"):
    if not os.path.exists(db_path):
        return None, set()
    try:
        df = pd.read_csv(db_path)
        # Oczyszczanie bazy
        banned_words = {"Phage", "Bacteriophage", "Virus", "Unidentified", "Unknown"}
        if 'host_genus' in df.columns:
            genera = set(df['host_genus'].dropna().str.capitalize().unique())
            genera = {g for g in genera if len(g) > 2 and g not in banned_words}
        else:
            genera = set()

        if 'accession' in df.columns:
            df.set_index('accession', inplace=True)

        return df, genera
    except Exception as e:
        return None, set()


def get_host_by_sseqid(sseqid, db_df):
    # Mapowanie accession -> host
    if db_df is not None and sseqid in db_df.index:
        # Zwracamy rodzaj gospodarza (np. Escherichia)
        return db_df.loc[sseqid, 'host_genus']
    return None