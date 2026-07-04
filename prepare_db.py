import pandas as pd
import os

# 1. Wczytaj dane z sequences.csv
if not os.path.exists("sequences.csv"):
    print("BŁĄD: Nie znaleziono pliku sequences.csv")
else:
    # low_memory=False naprawia ostrzeżenie DtypeWarning
    df_refseq = pd.read_csv("sequences.csv", low_memory=False)
    print("Wczytano sequences.csv")

# 2. Wczytaj dane z PHROG
if not os.path.exists("phrog_annotations.csv"):
    print("BŁĄD: Nie znaleziono pliku phrog_annotations.csv! Pobierz go i umieść w folderze projektu.")
else:
    df_phrog = pd.read_csv("phrog_annotations.csv", low_memory=False)
    print("Wczytano phrog_annotations.csv")

    # 3. Połącz (Etap 3: scalanie danych)
    # Tu dopisz resztę swojego kodu, np. merge...
    # master_db = pd.merge(df_refseq, df_phrog, on='Accession', how='left')
    # ...