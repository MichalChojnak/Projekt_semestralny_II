import pandas as pd
import numpy as np
import os
import math


class EvidenceEngine:
    def __init__(self, db_path="data/processed/phage_host_database.parquet"):
        if not os.path.exists(db_path):
            raise FileNotFoundError(f"Nie znaleziono bazy danych: {db_path}")

        self.db = pd.read_parquet(db_path)

        # 1. Automatyczne mapowanie nazw
        if 'category' in self.db.columns and 'Kategoria' not in self.db.columns:
            self.db.rename(columns={'category': 'Kategoria'}, inplace=True)

        # 2. Bezpieczna Normalizacja Accession
        self.db["Accession"] = self.db["Accession"].astype(str).fillna("")
        self.db["Accession"] = self.db["Accession"].apply(
            lambda x: x.split("|")[-2].replace(".1", "") if "|" in x else x.replace(".1", "")
        )
        self.db["Genus"] = self.db["Host"].astype(str).str.split().str[0]

        # Statystyki do Specificity Index
        self.cat_genus_counts = self.db.groupby(['Kategoria', 'Genus']).size()
        self.cat_total_counts = self.db.groupby('Kategoria').size()

        self.weights = {
            "RBP": 20, "Tail Fiber": 18, "Tail Spike": 17, "Depolymerase": 18,
            "Baseplate": 10, "Endolysin": 4, "Holin": 2,
            "Capsid": 0, "Portal": 0, "Terminase": 0, "Large terminase": 0, "DNA polymerase": 0
        }
        self.ADSORPTION_MODULE = {"Baseplate", "Tail Fiber", "Tail Spike", "RBP"}

    def _get_specificity(self, category, genus):
        cat_in_genus = self.cat_genus_counts.get((category, genus), 0)
        total_cat = self.cat_total_counts.get(category, 1)
        return cat_in_genus / total_cat

    def predict_host(self, annot_df, orf_df):
        print(f"DEBUG: Rozpoczynam predykcję. Wejście DIAMOND: {len(annot_df)} wierszy.")

        # Normalizacja DIAMOND subject
        annot_df = annot_df.copy()
        annot_df["subject"] = annot_df["subject"].astype(str).apply(
            lambda x: x.split("|")[-2].replace(".1", "") if "|" in x else x.replace(".1", "")
        )

        # Merge
        df = annot_df.merge(self.db, left_on='subject', right_on='Accession', how='left')
        df = df.merge(orf_df[['ID', 'Długość (aa)']], left_on='query', right_on='ID', how='left')

        if 'Kategoria_y' in df.columns: df.rename(columns={'Kategoria_y': 'Kategoria'}, inplace=True)

        print(f"DEBUG: Po merge: {len(df)} wierszy. Czy są dopasowania Host: {df['Host'].notna().sum()}")

        genus_results = {}
        rejected_hits = 0

        for _, row in df.iterrows():
            if pd.isna(row['Host']):
                continue  # Pomiń, jeśli brak dopasowania w bazie

            # Coverage
            prot_len = row.get('Długość (aa)', 0)
            if prot_len == 0: prot_len = 1
            coverage = min((row['qend'] - row['qstart'] + 1) / prot_len, 1.0)

            # FILTRY - Jeśli lista jest pusta, zakomentuj te warunki lub je poluzuj!
            if row['pident'] < 70 or row['bitscore'] < 120 or coverage < 0.45:
                rejected_hits += 1
                continue

            genus = row['Genus']
            cat = row.get('Kategoria', 'Other')

            # Scoring
            ev_w = 1.5 if cat in ["RBP", "Tail Fiber", "Tail Spike", "Depolymerase"] else (
                0.8 if cat in ["Endolysin", "Holin"] else 0.25)
            phrog_conf = 1.3 if "PHROG" in str(row.get('stitle', '')) else 0.4
            spec_idx = self._get_specificity(cat, genus)

            weight = self.weights.get(cat, 1.0)
            id_factor = 1.0 if row['pident'] > 95 else (
                0.95 if row['pident'] > 90 else (0.8 if row['pident'] > 80 else 0.6))
            bit_factor = min(row['bitscore'] / 300.0, 1.0)

            score = weight * id_factor * coverage * bit_factor * ev_w * phrog_conf * spec_idx

            if genus not in genus_results:
                genus_results[genus] = {"score": 0, "categories": set(), "evidence": [], "hits": 0}

            genus_results[genus]["score"] += score
            genus_results[genus]["categories"].add(cat)
            genus_results[genus]["hits"] += 1
            genus_results[genus]["evidence"].append(f"{cat} (id: {row['pident']:.1f}%)")

        print(f"DEBUG: Odrzucono filtrów: {rejected_hits}. Znaleziono kandydatów: {len(genus_results)}")

        final_results = []
        for genus, data in genus_results.items():
            cat_bonus = 1 + math.log(len(data["categories"]) + 1)
            prot_bonus = 1 + math.sqrt(data["hits"]) / 5
            mod_bonus = 1.25 if len(self.ADSORPTION_MODULE.intersection(data["categories"])) >= 3 else 1.0
            total_score = data["score"] * cat_bonus * prot_bonus * mod_bonus

            final_results.append({"Host": genus, "Score": total_score, "Evidence": data["evidence"][:10]})

        return sorted(final_results, key=lambda x: x["Score"], reverse=True)