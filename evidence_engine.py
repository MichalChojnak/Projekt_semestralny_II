# evidence_engine.py
import pandas as pd
from protein_classifier import PROTEIN_WEIGHTS
from host_database import get_host_by_sseqid


def calculate_genome_features(orf_df):
    # Punkt 14: Genome Features
    length = orf_df["Stop"].max() if not orf_df.empty else 0
    total_coding = (orf_df["Stop"] - orf_df["Start"]).sum()
    coding_density = (total_coding / length) * 100 if length > 0 else 0
    avg_orf = orf_df["Długość (aa)"].mean() if not orf_df.empty else 0

    return {
        "genome_length": length,
        "coding_density": coding_density,
        "average_orf_len": avg_orf
    }


def predict_host(annot_df, orf_df, host_db_df, known_genera):
    # Punkt 7: Obliczanie coverage. Potrzebujemy długości zapytania z orf_df.
    annot_df = annot_df.merge(orf_df[["ID", "Długość (aa)"]], left_on="query", right_on="ID", how="left")
    annot_df["coverage"] = ((annot_df["qend"] - annot_df["qstart"] + 1) / annot_df["Długość (aa)"]) * 100

    # Słownik do trzymania dowodów XAI (Punkt 15)
    host_evidence = {}

    # Sortowanie po bitscore, aby najlepsze trafienia przetwarzać pierwsze
    annot_df = annot_df.sort_values(by="bitscore", ascending=False)

    # Punkt 8: Śledzenie liczby wystąpień danego białka u konkretnego hosta
    protein_counts_per_host = {}

    for _, row in annot_df.iterrows():
        # Punkt 6: Filtr E-value i Coverage
        if float(row["evalue"]) > 1e-3 or row["coverage"] < 40:
            continue

        kategoria = row["Kategoria"]
        waga_bazowa = PROTEIN_WEIGHTS.get(kategoria, 0.1)

        # Punkt 4 i 5: Oznaczanie hosta primarnie z sseqid, awaryjnie z tytułu
        host = get_host_by_sseqid(row["subject"], host_db_df)
        if not host:
            # Fallback na stitle
            for g in known_genera:
                if g.lower() in str(row["stitle"]).lower():
                    host = g
                    break

        if host:
            host_cap = host.capitalize()

            if host_cap not in host_evidence:
                host_evidence[host_cap] = {
                    "score": 0.0,
                    "proteins_found": {},
                    "details": []
                }
                protein_counts_per_host[host_cap] = {}

            # Punkt 8: Malejące znaczenie (Decay function)
            count = protein_counts_per_host[host_cap].get(kategoria, 0)
            decay_factor = 0.5 ** count  # 1sze trafienie 100%, 2gie 50%, 3cie 25%
            protein_counts_per_host[host_cap][kategoria] = count + 1

            score = (row["bitscore"] / 500) * (row["pident"] / 100) * waga_bazowa * decay_factor

            host_evidence[host_cap]["score"] += score
            host_evidence[host_cap]["proteins_found"][kategoria] = host_evidence[host_cap]["proteins_found"].get(
                kategoria, 0) + 1

            # Punkt 15: Zbieranie danych do Explainable AI
            host_evidence[host_cap]["details"].append({
                "protein": kategoria,
                "identity": row["pident"],
                "bitscore": row["bitscore"],
                "coverage": row["coverage"],
                "subject": row["subject"]
            })

    # Formatowanie wyników
    final_results = {}
    for h, data in host_evidence.items():
        if data["score"] > 0:
            final_results[h] = data

    # Zwraca posortowaną listę krotek: (host, dane_xai)
    return sorted(final_results.items(), key=lambda x: x[1]["score"], reverse=True)