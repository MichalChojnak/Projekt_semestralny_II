import pandas as pd


def calculate_genome_features(orf_df):
    """Oblicza podstawowe statystyki genomu potrzebne do raportu."""
    length = orf_df["Stop"].max() if not orf_df.empty else 0
    total_coding = (orf_df["Stop"] - orf_df["Start"]).sum()
    coding_density = (total_coding / length) * 100 if length > 0 else 0
    avg_orf = orf_df["Długość (aa)"].mean() if not orf_df.empty else 0
    return {
        "genome_length": length,
        "coding_density": coding_density,
        "average_orf_len": avg_orf
    }


def predict_host(annot_df, orf_df, tax_df):
    """
    Zaawansowany Evidence Engine:
    1. Oblicza coverage (length/qlen).
    2. Łączy dane z taksonomią (sseqid -> species).
    3. Nadaje wagi: bakterie (1.0) vs fagi (0.5 jako dowód pośredni).
    4. Zwraca listę krotek (gatunek, {'score': X, 'confidence': Y}).
    """

    # 1. Obliczanie Coverage
    annot_df["coverage"] = (annot_df["length"] / annot_df["qlen"]) * 100

    # 2. Filtracja surowa
    hits = annot_df[(annot_df["evalue"] <= 1e-5) & (annot_df["coverage"] >= 30)].copy()

    # 3. Łączenie z bazą taksonomiczną
    hits = hits.merge(tax_df, left_on="subject", right_on="sseqid", how="left")

    # 4. Nadawanie wag (Evidence-Based Phage Homology Transfer)
    def calculate_weight(row):
        species = str(row["species"]).lower()
        # Jeśli trafienie jest wirusem/fagiem, traktujemy to jako dowód pośredni (waga 0.5)
        if any(kw in species for kw in ["virus", "phage", "tequintavirus"]):
            return 0.5
        return 1.0  # Bezpośrednie dopasowanie do bakterii to twardy dowód (waga 1.0)

    hits["weight"] = hits.apply(calculate_weight, axis=1)

    # Obliczanie Score
    hits["score"] = hits["bitscore"] * (hits["pident"] / 100) * (hits["coverage"] / 100) * hits["weight"]

    # 5. Agregacja
    host_summary = hits.groupby("species").agg({
        "score": "sum",
        "query": "nunique"
    }).rename(columns={"query": "protein_count"})

    # Filtracja słabych sygnałów
    host_summary = host_summary[host_summary["score"] > 0.5]

    if host_summary.empty:
        return []

    # Obliczanie pewności
    total_score = host_summary["score"].sum()
    host_summary["confidence"] = (host_summary["score"] / total_score) * 100

    # Sortowanie
    host_summary = host_summary.sort_values(by="score", ascending=False).reset_index()

    # Konwersja na listę krotek (wymagane przez app.py i report.py)
    final_results = []
    for _, row in host_summary.iterrows():
        final_results.append((
            row["species"],
            {
                "score": row["score"],
                "confidence": row["confidence"],
                "protein_count": row["protein_count"]
            }
        ))

    return final_results