import streamlit as st
import pandas as pd
import plotly.graph_objects as go
import math
import random

# --- KONFIGURACJA ---
st.set_page_config(page_title="BioAnalyzer PRO v6.2", layout="wide")

KNOWN_GENERA = ["Escherichia", "Salmonella", "Shigella", "Pseudomonas", "Staphylococcus", "Bacillus"]


# --- SILNIK ANALITYCZNY ---

def get_category(opis):
    opis = str(opis).lower()
    if any(x in opis for x in ["rbp", "receptor binding"]): return "Receptor Binding Protein (RBP)"
    if any(x in opis for x in ["tail fiber", "tailfibre"]): return "Tail Fiber"
    if any(x in opis for x in ["tail spike"]): return "Tail Spike"
    if any(x in opis for x in ["endolysin", "lysin"]): return "Endolizyny"
    if any(x in opis for x in ["holin"]): return "Holiny"
    if any(x in opis for x in ["depolymerase"]): return "Depolimerazy"
    if any(x in opis for x in ["baseplate"]): return "Baseplate"
    if any(x in opis for x in ["capsid", "major capsid"]): return "Kapsyd"
    if any(x in opis for x in ["portal"]): return "Białka portalu"
    if any(x in opis for x in ["polymerase"]): return "Polimerazy"
    if any(x in opis for x in ["terminase"]): return "Terminazy"
    if any(x in opis for x in ["integrase"]): return "Integrazy"
    return "Pozostałe"


def mock_detect_orfs(fasta_content):
    """Symuluje wykrywanie ORF-ów w pliku FASTA."""
    lines = fasta_content.decode("utf-8").splitlines()
    sequence_len = len("".join([l for l in lines if not l.startswith(">")]))

    # Symulacja ORF-ów (co 2000 bp)
    orfs = []
    for i in range(0, sequence_len - 500, 2000):
        orfs.append({
            "ORF": f"ORF_{i // 2000 + 1}",
            "Start": i,
            "Stop": i + 1500,
            "Strand": random.choice(["+", "-"]),
            "Length": 1500
        })
    return pd.DataFrame(orfs)


def plot_genomic_map(df):
    fig = go.Figure()
    for _, row in df.iterrows():
        y = 1 if row["Strand"] == "+" else -1
        fig.add_trace(go.Scatter(
            x=[row["Start"], row["Stop"]], y=[y, y],
            mode="lines+markers", name=row["ORF"],
            line=dict(width=10, color="blue" if y == 1 else "red")
        ))
    fig.update_layout(title="Mapa Genetyczna (ORF)", yaxis=dict(tickvals=[-1, 1], ticktext=["Minus", "Plus"]))
    return fig


def calculate_evidence_v5(df, known_genera, weights):
    host_stats = {
        genus: {"score": 0.0, "protein_hits": 0, "categories": set(), "best_identity": 0.0, "best_bitscore": 0.0} for
        genus in known_genera}
    related_genera = {"escherichia": ["salmonella", "shigella"], "salmonella": ["escherichia", "shigella"],
                      "shigella": ["escherichia", "salmonella"]}

    for _, row in df.iterrows():
        title = str(row["stitle"]).lower()
        identity, bitscore = float(row["pident"]), float(row["bitscore"])
        kategoria = row["Kategoria"]
        if identity < 35 or bitscore < 50: continue
        waga = weights.get(kategoria, 0.1)
        identity_score, bitscore_score = min(identity / 100, 1), min(bitscore / 400, 1)
        quality_bonus = 1.25 if identity >= 95 else (1.15 if identity >= 85 else (1.05 if identity >= 70 else 1.0))
        score = identity_score * bitscore_score * waga * quality_bonus
        for genus in known_genera:
            if genus.lower() in title or f"{genus.lower()} phage" in title:
                host_stats[genus]["score"] += score
                host_stats[genus]["protein_hits"] += 1
                host_stats[genus]["categories"].add(kategoria)
                host_stats[genus]["best_identity"] = max(host_stats[genus]["best_identity"], identity)
                host_stats[genus]["best_bitscore"] = max(host_stats[genus]["best_bitscore"], bitscore)
                if genus.lower() in related_genera:
                    for cousin in related_genera[genus.lower()]:
                        if cousin.capitalize() in host_stats:
                            host_stats[cousin.capitalize()]["score"] += (score * 0.4)
                            host_stats[cousin.capitalize()]["categories"].add(kategoria)

    results = {}
    for host, data in host_stats.items():
        diversity = len(data["categories"])
        protein_bonus = min(data["protein_hits"] * 0.05, 0.40)
        confidence = (1 + math.log(diversity + 1)) * (1 + protein_bonus) * (
                    0.7 + 0.3 * (data["best_identity"] / 100)) * (0.7 + 0.3 * min(data["best_bitscore"] / 500, 1))
        results[host] = {"final_score": data["score"] * confidence, "data": data}
    return results


# --- INTERFEJS ---

st.title("🧬 BioAnalyzer PRO: Profesjonalna Analiza Genomowa")

with st.sidebar:
    st.header("Opcje i Status")
    uploaded_file = st.file_uploader("Wgraj genom (FASTA)", type=["fasta", "fna"])

tab1, tab2, tab3, tab4 = st.tabs(["1. Wykrywanie ORF", "2. Adnotacja", "3. Predykcja Hosta", "4. Raport PDF"])

with tab1:
    st.subheader("Wykrywanie ORF i Mapa Genetyczna")
    if uploaded_file:
        if st.button("Uruchom wykrywanie ORF"):
            orfs = mock_detect_orfs(uploaded_file.getvalue())
            st.session_state["orf_df"] = orfs
            st.success(f"Wykryto {len(orfs)} ORF-ów.")

        if "orf_df" in st.session_state:
            st.plotly_chart(plot_genomic_map(st.session_state["orf_df"]), use_container_width=True)
            st.dataframe(st.session_state["orf_df"])
    else:
        st.warning("Wgraj plik FASTA.")

with tab2:
    st.subheader("Adnotacja i Funkcje")
    if st.button("Uruchom Adnotację"):
        data = {"query": ["orf1", "orf2", "orf3"],
                "stitle": ["Receptor binding protein Escherichia phage", "Endolysin Salmonella phage",
                           "Holin Escherichia phage"], "pident": [98.5, 92.0, 95.0], "bitscore": [450, 300, 120]}
        df = pd.DataFrame(data)
        df["Kategoria"] = df["stitle"].apply(get_category)
        st.session_state["annot_df"] = df
        st.success("Adnotacja zakończona!")
    if "annot_df" in st.session_state:
        st.dataframe(st.session_state["annot_df"])

with tab3:
    st.subheader("Host Evidence Engine V5")
    if "annot_df" in st.session_state:
        weights = {"Receptor Binding Protein (RBP)": 15, "Tail Fiber": 12, "Tail Spike": 12, "Depolimerazy": 10,
                   "Endolizyny": 8, "Baseplate": 6, "Holiny": 2, "Kapsyd": 3, "Białka portalu": 3, "Polimerazy": 2,
                   "Terminazy": 2, "Integrazy": 2, "Pozostałe": 0.5, "Brak adnotacji (ORFan)": 0.1}
        results = calculate_evidence_v5(st.session_state["annot_df"], KNOWN_GENERA, weights)
        total_score = sum(res["final_score"] for res in results.values())
        sorted_keys = sorted(results.keys(), key=lambda k: results[k]["final_score"], reverse=True)
        for host in sorted_keys:
            res, data = results[host], results[host]["data"]
            prob = (res["final_score"] / total_score * 100) if total_score > 0 else 0
            c1, c2, c3 = st.columns([2, 5, 1])
            with c1: st.write(f"**{host}**")
            with c2: st.progress(prob / 100)
            with c3: st.write(f"{prob:.1f}%")
            with st.expander("Szczegóły"):
                st.write(f"Liczba trafień: {data['protein_hits']} | Best Id: {data['best_identity']:.1f}%")
    else:
        st.info("Najpierw wykonaj adnotację.")

with tab4:
    st.subheader("Raport PDF")
    st.write("Funkcjonalność w budowie.")