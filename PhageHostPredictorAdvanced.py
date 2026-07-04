import streamlit as st
import pandas as pd
from Bio import SeqIO
import subprocess
import os
import plotly.graph_objects as go
from fpdf import FPDF
import math

st.set_page_config(layout="wide", page_title="BioAnalyzer PRO v6.1")


# --- MODUŁ BAZY WIEDZY ---
@st.cache_data
def load_host_database():
    db_path = "host_database.csv"
    if os.path.exists(db_path):
        try:
            df = pd.read_csv(db_path)
            # Pobieramy unikalne nazwy rodzajów bakterii (np. Escherichia, Pseudomonas)
            genera = set(df['host_genus'].dropna().str.capitalize().unique())
            # Oczyszczamy z przypadkowych słów, które mogły się pobrać z NCBI
            banned_words = {"Phage", "Bacteriophage", "Virus", "Unidentified", "Unknown"}
            genera = {g for g in genera if len(g) > 2 and g not in banned_words}
            return df, genera
        except Exception as e:
            return None, set()
    return None, set()


# Ładowanie bazy przy starcie
HOST_DB_DF, KNOWN_GENERA = load_host_database()

# --- Moduł 2: Protein Importance ---
PROTEIN_WEIGHTS = {
    "Receptor Binding Protein (RBP)": 12,
    "Tail Fiber": 10,
    "Tail Spike": 10,
    "Depolimerazy": 9,
    "Adsorpcja": 8,
    "Baseplate": 6,
    "Lizyny i Holiny": 2,
    "Kapsyd": 0,
    "Białka portalu": 0,
    "Polimerazy": 0,
    "Terminazy": 0,
    "Integrazy": 0,
    "Pozostałe": 1
}

CATEGORY_COLORS = {
    "Receptor Binding Protein (RBP)": "#ff7f0e",
    "Tail Fiber": "#ffbb78",
    "Tail Spike": "#ff9896",
    "Depolimerazy": "#8c564b",
    "Adsorpcja": "#e377c2",
    "Baseplate": "#1f77b4",
    "Lizyny i Holiny": "#d62728",
    "Kapsyd": "#2ca02c",
    "Białka portalu": "#17becf",
    "Polimerazy": "#9467bd",
    "Terminazy": "#c5b0d5",
    "Integrazy": "#f7b6d2",
    "Pozostałe": "#7f7f7f",
    "Brak adnotacji (ORFan)": "#c7c7c7"
}


def kategoryzuj_bialko(opis):
    opis = str(opis).lower()
    if any(x in opis for x in ["receptor binding", "rbp"]): return "Receptor Binding Protein (RBP)"
    if "tail fiber" in opis: return "Tail Fiber"
    if "tail spike" in opis: return "Tail Spike"
    if "depolymerase" in opis: return "Depolimerazy"
    if "adsorption" in opis: return "Adsorpcja"
    if "baseplate" in opis: return "Baseplate"
    if any(x in opis for x in ["lysin", "endolysin", "holin"]): return "Lizyny i Holiny"
    if any(x in opis for x in ["capsid", "coat", "head"]): return "Kapsyd"
    if "portal" in opis: return "Białka portalu"
    if "polymerase" in opis: return "Polimerazy"
    if "terminase" in opis: return "Terminazy"
    if "integrase" in opis: return "Integrazy"
    return "Pozostałe"


def usun_pl(tekst):
    zamienniki = {'ą': 'a', 'ć': 'c', 'ę': 'e', 'ł': 'l', 'ń': 'n', 'ó': 'o', 'ś': 's', 'ź': 'z', 'ż': 'z',
                  'Ą': 'A', 'Ć': 'C', 'Ę': 'E', 'Ł': 'L', 'Ń': 'N', 'Ó': 'O', 'Ś': 'S', 'Ź': 'Z', 'Ż': 'Z'}
    for pl, asc in zamienniki.items():
        tekst = tekst.replace(pl, asc)
    return tekst


# --- Interfejs (Pasek boczny) ---
with st.sidebar:
    st.title("⚙️ Opcje i Status")

    # Wyświetlanie statusu nowej bazy
    st.markdown("### Status Bazy Wiedzy")
    if HOST_DB_DF is not None:
        st.success(
            f"✅ Baza aktywna\n- **{len(HOST_DB_DF):,}** białek fagowych\n- **{len(KNOWN_GENERA)}** unikalnych rodzajów gospodarzy")
    else:
        st.error("❌ Brak pliku `host_database.csv`. Engine będzie działał w trybie ograniczonym.")

    st.markdown("---")
    uploaded_file = st.file_uploader("Wgraj genom (FASTA)", type=["fasta", "fa"])

st.title("🧬 BioAnalyzer PRO: Profesjonalna Analiza Genomowa")

if uploaded_file:
    with open("input.fasta", "wb") as f:
        f.write(uploaded_file.getbuffer())

    tab1, tab2, tab3, tab4 = st.tabs(
        ["1. Wykrywanie ORF i Mapa", "2. Adnotacja i Funkcje", "3. Predykcja Hosta (Evidence V2)", "4. Raport PDF"])

    # ZAKŁADKA 1: ORF i WIZUALIZACJA
    with tab1:
        if st.button("Uruchom wykrywanie ORF (ATG, GTG, TTG)"):
            with st.spinner("Prodigal analizuje..."):
                exe_path = os.path.join(os.getcwd(), "prodigal.exe")
                subprocess.run(
                    [exe_path, "-i", "input.fasta", "-a", "proteins.faa", "-o", "genes.gff", "-q", "-g", "11"],
                    check=True)
                records = list(SeqIO.parse("proteins.faa", "fasta"))

                # ZMODYFIKOWANA LINIJKA - dodane pole "Sekwencja (aa)"
                data = [{"ID": r.id, "Start": int(r.description.split(" # ")[1]),
                         "Stop": int(r.description.split(" # ")[2]),
                         "Nić": "+" if r.description.split(" # ")[3] == "1" else "-", "Długość (aa)": len(r.seq),
                         "Sekwencja (aa)": str(r.seq)} for r in records]

                st.session_state["orf_df"] = pd.DataFrame(data)

        if "orf_df" in st.session_state:
            has_annot = "annot_df" in st.session_state
            st.subheader("Interaktywna Mapa Genomu")

            fig = go.Figure()
            max_stop = st.session_state["orf_df"]["Stop"].max()
            fig.add_trace(
                go.Scatter(x=[0, max_stop], y=[0, 0], mode="lines", line=dict(color="black", width=2), hoverinfo="skip",
                           showlegend=False))

            if has_annot:
                for cat, color in CATEGORY_COLORS.items():
                    fig.add_trace(
                        go.Scatter(x=[None], y=[None], mode='markers', marker=dict(size=10, color=color), name=cat))

            arrow_len = max_stop * 0.015

            for _, row in st.session_state["orf_df"].iterrows():
                start, stop, strand, g_id = row["Start"], row["Stop"], row["Nić"], row["ID"]
                y_base = 0.4 if strand == "+" else -0.4
                color = "rgba(31, 119, 180, 0.8)" if strand == "+" else "rgba(214, 39, 40, 0.8)"
                cat_text = ""

                if has_annot:
                    annot_match = st.session_state["annot_df"][st.session_state["annot_df"]["query"] == g_id]
                    if not annot_match.empty:
                        kategoria = annot_match.iloc[0]["Kategoria"]
                        color = CATEGORY_COLORS.get(kategoria, CATEGORY_COLORS["Pozostałe"])
                        cat_text = f"<br>Kategoria: {kategoria}"
                    else:
                        color = CATEGORY_COLORS["Brak adnotacji (ORFan)"]
                        cat_text = "<br>Kategoria: Brak adnotacji (ORFan)"

                al = min(arrow_len, (stop - start) * 0.4)
                x_poly = [start, stop - al, stop, stop - al, start, start] if strand == '+' else [stop, start + al,
                                                                                                  start, start + al,
                                                                                                  stop, stop]
                y_poly = [y_base - 0.25, y_base - 0.25, y_base, y_base + 0.25, y_base + 0.25, y_base - 0.25]

                hover_text = f"<b>{g_id}</b><br>Start: {start}<br>Stop: {stop}<br>Nić: {strand}<br>Długość: {row['Długość (aa)']} aa{cat_text}"
                fig.add_trace(
                    go.Scatter(x=x_poly, y=y_poly, fill="toself", fillcolor=color, line=dict(color="black", width=1),
                               name=g_id, text=hover_text, hoverinfo="text", showlegend=False))

            fig.update_layout(yaxis=dict(showticklabels=False, range=[-1, 1], zeroline=False),
                              xaxis=dict(title="Pozycja (nt)", showgrid=True), height=400, plot_bgcolor="white",
                              legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1))
            st.plotly_chart(fig, use_container_width=True)

            st.subheader("Tabela zidentyfikowanych ORF")
            st.dataframe(st.session_state["orf_df"], use_container_width=True)

    # ZAKŁADKA 2: ADNOTACJA
    with tab2:
        if st.button("Uruchom Adnotację"):
            with st.spinner("DIAMOND przeszukuje Swiss-Prot..."):
                diamond_path = os.path.join(os.getcwd(), "diamond.exe")
                cmd = [diamond_path, "blastp", "-d", "sprot_db.dmnd", "-q", "proteins.faa", "-o", "annot.tsv", "-k",
                       "1", "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
                       "qend", "sstart", "send", "evalue", "bitscore", "stitle"]
                subprocess.run(cmd, check=True)
                df = pd.read_csv("annot.tsv", sep="\t",
                                 names=["query", "subject", "pident", "length", "mismatch", "gap", "qstart", "qend",
                                        "sstart", "send", "evalue", "bitscore", "stitle"])
                df["evalue_fmt"] = df["evalue"].map('{:.2e}'.format)
                df["Kategoria"] = df["subject"].apply(kategoryzuj_bialko)
                st.session_state["annot_df"] = df

        if "annot_df" in st.session_state:
            st.dataframe(st.session_state["annot_df"][
                             ["query", "subject", "stitle", "Kategoria", "pident", "bitscore", "evalue_fmt"]],
                         use_container_width=True)
        if "annot_df" in st.session_state:
            # --- DODAJ TO: Wykres statystyk kategorii ---
            st.subheader("Statystyki Adnotacji")
            cat_counts = st.session_state["annot_df"]["Kategoria"].value_counts()

            fig_pie = go.Figure(data=[go.Pie(labels=cat_counts.index, values=cat_counts.values, hole=0.3)])
            fig_pie.update_layout(height=400, margin=dict(t=0, b=0, l=0, r=0))
            st.plotly_chart(fig_pie, use_container_width=True)
            # -------------------------------------------
            # ZAKŁADKA 3: PREDYKCJA HOSTA (SYSTEM ZINTEGROWANY)
            with tab3:
                st.subheader("Host Evidence Dashboard (Wieloskładnikowy V4)")

                # 1. Konfiguracja wag (z podziałem na endolizyny i holiny)
                PROTEIN_WEIGHTS = {
                    "Receptor Binding Protein (RBP)": 15,
                    "Tail Fiber": 12,
                    "Tail Spike": 12,
                    "Depolimerazy": 10,
                    "Endolizyny": 8,
                    "Baseplate": 6,
                    "Holiny": 2,
                    "Kapsyd": 3,
                    "Białka portalu": 3,
                    "Polimerazy": 2,
                    "Terminazy": 2,
                    "Integrazy": 2,
                    "Pozostałe": 0.5,
                    "Brak adnotacji (ORFan)": 0.1
                }

                # 2. Definicja pokrewieństwa (dziedziczenie dowodów)
                RELATED_GENERA = {
                    "escherichia": ["salmonella", "shigella"],
                    "salmonella": ["escherichia", "shigella"],
                    "shigella": ["escherichia", "salmonella"]
                }


                # Funkcja pomocnicza do kategoryzacji (z poprawką na Endolizyny)
                def get_category(opis):
                    opis = str(opis).lower()
                    if any(x in opis for x in ["rbp", "receptor binding"]): return "Receptor Binding Protein (RBP)"
                    if any(x in opis for x in ["tail fiber", "tailfibre"]): return "Tail Fiber"
                    if any(x in opis for x in ["tail spike"]): return "Tail Spike"
                    if any(x in opis for x in ["endolysin", "lysin"]): return "Endolizyny"  # Rozdzielone!
                    if any(x in opis for x in ["holin"]): return "Holiny"
                    if any(x in opis for x in ["depolymerase"]): return "Depolimerazy"
                    if any(x in opis for x in ["baseplate"]): return "Baseplate"
                    if any(x in opis for x in ["capsid", "major capsid"]): return "Kapsyd"
                    if any(x in opis for x in ["portal"]): return "Białka portalu"
                    if any(x in opis for x in ["polymerase"]): return "Polimerazy"
                    if any(x in opis for x in ["terminase"]): return "Terminazy"
                    if any(x in opis for x in ["integrase"]): return "Integrazy"
                    return "Pozostałe"


                if "annot_df" in st.session_state:
                    # Inicjalizacja statystyk dla każdego znanego hosta
                    host_stats = {genus: {'score': 0, 'categories': set()} for genus in KNOWN_GENERA}

                    for _, row in st.session_state["annot_df"].iterrows():
                        title = str(row["stitle"]).lower()
                        kategoria = get_category(title)
                        waga = PROTEIN_WEIGHTS.get(kategoria, 0.1)

                        # Obliczanie surowego wyniku białka
                        score = (row["bitscore"] / 500) * (row["pident"] / 100) * waga

                        # 3. Logika przypisywania do hosta
                        for genus in KNOWN_GENERA:
                            # Dopasowanie bezpośrednie
                            if genus.lower() in title or f"{genus.lower()} phage" in title:
                                host_stats[genus]['score'] += score
                                host_stats[genus]['categories'].add(kategoria)

                                # Dziedziczenie dla kuzynów (40% punktów)
                                base_genus = genus.lower()
                                if base_genus in RELATED_GENERA:
                                    for cousin in RELATED_GENERA[base_genus]:
                                        cousin_cap = cousin.capitalize()
                                        if cousin_cap in host_stats:
                                            host_stats[cousin_cap]['score'] += (score * 0.4)
                                            host_stats[cousin_cap]['categories'].add(kategoria)

                    # 4. Obliczanie ostatecznego wyniku z bonusem za różnorodność
                    final_results = {}
                    for host, data in host_stats.items():
                        if data['score'] > 0:
                            diversity_bonus = 1 + math.log(len(data['categories']) + 1)
                            final_results[host] = data['score'] * diversity_bonus

                    # 5. Wyświetlanie wyników
                    if final_results:
                        total_points = sum(final_results.values())
                        sorted_hosts = sorted(final_results.items(), key=lambda item: item[1], reverse=True)

                        st.markdown("### Rozkład prawdopodobieństwa hosta")

                        for host, score in sorted_hosts:
                            prob = (score / total_points) * 100
                            if prob > 3:  # Filtrujemy nieistotne wyniki poniżej 3%
                                col1, col2, col3 = st.columns([2, 5, 1])
                                with col1: st.write(f"**{host}**")
                                with col2: st.progress(prob / 100)
                                with col3: st.write(f"{prob:.1f}%")
                    else:
                        st.warning("Brak dopasowań w bazie danych.")
                else:
                    st.info("Najpierw wykonaj adnotację w Zakładce 2.")

    # ZAKŁADKA 4: GENEROWANIE RAPORTU
    with tab4:
        st.subheader("Eksport Wyników Analizy")
        st.write("Wygeneruj podsumowanie wszystkich operacji w formie pliku PDF.")

        if st.button("Generuj PDF"):
            if "orf_df" not in st.session_state:
                st.error("Najpierw wygeneruj ORF (Zakładka 1)!")
            else:
                pdf = FPDF()
                pdf.add_page()
                pdf.set_font("Arial", 'B', 16)
                pdf.cell(0, 10, "Raport Analizy Genomowej - BioAnalyzer PRO", ln=True, align="C")
                pdf.ln(10)

                pdf.set_font("Arial", '', 12)
                pdf.cell(0, 10, usun_pl(f"Wgrany plik: {uploaded_file.name}"), ln=True)
                pdf.cell(0, 10, usun_pl(f"Liczba zidentyfikowanych ORF: {len(st.session_state['orf_df'])}"), ln=True)

                if HOST_DB_DF is not None:
                    pdf.cell(0, 10, usun_pl(f"Silnik Predykcji: Baza NCBI ({len(HOST_DB_DF):,} bialek)"), ln=True)

                if "annot_df" in st.session_state:
                    pdf.ln(5)
                    pdf.set_font("Arial", 'B', 12)
                    pdf.cell(0, 10, usun_pl("Adnotacja funkcjonalna - Kategorie:"), ln=True)
                    pdf.set_font("Arial", '', 12)

                    counts = st.session_state["annot_df"]["Kategoria"].value_counts()
                    for cat, count in counts.items():
                        pdf.cell(0, 8, usun_pl(f"- {cat}: {count} bialek"), ln=True)

                if "host_scores" in st.session_state and st.session_state["host_scores"]:
                    pdf.ln(5)
                    pdf.set_font("Arial", 'B', 12)
                    pdf.cell(0, 10, usun_pl("Prawdopodobny gospodarz (Wyniki z Evidence V2):"), ln=True)
                    pdf.set_font("Arial", '', 12)

                    scores = st.session_state["host_scores"]
                    total = sum(scores.values())
                    for host, pts in scores.items():
                        prob = (pts / total) * 100
                        pdf.cell(0, 8, usun_pl(f"- {host}: {prob:.1f}% dopasowan"), ln=True)

                pdf.output("raport_bioanalyzer.pdf")
                with open("raport_bioanalyzer.pdf", "rb") as pdf_file:
                    st.download_button(
                        label="Pobierz Raport PDF",
                        data=pdf_file,
                        file_name="Raport_BioAnalyzer.pdf",
                        mime="application/pdf"
                    )
                st.success("Raport został wygenerowany! Kliknij przycisk wyżej, aby go pobrać.")

else:
    st.info("Wgraj plik FASTA w panelu bocznym, aby rozpocząć analizę.")