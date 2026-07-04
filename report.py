# report.py
from fpdf import FPDF


def usun_pl(tekst):
    zamienniki = {'ą': 'a', 'ć': 'c', 'ę': 'e', 'ł': 'l', 'ń': 'n', 'ó': 'o', 'ś': 's', 'ź': 'z', 'ż': 'z',
                  'Ą': 'A', 'Ć': 'C', 'Ę': 'E', 'Ł': 'L', 'Ń': 'N', 'Ó': 'O', 'Ś': 'S', 'Ź': 'Z', 'Ż': 'Z'}
    for pl, asc in zamienniki.items():
        tekst = tekst.replace(pl, asc)
    return tekst


def generate_pdf_report(filename, orf_df, annot_df, host_results):
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", 'B', 16)
    pdf.cell(0, 10, "Raport Analizy Genomowej - BioAnalyzer PRO", ln=True, align="C")
    pdf.ln(10)

    pdf.set_font("Arial", '', 12)
    pdf.cell(0, 10, usun_pl(f"Wgrany plik: {filename}"), ln=True)
    pdf.cell(0, 10, usun_pl(f"Liczba zidentyfikowanych ORF: {len(orf_df)}"), ln=True)

    if annot_df is not None:
        pdf.ln(5)
        pdf.set_font("Arial", 'B', 12)
        pdf.cell(0, 10, usun_pl("Podsumowanie Adnotacji:"), ln=True)
        pdf.set_font("Arial", '', 12)
        counts = annot_df["Kategoria"].value_counts()
        for cat, count in counts.items():
            pdf.cell(0, 8, usun_pl(f"- {cat}: {count} bialek"), ln=True)

    if host_results:
        pdf.ln(5)
        pdf.set_font("Arial", 'B', 12)
        pdf.cell(0, 10, usun_pl("Predykcja gospodarza (gatunek):"), ln=True)
        pdf.set_font("Arial", '', 12)

        total_score = sum(data["score"] for host, data in host_results)
        for host, data in host_results:
            prob = (data["score"] / total_score) * 100
            if prob > 1.0:  # Pokazujemy gospodarzy z prawdopodobienstwem > 1%
                pdf.cell(0, 8, usun_pl(f"- {host}: {prob:.1f}%"), ln=True)

    pdf.output("raport_bioanalyzer.pdf")
    return "raport_bioanalyzer.pdf"