# report.py
import os
from fpdf import FPDF


def generate_pdf_report(filename, orf_df, annot_df, host_results, chart_path=None):
    pdf = FPDF()
    pdf.add_page()

    # ŁADOWANIE CZCIONKI
    pdf.add_font("ArialPL", "", r"C:\Windows\Fonts\arial.ttf")

    has_bold = False
    if os.path.exists(r"C:\Windows\Fonts\arialbd.ttf"):
        pdf.add_font("ArialPL", "B", r"C:\Windows\Fonts\arialbd.ttf")
        has_bold = True

    bold_style = "B" if has_bold else ""

    # 1. NAGŁÓWEK
    pdf.set_font("ArialPL", style=bold_style, size=18)
    pdf.cell(0, 10, "Raport Analizy Genomowej - Fagometr", ln=True, align="C")
    pdf.line(10, 20, 200, 20)
    pdf.ln(10)

    # 2. INFORMACJE OGÓLNE
    pdf.set_font("ArialPL", style=bold_style, size=12)
    pdf.cell(0, 10, "1. Podsumowanie pliku i detekcji ORF:", ln=True)

    pdf.set_font("ArialPL", style="", size=11)
    pdf.cell(0, 8, f"Wgrany plik: {filename}", ln=True)

    liczba_orf = len(orf_df) if orf_df is not None else 0
    pdf.cell(0, 8, f"Liczba zidentyfikowanych ORF (potencjalnych genów): {liczba_orf}", ln=True)

    liczba_hitow = len(annot_df) if annot_df is not None else 0
    pdf.cell(0, 8, f"Liczba wszystkich dopasowań w bazie (DIAMOND): {liczba_hitow}", ln=True)
    pdf.ln(5)

    # 3. PODSUMOWANIE ADNOTACJI
    if annot_df is not None:
        pdf.set_font("ArialPL", style=bold_style, size=12)
        pdf.cell(0, 10, "2. Adnotacje Funkcjonalne (kategorie białek):", ln=True)

        pdf.set_font("ArialPL", style="", size=11)
        counts = annot_df["Kategoria"].value_counts()
        for cat, count in counts.items():
            pdf.cell(0, 6, f" - {cat}: {count} dopasowań", ln=True)
        pdf.ln(5)

        if chart_path and os.path.exists(chart_path):
            pdf.image(chart_path, x=25, w=160)
            pdf.ln(5)

    # 4. TABELA: PREDYKCJA GOSPODARZA
    if host_results:
        pdf.add_page()
        pdf.set_font("ArialPL", style=bold_style, size=12)
        pdf.cell(0, 10, "3. Predykcja Gospodarza (Ranking):", ln=True)
        pdf.ln(5)

        pdf.set_font("ArialPL", style=bold_style, size=11)
        pdf.cell(90, 10, "Gospodarz (Host)", border=1, align="C")
        pdf.cell(45, 10, "Zgodność", border=1, align="C")
        pdf.cell(45, 10, "E-value", border=1, ln=True, align="C")

        pdf.set_font("ArialPL", style="", size=11)
        total_score = sum(item.get("Score", 0) for item in host_results)

        for item in host_results:
            host = item.get("Host", "Nieznany")
            score = item.get("Score", 0)
            e_val = item.get("E-value", 0)

            prob = (score / total_score) * 100 if total_score > 0 else 0

            if prob > 0.5:
                pdf.cell(90, 10, host, border=1)
                pdf.cell(45, 10, f"{prob:.1f}%", border=1, align="C")
                pdf.cell(45, 10, f"{e_val:.2e}", border=1, ln=True, align="C")

    # STOPKA DOKUMENTU
    pdf.set_y(-15)
    pdf.set_font("ArialPL", style="", size=8)
    pdf.cell(0, 10, "Dokument wygenerowany automatycznie przez system Fagometr", align="C")

    # ZAPIS
    output_filename = "raport_analizy.pdf"
    pdf.output(output_filename)
    return output_filename