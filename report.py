# report.py
from fpdf import FPDF


def usun_pl(tekst):
    """Pomocnicza funkcja usuwająca polskie znaki, których standardowy FPDF nie obsługuje."""
    zamienniki = {'ą': 'a', 'ć': 'c', 'ę': 'e', 'ł': 'l', 'ń': 'n', 'ó': 'o', 'ś': 's', 'ź': 'z', 'ż': 'z',
                  'Ą': 'A', 'Ć': 'C', 'Ę': 'E', 'Ł': 'L', 'Ń': 'N', 'Ó': 'O', 'Ś': 'S', 'Ź': 'Z', 'Ż': 'Z'}
    for pl, asc in zamienniki.items():
        tekst = tekst.replace(pl, asc)
    return tekst


def generate_pdf_report(filename, orf_df, annot_df, host_results):
    pdf = FPDF()
    pdf.add_page()

    # 1. Nagłówek
    pdf.set_font("Arial", 'B', 18)
    pdf.cell(0, 10, usun_pl("Raport Analizy Genomowej - Fagometr"), ln=True, align="C")
    pdf.line(10, 20, 200, 20)  # Linia pozioma
    pdf.ln(10)

    # 2. Metadane
    pdf.set_font("Arial", 'B', 12)
    pdf.cell(0, 10, usun_pl("Informacje ogólne:"), ln=True)
    pdf.set_font("Arial", '', 11)
    pdf.cell(0, 8, usun_pl(f"Wgrany plik: {filename}"), ln=True)
    pdf.cell(0, 8, usun_pl(f"Liczba zidentyfikowanych ORF: {len(orf_df)}"), ln=True)
    pdf.ln(5)

    # 3. Podsumowanie adnotacji
    if annot_df is not None:
        pdf.set_font("Arial", 'B', 12)
        pdf.cell(0, 10, usun_pl("Podsumowanie Adnotacji Funkcjonalnych:"), ln=True)
        pdf.set_font("Arial", '', 11)
        counts = annot_df["Kategoria"].value_counts()
        for cat, count in counts.items():
            pdf.cell(0, 6, usun_pl(f" - {cat}: {count} bialek"), ln=True)
        pdf.ln(5)

    # 4. Predykcja gospodarza (NAPRAWIONA PĘTLA)
    if host_results:
        pdf.set_font("Arial", 'B', 12)
        pdf.cell(0, 10, usun_pl("Predykcja gospodarza:"), ln=True)
        pdf.set_font("Arial", '', 11)

        # Obliczenie sumy wyników (Score)
        total_score = sum(item.get("Score", 0) for item in host_results)

        # Iteracja po liście słowników
        for item in host_results:
            host = item.get("Host", "Nieznany")
            score = item.get("Score", 0)
            e_val = item.get("E-value", 0)

            prob = (score / total_score) * 100 if total_score > 0 else 0

            if prob > 0.5:  # Filtrowanie mniej istotnych wyników
                tekst = usun_pl(f" - {host}: {prob:.1f}% (E-value: {e_val:.2e})")
                pdf.cell(0, 7, tekst, ln=True)

    # Stopka
    pdf.set_y(-20)
    pdf.set_font("Arial", 'I', 8)
    pdf.cell(0, 10, usun_pl("Dokument wygenerowany automatycznie przez system Fagometr"), align="C")

    # Zapis
    output_filename = "raport_analizy.pdf"
    pdf.output(output_filename)
    return output_filename