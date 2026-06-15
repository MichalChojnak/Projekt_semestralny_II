import tkinter as tk
from tkinter import filedialog, messagebox
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from fpdf import FPDF, XPos, YPos  # Dodaliśmy importy XPos i YPos dla nowego FPDF
from pathlib import Path
import re


# --- 1. SILNIK BIOINFORMATYCZNY (Szukanie ORF) ---

def find_orfs_in_record(seq_record, min_aa_length):
    """
    Skanuje rekord FASTA w 6 ramkach odczytu (+ i -).
    Zwraca listę słowników z informacjami o znalezionych ORF.
    """
    orfs = []
    seq = seq_record.seq
    seq_len = len(seq)

    # Przeszukujemy nić wiodącą (+1) oraz opóźnioną (-1)
    for strand, nuc_seq in [(1, seq), (-1, seq.reverse_complement())]:
        # 3 ramki odczytu dla każdej nici (0, 1, 2)
        for frame in range(3):
            # ROZWIĄZANIE BiopythonWarning: Przycinamy sekwencję, aby jej długość była wielokrotnością 3
            seq_to_translate = nuc_seq[frame:]
            trim_len = (len(seq_to_translate) // 3) * 3

            # Tłumaczymy dokładnie dociętą sekwencję (Tabela 11)
            translated = str(seq_to_translate[:trim_len].translate(table=11))

            start_aa_pos = 0
            # Dzielimy białko kodonami stop (*)
            for match in re.finditer(r'([A-Z]*?)\*', translated):
                protein_chunk = match.group(1)

                # Szukamy pierwszego kodonu Start ('M')
                m_idx = protein_chunk.find('M')
                if m_idx != -1:
                    orf_protein = protein_chunk[m_idx:]

                    # Jeśli ORF jest wystarczająco długi, zapisujemy go
                    if len(orf_protein) >= min_aa_length:
                        # Obliczanie fizycznych współrzędnych w genomie
                        aa_start_in_frame = start_aa_pos + m_idx
                        nt_start_relative = frame + (aa_start_in_frame * 3)
                        nt_end_relative = nt_start_relative + (len(orf_protein) * 3) + 3  # +3 za stop kodon

                        if strand == 1:
                            real_start = nt_start_relative + 1
                            real_end = nt_end_relative
                        else:
                            # Przeliczenie współrzędnych dla nici opóźnionej
                            real_start = seq_len - nt_end_relative + 1
                            real_end = seq_len - nt_start_relative

                        orfs.append({
                            "Seq_ID": seq_record.id,
                            "Strand": "+" if strand == 1 else "-",
                            "Start": real_start,
                            "End": real_end,
                            "Length (AA)": len(orf_protein),
                            "Protein Sequence": orf_protein
                        })

                start_aa_pos = match.end()

    return sorted(orfs, key=lambda x: x["Start"])


# --- 2. GENEROWANIE RAPORTÓW ---

def generate_excel(orfs, output_path):
    df = pd.DataFrame(orfs)
    df.to_excel(output_path, index=False, sheet_name="Found ORFs")


def generate_pdf(orfs, output_path, min_len):
    """Generuje raport PDF, w 100% zgodny z najnowszym standardem fpdf2."""
    pdf = FPDF()
    pdf.add_page()

    # Zmieniono "Arial" na domyślną czcionkę "Helvetica"
    pdf.set_font("Helvetica", 'B', 16)
    # Zmieniono txt= na text= oraz ln=True na new_x/new_y
    pdf.cell(0, 10, text="Raport z poszukiwania ORF", new_x=XPos.LMARGIN, new_y=YPos.NEXT, align='C')

    pdf.set_font("Helvetica", size=11)
    pdf.cell(0, 10, text=f"Liczba znalezionych ramek: {len(orfs)}", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.cell(0, 10, text=f"Minimalna dlugosc (aminokwasy): {min_len}", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.ln(5)

    pdf.set_font("Courier", size=9)
    for i, orf in enumerate(orfs, 1):
        info = f"#{i} | ID: {orf['Seq_ID']} | Nic: {orf['Strand']} | Start: {orf['Start']} | Stop: {orf['End']} | Dlugosc: {orf['Length (AA)']} AA"
        pdf.cell(0, 6, text=info, new_x=XPos.LMARGIN, new_y=YPos.NEXT)

        seq_preview = orf['Protein Sequence'][:50] + "..." if len(orf['Protein Sequence']) > 50 else orf[
            'Protein Sequence']
        pdf.cell(0, 6, text=f"Seq: {seq_preview}", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        pdf.ln(2)

    pdf.output(output_path)


# --- 3. INTERFEJS UŻYTKOWNIKA (GUI) ---

def run_analysis():
    fasta_path = entry_file.get()
    if not fasta_path:
        messagebox.showwarning("Błąd", "Najpierw wybierz plik FASTA!")
        return

    try:
        min_aa = int(entry_len.get())
        if min_aa < 10:
            raise ValueError
    except ValueError:
        messagebox.showwarning("Błąd", "Minimalna długość musi być liczbą całkowitą.")
        return

    try:
        records = list(SeqIO.parse(fasta_path, "fasta"))
        if not records:
            messagebox.showerror("Błąd", "Plik jest pusty lub uszkodzony.")
            return

        all_orfs = []
        for record in records:
            orfs = find_orfs_in_record(record, min_aa)
            all_orfs.extend(orfs)

        if not all_orfs:
            messagebox.showinfo("Wynik", "Nie znaleziono żadnych ORF.")
            return

        out_dir = Path("results_orf")
        out_dir.mkdir(exist_ok=True)

        base_name = Path(fasta_path).stem
        excel_path = out_dir / f"{base_name}_orfs.xlsx"
        pdf_path = out_dir / f"{base_name}_orfs.pdf"

        generate_excel(all_orfs, excel_path)
        generate_pdf(all_orfs, pdf_path, min_aa)

        messagebox.showinfo("Sukces",
                            f"Analiza zakończona!\nZnaleziono {len(all_orfs)} ORF.\nRaporty w folderze: results_orf/")

    except Exception as e:
        messagebox.showerror("Wystąpił błąd", str(e))


root = tk.Tk()
root.title("Phage ORF Explorer")
root.geometry("400x250")

tk.Label(root, text="Genom faga (.fasta):", font=("Arial", 10, "bold")).pack(pady=(15, 5))
frame_file = tk.Frame(root)
frame_file.pack()
entry_file = tk.Entry(frame_file, width=40)
entry_file.pack(side=tk.LEFT, padx=5)
tk.Button(frame_file, text="Wybierz", command=lambda: entry_file.insert(0, filedialog.askopenfilename(
    filetypes=[("FASTA files", "*.fasta *.fa")]))).pack(side=tk.LEFT)

tk.Label(root, text="Min. długość białka (aminokwasy):", font=("Arial", 10)).pack(pady=(15, 5))
entry_len = tk.Entry(root, width=10, justify='center')
entry_len.insert(0, "50")
entry_len.pack()

tk.Button(root, text="ZNAJDŹ ORF I GENERUJ RAPORTY", bg="#28a745", fg="white", font=("Arial", 10, "bold"),
          command=run_analysis).pack(pady=25)

root.mainloop()