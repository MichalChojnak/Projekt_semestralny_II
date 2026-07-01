import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from Bio import SeqIO
import pandas as pd
from fpdf import FPDF, XPos, YPos
from pathlib import Path
import re

# --- 1. SILNIK BIOINFORMATYCZNY ---

def find_orfs_in_record(seq_record, min_aa_length):
    orfs = []
    seq = seq_record.seq
    seq_len = len(seq)
    for strand, nuc_seq in [(1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            seq_to_translate = nuc_seq[frame:]
            trim_len = (len(seq_to_translate) // 3) * 3
            translated = str(seq_to_translate[:trim_len].translate(table=11))
            start_aa_pos = 0
            for match in re.finditer(r'([A-Z]*?)\*', translated):
                protein_chunk = match.group(1)
                m_idx = protein_chunk.find('M')
                if m_idx != -1:
                    orf_protein = protein_chunk[m_idx:]
                    if len(orf_protein) >= min_aa_length:
                        aa_start_in_frame = start_aa_pos + m_idx
                        nt_start_relative = frame + (aa_start_in_frame * 3)
                        nt_end_relative = nt_start_relative + (len(orf_protein) * 3) + 3
                        real_start = nt_start_relative + 1 if strand == 1 else seq_len - nt_end_relative + 1
                        real_end = nt_end_relative if strand == 1 else seq_len - nt_start_relative
                        orfs.append({
                            "Seq_ID": seq_record.id, "Strand": "+" if strand == 1 else "-",
                            "Start": real_start, "End": real_end,
                            "Length (AA)": len(orf_protein), "Protein Sequence": orf_protein
                        })
                start_aa_pos = match.end()
    return sorted(orfs, key=lambda x: x["Start"])

# --- 2. GENEROWANIE PLIKÓW ---

def save_orfs_to_fasta(orfs, output_path):
    with open(output_path, "w") as f:
        for i, orf in enumerate(orfs):
            header = f">ORF_{i+1}|{orf['Seq_ID']}|{orf['Strand']}|{orf['Start']}-{orf['End']}"
            f.write(f"{header}\n{orf['Protein Sequence']}\n")

def generate_excel(orfs, output_path):
    pd.DataFrame(orfs).to_excel(output_path, index=False)

def generate_pdf(orfs, output_path, min_len):
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Helvetica", 'B', 16)
    pdf.cell(0, 10, text="Raport z poszukiwania ORF", new_x=XPos.LMARGIN, new_y=YPos.NEXT, align='C')
    pdf.set_font("Helvetica", size=11)
    pdf.cell(0, 10, text=f"Liczba ORF: {len(orfs)}", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    for i, orf in enumerate(orfs, 1):
        pdf.cell(0, 6, text=f"#{i} | {orf['Seq_ID']} | {orf['Strand']} | {orf['Start']}-{orf['End']} | {orf['Length (AA)']} AA", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.output(output_path)

# --- 3. GUI ---

class ORFFinderApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Phage ORF Explorer Pro")
        self.root.geometry("800x600")
        self.results = []

        tk.Label(root, text="Genom faga (.fasta):").pack(pady=5)
        self.entry_file = tk.Entry(root, width=50)
        self.entry_file.pack()
        tk.Button(root, text="Wybierz plik", command=self.browse).pack()

        tk.Label(root, text="Min. długość AA:").pack(pady=5)
        self.entry_len = tk.Entry(root, width=10); self.entry_len.insert(0, "50"); self.entry_len.pack()

        tk.Button(root, text="URUCHOM ANALIZĘ", bg="#007bff", fg="white", command=self.run).pack(pady=10)

        # Tabela
        self.tree = ttk.Treeview(root, columns=("id", "strand", "start", "end", "len"), show="headings")
        for col in ("id", "strand", "start", "end", "len"): self.tree.heading(col, text=col.upper())
        self.tree.pack(fill=tk.BOTH, expand=True, padx=10)

        # Eksport
        btn_frame = tk.Frame(root); btn_frame.pack(pady=10)
        tk.Button(btn_frame, text="Eksportuj Excel", command=self.save_excel).pack(side=tk.LEFT, padx=5)
        tk.Button(btn_frame, text="Eksportuj PDF", command=self.save_pdf).pack(side=tk.LEFT, padx=5)
        tk.Button(btn_frame, text="Eksportuj FASTA (.faa)", command=self.save_fasta).pack(side=tk.LEFT, padx=5)

    def browse(self):
        file = filedialog.askopenfilename()
        if file: self.entry_file.delete(0, tk.END); self.entry_file.insert(0, file)

    def run(self):
        try:
            records = list(SeqIO.parse(self.entry_file.get(), "fasta"))
            self.results = []
            for rec in records: self.results.extend(find_orfs_in_record(rec, int(self.entry_len.get())))
            for i in self.tree.get_children(): self.tree.delete(i)
            for orf in self.results: self.tree.insert("", tk.END, values=(orf["Seq_ID"], orf["Strand"], orf["Start"], orf["End"], orf["Length (AA)"]))
        except Exception as e: messagebox.showerror("Błąd", str(e))

    def save_excel(self): generate_excel(self.results, filedialog.asksaveasfilename(defaultextension=".xlsx"))
    def save_pdf(self): generate_pdf(self.results, filedialog.asksaveasfilename(defaultextension=".pdf"), self.entry_len.get())
    def save_fasta(self): save_orfs_to_fasta(self.results, filedialog.asksaveasfilename(defaultextension=".faa"))

root = tk.Tk()
app = ORFFinderApp(root)
root.mainloop()