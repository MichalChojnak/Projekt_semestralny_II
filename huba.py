mport tkinter as tk
from tkinter import filedialog, messagebox, ttk
import pandas as pd
from fpdf import FPDF, XPos, YPos
from pathlib import Path


# --- 1. LOGIKA ---

def calculate_sequence_stats(seq):
    seq = seq.upper()
    length = len(seq)
    if length == 0: return 0, 0
    n_pct = (seq.count('N') / length) * 100
    gc_pct = ((seq.count('G') + seq.count('C')) / length) * 100
    return length, n_pct, gc_pct


def parse_fasta(file_path):
    sequences = []
    with open(file_path, 'r', encoding='utf-8') as f:
        header, seq = None, []
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if header: sequences.append({"id": header, "sequence": "".join(seq)})
                header, seq = line[1:], []
            else:
                seq.append(line)
        if header: sequences.append({"id": header, "sequence": "".join(seq)})
    return sequences


# --- 2. GUI I FUNKCJE ---

data_store = []  # Przechowuje wczytane dane


def load_file():
    path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta *.fa")])
    if not path: return

    global data_store
    data_store = parse_fasta(path)

    # Czyszczenie tabeli
    for i in tree.get_children(): tree.delete(i)

    # Dodawanie danych do tabeli
    for idx, r in enumerate(data_store):
        length, n_pct, gc_pct = calculate_sequence_stats(r['sequence'])
        tree.insert("", "end", values=(idx + 1, r['id'][:20], length, f"{n_pct:.2f}%", f"{gc_pct:.2f}%"))


def save_excel():
    if not data_store: return
    path = filedialog.asksaveasfilename(defaultextension=".xlsx", filetypes=[("Excel", "*.xlsx")])
    if path:
        df = pd.DataFrame(data_store)
        if 'sequence' in df.columns: df = df.drop(columns=['sequence'])
        df.to_excel(path, index=False)
        messagebox.showinfo("Sukces", "Zapisano do Excela")


def save_pdf():
    if not data_store: return
    path = filedialog.asksaveasfilename(defaultextension=".pdf", filetypes=[("PDF", "*.pdf")])
    if path:
        pdf = FPDF()
        pdf.add_page()
        pdf.set_font("helvetica", 'B', 16)
        pdf.cell(200, 10, text="Raport Sekwencji DNA", new_x=XPos.LMARGIN, new_y=YPos.NEXT, align='C')
        pdf.set_font("helvetica", size=10)
        for r in data_store:
            l, n, g = calculate_sequence_stats(r['sequence'])
            pdf.cell(200, 8, text=f"ID: {r['id']} | Dlugosc: {l} | N: {n:.2f}% | GC: {g:.2f}%", new_x=XPos.LMARGIN,
                     new_y=YPos.NEXT)
        pdf.output(path)
        messagebox.showinfo("Sukces", "Zapisano do PDF")


# --- 3. INTERFEJS ---

root = tk.Tk()
root.title("NGS Pipeline Pro")
root.geometry("600x500")

tk.Button(root, text="Wczytaj plik FASTA", command=load_file).pack(pady=10)

# Tabela
columns = ("nr", "id", "len", "n_pct", "gc_pct")
tree = ttk.Treeview(root, columns=columns, show="headings")
tree.heading("nr", text="Nr")
tree.heading("id", text="ID Sekwencji")
tree.heading("len", text="Dlugosc")
tree.heading("n_pct", text="% N")
tree.heading("gc_pct", text="% GC")
tree.column("nr", width=30);
tree.column("len", width=60)
tree.pack(fill=tk.BOTH, expand=True, padx=10)

# Przyciski zapisu
btn_frame = tk.Frame(root)
btn_frame.pack(pady=10)
tk.Button(btn_frame, text="Zapisz Excel", command=save_excel, bg="lightblue").pack(side=tk.LEFT, padx=5)
tk.Button(btn_frame, text="Zapisz PDF", command=save_pdf, bg="lightgreen").pack(side=tk.LEFT, padx=5)

root.mainloop()