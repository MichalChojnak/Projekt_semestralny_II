import tkinter as tk
from tkinter import filedialog, messagebox, ttk  # Dodano import ttk dla tabeli
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from fpdf import FPDF, XPos, YPos
from pathlib import Path
import re


# --- 1. SILNIK BIOINFORMATYCZNY (Szukanie ORF) ---

def find_orfs_in_record(seq_record, min_aa_length):
    """Skanuje rekord FASTA w 6 ramkach odczytu (+ i -)."""
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

                        if strand == 1:
                            real_start = nt_start_relative + 1
                            real_end = nt_end_relative
                        else:
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


# --- 2. FUNKCJE GENERUJĄCE ---

def generate_excel(orfs, output_path):
    df = pd.DataFrame(orfs)
    df.to_excel(output_path, index=False, sheet_name="Found ORFs")


def generate_pdf(orfs, output_path, min_len):
    pdf = FPDF()
    pdf.add_page()

    pdf.set_font("Helvetica", 'B', 16)
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


# --- 3. LOGIKA INTEGRACJI I WYŚWIETLANIA ---

def get_analyzed_orfs():
    """Wspólna funkcja pomocnicza walidująca wejście i uruchamiająca analizę."""
    fasta_path = entry_file.get()
    if not fasta_path:
        messagebox.showwarning("Błąd", "Najpierw wybierz plik FASTA!")
        return None, None

    try:
        min_aa = int(entry_len.get())
        if min_aa < 10:
            raise ValueError
    except ValueError:
        messagebox.showwarning("Błąd", "Minimalna długość musi być liczbą całkowitą (minimum 10).")
        return None, None

    try:
        records = list(SeqIO.parse(fasta_path, "fasta"))
        if not records:
            messagebox.showerror("Błąd", "Plik jest pusty lub uszkodzony.")
            return None, None

        all_orfs = []
        for record in records:
            orfs = find_orfs_in_record(record, min_aa)
            all_orfs.extend(orfs)

        if not all_orfs:
            messagebox.showinfo("Wynik", "Nie znaleziono żadnych ORF spełniających kryteria.")
            return None, None

        return all_orfs, min_aa
    except Exception as e:
        messagebox.showerror("Wystąpił błąd analizy", str(e))
        return None, None


def display_results_in_gui():
    """Pobiera wyniki analizy i wprowadza je do tabeli w GUI."""
    # Czyszczenie starych wyników z tabeli
    for item in tree.get_children():
        tree.delete(item)

    orfs, _ = get_analyzed_orfs()
    if not orfs:
        lbl_status.config(text="Status: Brak wyników do wyświetlenia", fg="red")
        return

    # Wstawianie nowych wierszy do tabeli
    for idx, orf in enumerate(orfs, 1):
        # Skracamy podgląd sekwencji do tabeli
        seq_short = orf["Protein Sequence"][:40] + "..." if len(orf["Protein Sequence"]) > 40 else orf[
            "Protein Sequence"]

        tree.insert("", tk.END, values=(
            idx,
            orf["Seq_ID"],
            orf["Strand"],
            orf["Start"],
            orf["End"],
            orf["Length (AA)"],
            seq_short
        ))

    lbl_status.config(text=f"Status: Znaleziono i wyświetlono {len(orfs)} ORF-ów.", fg="green")


def export_to_excel():
    """Uruchamia analizę i pyta o miejsce zapisu pliku Excel."""
    orfs, _ = get_analyzed_orfs()
    if not orfs:
        return

    base_name = Path(entry_file.get()).stem
    save_path = filedialog.asksaveasfilename(
        defaultextension=".xlsx",
        filetypes=[("Skoroszyt Excel", "*.xlsx")],
        initialfile=f"{base_name}_orfs.xlsx",
        title="Wybierz miejsce zapisu raportu Excel"
    )

    if save_path:
        generate_excel(orfs, save_path)
        messagebox.showinfo("Sukces", f"Raport Excel został pomyślnie zapisany w:\n{save_path}")


def export_to_pdf():
    """Uruchamia analizę i pyta o miejsce zapisu pliku PDF."""
    orfs, min_aa = get_analyzed_orfs()
    if not orfs:
        return

    base_name = Path(entry_file.get()).stem
    save_path = filedialog.asksaveasfilename(
        defaultextension=".pdf",
        filetypes=[("Dokument PDF", "*.pdf")],
        initialfile=f"{base_name}_orfs.pdf",
        title="Wybierz miejsce zapisu raportu PDF"
    )

    if save_path:
        generate_pdf(orfs, save_path, min_aa)
        messagebox.showinfo("Sukces", f"Raport PDF został pomyślnie zapisany w:\n{save_path}")


# --- 4. INTERFEJS UŻYTKOWNIKA ---

root = tk.Tk()
root.title("Phage ORF Explorer Pro")
root.geometry("750x550")  # Powiększono okno, aby zmieścić tabelę z wynikami

# Sekcja wyboru pliku
tk.Label(root, text="Genom faga (.fasta):", font=("Arial", 10, "bold")).pack(pady=(15, 5))
frame_file = tk.Frame(root)
frame_file.pack()
entry_file = tk.Entry(frame_file, width=50)
entry_file.pack(side=tk.LEFT, padx=5)


def select_input_file():
    selected = filedialog.askopenfilename(filetypes=[("Pliki FASTA", "*.fasta *.fa")])
    if selected:
        entry_file.delete(0, tk.END)
        entry_file.insert(0, selected)


tk.Button(frame_file, text="Wybierz plik", command=select_input_file).pack(side=tk.LEFT)

# Sekcja parametru długości
tk.Label(root, text="Min. długość białka (aminokwasy):", font=("Arial", 10)).pack(pady=(10, 5))
entry_len = tk.Entry(root, width=10, justify='center')
entry_len.insert(0, "50")
entry_len.pack()

# Przycisk głównej analizy
btn_analyze = tk.Button(
    root,
    text="Uruchom analizę i pokaż wyniki",
    bg="#007bff",
    fg="white",
    font=("Arial", 10, "bold"),
    padx=15, pady=6,
    command=display_results_in_gui
)
btn_analyze.pack(pady=15)

# --- SEKCJA WIZUALIZACJI WYNIKÓW (TABELA) ---
tk.Label(root, text="Podgląd znalezionych ORF-ów:", font=("Arial", 10, "bold")).pack(anchor="w", padx=20)

frame_table = tk.Frame(root)
frame_table.pack(fill=tk.BOTH, expand=True, padx=20, pady=5)

# Paski przewijania (Scrollbars)
scroll_y = tk.Scrollbar(frame_table, orient=tk.VERTICAL)
scroll_y.pack(side=tk.RIGHT, fill=tk.Y)

# Konfiguracja kolumn tabeli
columns = ("lp", "id", "strand", "start", "end", "len", "seq")
tree = ttk.Treeview(frame_table, columns=columns, show="headings", yscrollcommand=scroll_y.set)
scroll_y.config(command=tree.yview)

# Nagłówki i szerokości kolumn
tree.heading("lp", text="#")
tree.heading("id", text="Seq ID")
tree.heading("strand", text="Nić")
tree.heading("start", text="Start")
tree.heading("end", text="Stop")
tree.heading("len", text="Długość (AA)")
tree.heading("seq", text="Sekwencja (podgląd)")

tree.column("lp", width=40, anchor="center")
tree.column("id", width=100, anchor="w")
tree.column("strand", width=50, anchor="center")
tree.column("start", width=70, anchor="center")
tree.column("end", width=70, anchor="center")
tree.column("len", width=90, anchor="center")
tree.column("seq", width=280, anchor="w")

tree.pack(fill=tk.BOTH, expand=True)

# Etykieta statusu pod tabelą
lbl_status = tk.Label(root, text="Status: Oczekiwanie na uruchomienie analizy", font=("Arial", 9, "italic"), fg="gray")
lbl_status.pack(anchor="w", padx=20, pady=2)

# Sekcja akcji eksportu (Przeniesiona na sam dół)
frame_actions = tk.Frame(root)
frame_actions.pack(pady=15)

btn_excel = tk.Button(
    frame_actions,
    text="Eksportuj do EXCEL",
    bg="#28a745",
    fg="white",
    font=("Arial", 9, "bold"),
    padx=10, pady=5,
    command=export_to_excel
)
btn_excel.pack(side=tk.LEFT, padx=15)

btn_pdf = tk.Button(
    frame_actions,
    text="Eksportuj do PDF",
    bg="#dc3545",
    fg="white",
    font=("Arial", 9, "bold"),
    padx=10, pady=5,
    command=export_to_pdf
)
btn_pdf.pack(side=tk.LEFT, padx=15)

root.mainloop()