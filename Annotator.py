import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import pandas as pd
import subprocess
import os
import threading


class PhageAnnotatorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Phage Annotator Pipeline")
        self.root.geometry("600x400")

        # GUI: Wybór pliku Excel
        tk.Label(root, text="Wczytaj plik z wynikami ORF (Excel):", font=("Arial", 10, "bold")).pack(pady=10)
        self.entry_file = tk.Entry(root, width=50)
        self.entry_file.pack(pady=5)
        tk.Button(root, text="Wybierz plik", command=self.browse_file).pack()

        # Przyciski akcji
        self.btn_annotate = tk.Button(root, text="URUCHOM ADNOTACJĘ (MMseqs2 + HMMER)",
                                      bg="#007bff", fg="white", font=("Arial", 10, "bold"),
                                      command=self.start_annotation_thread)
        self.btn_annotate.pack(pady=20)

        self.status_label = tk.Label(root, text="Gotowy", fg="gray")
        self.status_label.pack()

    def browse_file(self):
        filename = filedialog.askopenfilename(filetypes=[("Excel files", "*.xlsx")])
        if filename:
            self.entry_file.delete(0, tk.END)
            self.entry_file.insert(0, filename)

    def start_annotation_thread(self):
        file_path = self.entry_file.get()
        if not file_path:
            messagebox.showwarning("Błąd", "Wybierz plik!")
            return

        # Uruchamiamy w osobnym wątku, żeby nie zawiesić GUI
        threading.Thread(target=self.run_annotation, args=(file_path,), daemon=True).start()

    def run_annotation(self, file_path):
        self.status_label.config(text="Wczytywanie danych...", fg="blue")
        df = pd.read_excel(file_path)

        # 1. Zapisujemy białka do pliku tymczasowego
        fasta_file = "to_annotate.fasta"
        with open(fasta_file, "w") as f:
            for idx, row in df.iterrows():
                f.write(f">orf_{idx}\n{row['Protein Sequence']}\n")

        # 2. MMseqs2 (wymaga zainstalowanego mmseqs w PATH)
        self.status_label.config(text="Uruchamianie MMseqs2...")
        subprocess.run(f"mmseqs easy-search {fasta_file} database/phrogs_db result_mmseqs.tsv tmp", shell=True)

        # 3. HMMER (wymaga zainstalowanego hmmer w PATH)
        self.status_label.config(text="Uruchamianie HMMER...")
        subprocess.run(f"hmmsearch --tblout result_hmmer.tsv database/pfam.hmm {fasta_file}", shell=True)

        self.status_label.config(text="Zakończono! Wyniki zapisano jako annotated_results.xlsx", fg="green")
        messagebox.showinfo("Sukces", "Adnotacja zakończona. Sprawdź folder programu.")


# Uruchomienie aplikacji
if __name__ == "__main__":
    root = tk.Tk()
    app = PhageAnnotatorApp(root)
    root.mainloop()