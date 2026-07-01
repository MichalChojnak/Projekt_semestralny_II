import tkinter as tk
from tkinter import filedialog, messagebox
import subprocess
import os
import threading

# --- KONFIGURACJA ŚCIEŻEK (Podaj tutaj swoje ścieżki do plików .exe) ---
# Użyj "r" przed cudzysłowem, aby Windows nie traktował backslashy jako znaków specjalnych
MMSEQS_PATH = r"C:\BioTools\mmseqs\bin"
HMMER_PATH = r"C:\BioTools\hmmer"


# ----------------------------------------------------------------------

class AnnotatorApp:
    # ... (kod klasy pozostaje taki sam, zmieniamy tylko metodę process) ...

    def process(self, faa_path):
        self.status.config(text="Status: Uruchamiam MMseqs2...", fg="blue")

        # Budujemy pełną komendę
        cmd_mmseqs = f'"{MMSEQS_PATH}" easy-search "{faa_path}" "database/phrogs_db" "results_mmseqs.tsv" "tmp"'
        subprocess.run(cmd_mmseqs, shell=True)

        self.status.config(text="Status: Uruchamiam HMMER...", fg="blue")

        cmd_hmmer = f'"{HMMER_PATH}" --tblout "results_hmmer.tsv" "database/Pfam-A.hmm" "{faa_path}"'
        subprocess.run(cmd_hmmer, shell=True)

        self.status.config(text="Status: Zakończono!", fg="green")
        messagebox.showinfo("Sukces", "Adnotacja zakończona!")