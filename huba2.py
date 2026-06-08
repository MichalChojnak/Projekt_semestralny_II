import sys
import os
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout,
                             QHBoxLayout, QPushButton, QLabel, QSpinBox,
                             QDoubleSpinBox, QFileDialog, QTableWidget,
                             QTableWidgetItem, QMessageBox, QHeaderView, QGroupBox)
from PyQt6.QtCore import Qt
from Bio import SeqIO


# 1. HUBA
# zawartość %N
def n_content(seq: str) -> float:
    return seq.count('N') / len(seq) if len(seq) > 0 else 0.0


# zawartość %GC
def gc_content(seq: str) -> float:
    seq_upper = seq.upper()
    g, c = seq_upper.count('G'), seq_upper.count('C')
    return ((g + c) / len(seq_upper)) * 100 if len(seq_upper) > 0 else 0.0


# wczytywanie pliku
def process_file(file_path: str) -> list[dict]:
    records = []
    ext = file_path.lower()

    if ext.endswith((".fasta", ".fa", ".fna")):
        fmt = "fasta"
    elif ext.endswith((".fastq", ".fq")):
        fmt = "fastq"
    elif ext.endswith((".gbk", ".gb")):
        fmt = "genbank"
    else:
        return records

    try:
        for record in SeqIO.parse(file_path, fmt):
            clean_seq = str(record.seq).upper().replace("-", "")
            records.append({
                "id": record.id,
                "length": len(clean_seq),
                "gc_pct": gc_content(clean_seq),
                "n_pct": n_content(clean_seq),
                "sequence": clean_seq,
                "source": os.path.basename(file_path)
            })
    except Exception as e:
        print(f"Błąd czytania {file_path}: {e}")
    return records


# 2. INTERFEJS UŻYTKOWNIKA (PyQt6)

class PhageQCApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("🧬 Phage Host Predictor - Moduł QC")
        self.resize(900, 600)

        # Zmienne do przechowywania wyników
        self.raw_records = []
        self.accepted_records = []
        self.rejected_records = []

        # Główny widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)

        # --- PANEL GÓRNY (Kontrolki) ---
        top_panel = QHBoxLayout()

        # Grupa: Parametry
        params_group = QGroupBox("⚙️ Parametry Filtrowania")
        params_layout = QHBoxLayout()

        params_layout.addWidget(QLabel("Min. długość (bp):"))
        self.min_len_spin = QSpinBox()
        self.min_len_spin.setRange(1, 1000000)
        self.min_len_spin.setValue(5000)
        self.min_len_spin.setSingleStep(500)
        params_layout.addWidget(self.min_len_spin)

        params_layout.addWidget(QLabel("Maks. 'N' (%):"))
        self.max_n_spin = QDoubleSpinBox()
        self.max_n_spin.setRange(0.0, 100.0)
        self.max_n_spin.setValue(5.0)
        params_layout.addWidget(self.max_n_spin)

        params_group.setLayout(params_layout)
        top_panel.addWidget(params_group)

        # Grupa: Akcje
        actions_group = QGroupBox("📂 Akcje")
        actions_layout = QHBoxLayout()

        self.btn_load = QPushButton("Wgraj pliki...")
        self.btn_load.clicked.connect(self.load_files)
        actions_layout.addWidget(self.btn_load)

        # DODANY PRZYCISK ANALIZY
        self.btn_run = QPushButton("Analizuj (Run QC)")
        self.btn_run.setEnabled(False)
        self.btn_run.clicked.connect(self.run_qc)
        actions_layout.addWidget(self.btn_run)

        self.btn_export = QPushButton("Eksportuj poprawne (FASTA)")
        self.btn_export.setEnabled(False)
        self.btn_export.clicked.connect(self.export_fasta)
        actions_layout.addWidget(self.btn_export)

        actions_group.setLayout(actions_layout)
        top_panel.addWidget(actions_group)

        main_layout.addLayout(top_panel)

        # --- PANEL ŚRODKOWY (Statystyki) ---
        stats_layout = QHBoxLayout()
        self.lbl_total = QLabel("Wszystkie: 0")
        self.lbl_accepted = QLabel("Zaakceptowane: 0")
        self.lbl_rejected = QLabel("Odrzucone: 0")

        self.lbl_accepted.setStyleSheet("color: green; font-weight: bold;")
        self.lbl_rejected.setStyleSheet("color: red; font-weight: bold;")

        stats_layout.addWidget(self.lbl_total)
        stats_layout.addWidget(self.lbl_accepted)
        stats_layout.addWidget(self.lbl_rejected)
        stats_layout.addStretch()
        main_layout.addLayout(stats_layout)

        # --- PANEL DOLNY (Tabela wyników) ---
        self.table = QTableWidget()
        self.table.setColumnCount(6)
        self.table.setHorizontalHeaderLabels(["ID Sekwencji", "Długość", "GC (%)", "N (%)", "Status", "Powód / Źródło"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        main_layout.addWidget(self.table)

    # 3. AKCJE (Obsługa przycisków)

    def load_files(self):
        files, _ = QFileDialog.getOpenFileNames(
            self,
            "Wybierz pliki sekwencji",
            "",
            "Pliki sekwencji (*.fasta *.fa *.fna *.fastq *.fq *.gbk);;Wszystkie pliki (*)"
        )
        if not files:
            return

        self.raw_records = []
        for file_path in files:
            self.raw_records.extend(process_file(file_path))

        if not self.raw_records:
            QMessageBox.warning(self, "Błąd", "Nie znaleziono poprawnych sekwencji.")
            self.btn_run.setEnabled(False)
            return

        self.btn_run.setEnabled(True)
        QMessageBox.information(self, "Wczytano",
                                f"Wczytano {len(self.raw_records)} sekwencji. Kliknij 'Analizuj', aby wykonać QC.")

    def run_qc(self):
        """Dodatkowa funkcja do uruchamiania filtrowania ręcznie."""
        self.filter_and_display(self.raw_records)

    def filter_and_display(self, all_records):
        self.accepted_records = []
        self.rejected_records = []
        allowed_chars = set("ACGTN")

        min_len = self.min_len_spin.value()
        max_n_pct = self.max_n_spin.value() / 100.0

        for r in all_records:
            invalid_chars = set(r["sequence"]) - allowed_chars
            if invalid_chars:
                r["reject_reason"] = f"Złe znaki: {', '.join(invalid_chars)}"
                self.rejected_records.append(r)
            elif r["length"] < min_len:
                r["reject_reason"] = "Za krótka"
                self.rejected_records.append(r)
            elif r["n_pct"] > max_n_pct:
                r["reject_reason"] = "Zbyt dużo N"
                self.rejected_records.append(r)
            else:
                self.accepted_records.append(r)

        self.lbl_total.setText(f"Wszystkie: {len(all_records)}")
        self.lbl_accepted.setText(f"Zaakceptowane: {len(self.accepted_records)}")
        self.lbl_rejected.setText(f"Odrzucone: {len(self.rejected_records)}")

        self.table.setRowCount(len(all_records))
        row = 0
        for r in self.accepted_records:
            self._add_table_row(row, r, "Zaakceptowano", r["source"], Qt.GlobalColor.darkGreen)
            row += 1

        for r in self.rejected_records:
            self._add_table_row(row, r, "Odrzucono", r["reject_reason"], Qt.GlobalColor.red)
            row += 1

        self.btn_export.setEnabled(len(self.accepted_records) > 0)

    def _add_table_row(self, row, record, status_text, reason_text, color):
        self.table.setItem(row, 0, QTableWidgetItem(record["id"]))
        self.table.setItem(row, 1, QTableWidgetItem(str(record["length"])))
        self.table.setItem(row, 2, QTableWidgetItem(f"{record['gc_pct']:.2f}"))
        self.table.setItem(row, 3, QTableWidgetItem(f"{record['n_pct'] * 100:.2f}"))

        status_item = QTableWidgetItem(status_text)
        status_item.setForeground(color)
        self.table.setItem(row, 4, status_item)
        self.table.setItem(row, 5, QTableWidgetItem(reason_text))

    def export_fasta(self):
        save_path, _ = QFileDialog.getSaveFileName(
            self,
            "Zapisz oczyszczony plik FASTA",
            "zwalidowane_fagi.fasta",
            "FASTA Files (*.fasta);;All Files (*)"
        )

        if save_path:
            try:
                with open(save_path, "w", encoding="utf-8") as f:
                    for r in self.accepted_records:
                        f.write(f">{r['id']} source:{r['source']} len:{r['length']}\n")
                        for i in range(0, r["length"], 80):
                            f.write(r["sequence"][i:i + 80] + "\n")

                QMessageBox.information(self, "Sukces", f"Zapisano poprawnie {len(self.accepted_records)} sekwencji!")
            except Exception as e:
                QMessageBox.critical(self, "Błąd zapisu", str(e))


if __name__ == '__main__':
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    window = PhageQCApp()
    window.show()
    sys.exit(app.exec())