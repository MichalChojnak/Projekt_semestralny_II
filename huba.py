import sys
import os
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
                             QPushButton, QLabel, QSpinBox, QFileDialog, QTableWidget,
                             QTableWidgetItem, QMessageBox, QHeaderView, QGroupBox, QProgressBar)
from PyQt6.QtCore import Qt
from Bio import SeqIO
from reportlab.pdfgen import canvas


# --- KLASA STYLU Z TWOJEGO POPRZEDNIEGO PROJEKTU ---
class StatBar(QWidget):
    def __init__(self, label, default_color="#50E3C2"):
        super().__init__()
        layout = QVBoxLayout(self)
        layout.setContentsMargins(2, 2, 2, 2)
        self.lbl = QLabel(f"{label}: 0")
        self.lbl.setStyleSheet("color: #a0aabf; font-size:10px; font-weight:bold;")
        self.bar = QProgressBar()
        self.bar.setFixedHeight(6)
        self.bar.setTextVisible(False)
        self.set_color(default_color)
        layout.addWidget(self.lbl)
        layout.addWidget(self.bar)

    def set_color(self, hex_color):
        self.bar.setStyleSheet(
            f"QProgressBar {{ background-color:#1a1e24; border-radius:3px; border:1px solid #282f3a; }} QProgressBar::chunk {{ background-color:{hex_color}; border-radius:2px; }}")

    def update_val(self, val, max_val, label_text):
        self.bar.setRange(0, max(1, int(max_val)))
        self.bar.setValue(int(val))
        self.lbl.setText(label_text)


# --- GŁÓWNE OKNO ---
class PhageQCApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("🧬 Phage QC Tool")
        self.setStyleSheet("background-color: #0d1117; color: #c9d1d9;")
        self.resize(950, 700)
        self.accepted = []
        self.rejected = []

        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)

        # Panele statystyk w Twoim stylu
        stats_layout = QHBoxLayout()
        self.sb_total = StatBar("Wszystkie", "#58a6ff")
        self.sb_acc = StatBar("Zaakceptowane", "#50E3C2")
        self.sb_rej = StatBar("Odrzucone", "#ff7b72")
        for sb in [self.sb_total, self.sb_acc, self.sb_rej]: stats_layout.addWidget(sb)
        layout.addLayout(stats_layout)

        # Kontrolki
        ctrl_layout = QHBoxLayout()
        self.btn_load = QPushButton("Wgraj Pliki")
        self.btn_export_csv = QPushButton("Eksport CSV")
        self.btn_export_pdf = QPushButton("Eksport PDF")
        for btn in [self.btn_load, self.btn_export_csv, self.btn_export_pdf]:
            btn.setStyleSheet("background-color: #21262d; border: 1px solid #30363d; padding: 8px;")
            ctrl_layout.addWidget(btn)
        layout.addLayout(ctrl_layout)

        self.table = QTableWidget()
        self.table.setStyleSheet("background-color: #161b22; color: #c9d1d9; border: 1px solid #30363d;")
        self.table.setColumnCount(5)
        self.table.setHorizontalHeaderLabels(["ID", "Len", "GC%", "N%", "Status"])
        layout.addWidget(self.table)

        self.btn_load.clicked.connect(self.load_files)
        self.btn_export_csv.clicked.connect(self.export_csv)
        self.btn_export_pdf.clicked.connect(self.export_pdf)

    def load_files(self):
        files, _ = QFileDialog.getOpenFileNames(self, "Wybierz pliki")
        if not files: return
        all_recs = []
        for f in files:
            for rec in SeqIO.parse(f, "fasta"):  # Uproszczono dla przykładu
                seq = str(rec.seq).upper()
                gc = ((seq.count('G') + seq.count('C')) / len(seq)) * 100
                n = seq.count('N') / len(seq)
                r = {"id": rec.id, "len": len(seq), "gc": gc, "n": n}
                if len(seq) > 5000 and n < 0.05:
                    self.accepted.append(r)
                else:
                    self.rejected.append(r)
                all_recs.append(r)

        # Aktualizacja UI
        total = len(all_recs)
        self.sb_total.update_val(total, total, f"Wszystkie: {total}")
        self.sb_acc.update_val(len(self.accepted), total, f"OK: {len(self.accepted)}")
        self.sb_rej.update_val(len(self.rejected), total, f"Err: {len(self.rejected)}")
        self.update_table()

    def update_table(self):
        self.table.setRowCount(len(self.accepted) + len(self.rejected))
        # (Logika wstawiania danych do tabeli jak poprzednio)

    def export_csv(self):
        path, _ = QFileDialog.getSaveFileName(self, "Zapisz CSV", "", "*.csv")
        if path:
            with open(path, 'w') as f:
                f.write("ID,Length,GC,N\n")
                for r in self.accepted: f.write(f"{r['id']},{r['len']},{r['gc']:.2f},{r['n']:.4f}\n")

    def export_pdf(self):
        path, _ = QFileDialog.getSaveFileName(self, "Zapisz PDF", "", "*.pdf")
        if path:
            c = canvas.Canvas(path)
            c.drawString(100, 800, "Raport QC Sekwencji")
            for i, r in enumerate(self.accepted[:20]):  # Pierwsze 20 dla testu
                c.drawString(100, 780 - (i * 20), f"{r['id']} | Len: {r['len']} | GC: {r['gc']:.1f}%")
            c.save()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = PhageQCApp()
    window.show()
    sys.exit(app.exec())