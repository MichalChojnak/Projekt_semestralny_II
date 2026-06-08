import sys
import os
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout,
                             QHBoxLayout, QPushButton, QLabel, QSpinBox,
                             QDoubleSpinBox, QFileDialog, QTableWidget,
                             QTableWidgetItem, QMessageBox, QHeaderView, QGroupBox,
                             QScrollArea)
from PyQt6.QtCore import Qt
from Bio import SeqIO

# Importy dodane na potrzeby nowych funkcjonalności (wykresy i raporty)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np


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
        self.setWindowTitle("Phage Host Predictor - Moduł QC + Raportowanie")
        self.resize(1200, 750)

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
        params_group = QGroupBox("Parametry Filtrowania")
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

        # Grupa: Akcje i Eksport
        actions_group = QGroupBox("Akcje i Raporty")
        actions_layout = QHBoxLayout()

        self.btn_load = QPushButton("Wgraj pliki...")
        self.btn_load.clicked.connect(self.load_files)
        actions_layout.addWidget(self.btn_load)

        self.btn_run = QPushButton("Analizuj (Run QC)")
        self.btn_run.setEnabled(False)
        self.btn_run.clicked.connect(self.run_qc)
        actions_layout.addWidget(self.btn_run)

        self.btn_export = QPushButton("Eksportuj FASTA")
        self.btn_export.setEnabled(False)
        self.btn_export.clicked.connect(self.export_fasta)
        actions_layout.addWidget(self.btn_export)

        self.btn_excel = QPushButton("Raport Excel (.xlsx)")
        self.btn_excel.setEnabled(False)
        self.btn_excel.clicked.connect(self.export_excel)
        actions_layout.addWidget(self.btn_excel)

        self.btn_pdf = QPushButton("Raport PDF (.pdf)")
        self.btn_pdf.setEnabled(False)
        self.btn_pdf.clicked.connect(self.export_pdf)
        actions_layout.addWidget(self.btn_pdf)

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

        # --- PANEL DOLNY (Tabela wyników + Wykresy z boku) ---
        bottom_layout = QHBoxLayout()

        # Lewa strona - Tabela
        self.table = QTableWidget()
        self.table.setColumnCount(7)
        self.table.setHorizontalHeaderLabels(
            ["ID Sekwencji", "Długość", "GC (%)", "N (%)", "Quality Score", "Status", "Powód / Źródło"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        bottom_layout.addWidget(self.table, stretch=2)

        # Prawa strona - Wykresy
        charts_layout = QVBoxLayout()

        # Wykres 1: Wykres Kołowy (Akceptowane vs Odrzucone)
        self.pie_fig = Figure(figsize=(4, 3), dpi=100)
        self.pie_canvas = FigureCanvas(self.pie_fig)
        charts_layout.addWidget(self.pie_canvas, stretch=1)

        # Wykres 2: Heatmapa w obszarze przewijanym (żeby obsłużyć wszystkie wyniki)
        self.scroll_area = QScrollArea()
        self.scroll_area.setWidgetResizable(True)
        self.heatmap_fig = Figure(figsize=(4, 4), dpi=100)
        self.heatmap_canvas = FigureCanvas(self.heatmap_fig)
        self.scroll_area.setWidget(self.heatmap_canvas)
        charts_layout.addWidget(self.scroll_area, stretch=2)

        bottom_layout.addLayout(charts_layout, stretch=1)
        main_layout.addLayout(bottom_layout)

        # Pierwsze rysowanie pustych wykresów startowych
        self.update_pie_chart(0, 0)
        self.update_heatmap([])

    # 3. AKCJE (Obsługa przycisków)

    def load_files(self):
        files, _ = QFileDialog.getOpenFileNames(self, "Wybierz pliki sekwencji", "",
                                                "Pliki sekwencji (*.fasta *.fa *.fna *.fastq *.fq *.gbk);;Wszystkie pliki (*)")
        if not files:
            return

        self.raw_records = []
        for file_path in files:
            self.raw_records.extend(process_file(file_path))

        if not self.raw_records:
            QMessageBox.warning(self, "Błąd", "Nie znaleziono poprawnych sekwencji w wybranych plikach.")
            self.btn_run.setEnabled(False)
            return

        self.btn_run.setEnabled(True)
        QMessageBox.information(self, "Wczytano pliki",
                                f"Pomyślnie wczytano {len(self.raw_records)} sekwencji.\nKliknij przycisk 'Analizuj (Run QC)', aby wyliczyć parametry.")

    def run_qc(self):
        if not self.raw_records:
            return
        self.filter_and_display(self.raw_records)

    def filter_and_display(self, all_records):
        self.accepted_records = []
        self.rejected_records = []
        allowed_chars = set("ACGTN")

        min_len = self.min_len_spin.value()
        max_n_pct = self.max_n_spin.value() / 100.0

        for r in all_records:
            # Algorytm wyliczania Quality Score (0-100)
            n_factor = max(0.0, 100.0 - (r["n_pct"] * 100 * 10))
            gc = r["gc_pct"]
            if 35 <= gc <= 65:
                gc_factor = 100.0
            else:
                deviation = min(abs(gc - 35), abs(gc - 65))
                gc_factor = max(0.0, 100.0 - (deviation * 4))

            len_factor = max(0.0, min(100.0, (r["length"] / 25000) * 100))

            r["quality_score"] = round((n_factor * 0.50) + (gc_factor * 0.25) + (len_factor * 0.25), 1)

            # Filtrowanie pod kątem statusów akceptacji
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

        # Aktualizacja etykiet
        self.lbl_total.setText(f"Wszystkie: {len(all_records)}")
        self.lbl_accepted.setText(f"Zaakceptowane: {len(self.accepted_records)}")
        self.lbl_rejected.setText(f"Odrzucone: {len(self.rejected_records)}")

        # Wypełnianie tabeli
        self.table.setRowCount(len(all_records))

        row = 0
        for r in self.accepted_records:
            self._add_table_row(row, r, "Zaakceptowano", r["source"], Qt.GlobalColor.darkGreen)
            row += 1

        for r in self.rejected_records:
            self._add_table_row(row, r, "Odrzucono", r["reject_reason"], Qt.GlobalColor.red)
            row += 1

        # Aktualizacja wykresów
        self.update_pie_chart(len(self.accepted_records), len(self.rejected_records))
        self.update_heatmap(all_records)

        # Aktywacja przycisków eksportu
        self.btn_export.setEnabled(len(self.accepted_records) > 0)
        self.btn_excel.setEnabled(len(all_records) > 0)
        self.btn_pdf.setEnabled(len(all_records) > 0)

    def _add_table_row(self, row, record, status_text, reason_text, color):
        self.table.setItem(row, 0, QTableWidgetItem(record["id"]))
        self.table.setItem(row, 1, QTableWidgetItem(str(record["length"])))
        self.table.setItem(row, 2, QTableWidgetItem(f"{record['gc_pct']:.2f}"))
        self.table.setItem(row, 3, QTableWidgetItem(f"{record['n_pct'] * 100:.2f}"))
        self.table.setItem(row, 4, QTableWidgetItem(str(record["quality_score"])))

        status_item = QTableWidgetItem(status_text)
        status_item.setForeground(color)
        self.table.setItem(row, 5, status_item)
        self.table.setItem(row, 6, QTableWidgetItem(reason_text))

    def update_pie_chart(self, accepted, rejected):
        self.pie_fig.clear()
        ax = self.pie_fig.add_subplot(111)

        if accepted == 0 and rejected == 0:
            ax.text(0.5, 0.5, "Brak danych.\nWgraj pliki.", ha='center', va='center', color='gray')
            ax.axis('off')
        else:
            labels_raw = ['Zaakceptowane', 'Odrzucone']
            sizes_raw = [accepted, rejected]
            colors_raw = ['#2ecc71', '#e74c3c']
            explode_raw = [0.1, 0]

            # Filtrowanie zerowych wartości, by wykres się nie psuł
            labels = [l for l, s in zip(labels_raw, sizes_raw) if s > 0]
            sizes = [s for s in sizes_raw if s > 0]
            colors = [c for c, s in zip(colors_raw, sizes_raw) if s > 0]
            explode = [e for e, s in zip(explode_raw, sizes_raw) if s > 0]

            ax.pie(sizes, explode=explode, labels=labels, colors=colors,
                   autopct='%1.1f%%', shadow=True, startangle=140, textprops={'fontsize': 9})
            ax.set_title("Podsumowanie Filtracji\n", fontsize=10, fontweight='bold')

        self.pie_fig.tight_layout()
        self.pie_canvas.draw()

    def update_heatmap(self, records):
        self.heatmap_fig.clear()
        ax = self.heatmap_fig.add_subplot(111)

        if not records:
            ax.text(0.5, 0.5, "Brak danych heatmapy.", ha='center', va='center', fontsize=10, color='gray')
            ax.axis('off')
            self.heatmap_canvas.draw()
            return

        subset = records

        # Dynamiczna wysokość wykresu
        fig_height = max(4.0, len(subset) * 0.25)
        self.heatmap_fig.set_size_inches(4.0, fig_height)

        labels = [r["id"][:12] for r in subset]

        gc_vals = [r["gc_pct"] for r in subset]
        n_vals = [r["n_pct"] * 100 for r in subset]
        qs_vals = [r["quality_score"] for r in subset]
        len_vals = [r["length"] for r in subset]  # Dodano Długość

        data_matrix = np.array([gc_vals, n_vals, qs_vals, len_vals]).T

        # Normalizacja kolumn (aby duża wartość długości nie zepsuła kolorów %GC czy %N)
        color_matrix = np.zeros_like(data_matrix, dtype=float)
        for j in range(4):
            cmax = data_matrix[:, j].max()
            color_matrix[:, j] = data_matrix[:, j] / cmax if cmax > 0 else 0

        cax = ax.imshow(color_matrix, cmap='YlGnBu', aspect='auto')
        ax.set_yticks(np.arange(len(labels)))
        ax.set_yticklabels(labels, fontsize=8)

        # Oś X przeniesiona na górę
        ax.set_xticks(np.arange(4))
        ax.set_xticklabels(['% GC', '% N', 'Quality Score', 'Długość'], fontsize=9)
        ax.xaxis.tick_top()

        font_size = 8 if len(subset) < 40 else 6

        # Rysowanie tekstu z prawdziwymi wartościami na znormalizowanym kolorze
        for i in range(len(labels)):
            for j in range(4):
                val = data_matrix[i, j]
                norm_val = color_matrix[i, j]
                text_color = 'white' if norm_val > 0.65 else 'black'

                # Formatowanie (dla długości liczba całkowita, reszta ułamki)
                text_str = f"{int(val)}" if j == 3 else f"{val:.1f}"

                ax.text(j, i, text_str, ha='center', va='center',
                        color=text_color, fontsize=font_size, weight='bold')

        # Usunięto colorbar, aby zrobić więcej miejsca dla 4 kolumn
        ax.set_title("Profil Parametrów QC\n\n", fontsize=10, fontweight='bold')
        self.heatmap_fig.tight_layout()
        self.heatmap_canvas.draw()

        self.heatmap_canvas.setMinimumSize(400, int(fig_height * 100))

    def export_excel(self):
        all_data = self.accepted_records + self.rejected_records
        if not all_data:
            return

        save_path, _ = QFileDialog.getSaveFileName(self, "Zapisz raport Excel", "Raport_Phage_QC.xlsx",
                                                   "Skoroszyt Excel (*.xlsx)")
        if not save_path:
            return

        try:
            import openpyxl
            from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
            from openpyxl.utils import get_column_letter

            wb = openpyxl.Workbook()
            ws = wb.active
            ws.title = "Analiza Jakości Phage QC"
            ws.views.sheetView[0].showGridLines = True

            headers = ["ID Sekwencji", "Długość (bp)", "GC (%)", "N (%)", "Quality Score", "Status",
                       "Komentarz / Plik Źródłowy"]
            ws.append(headers)

            header_font = Font(name="Segoe UI", size=11, bold=True, color="FFFFFF")
            header_fill = PatternFill(start_color="2C3E50", end_color="2C3E50", fill_type="solid")
            center_align = Alignment(horizontal="center", vertical="center")

            for col_idx, header in enumerate(headers, 1):
                cell = ws.cell(row=1, column=col_idx)
                cell.font = header_font
                cell.fill = header_fill
                cell.alignment = center_align

            fill_ok = PatternFill(start_color="E8F8F5", end_color="E8F8F5", fill_type="solid")
            fill_err = PatternFill(start_color="FDEDEC", end_color="FDEDEC", fill_type="solid")
            thin_border = Border(left=Side(style='thin', color='DDDDDD'), right=Side(style='thin', color='DDDDDD'),
                                 top=Side(style='thin', color='DDDDDD'), bottom=Side(style='thin', color='DDDDDD'))

            for row_idx, r in enumerate(all_data, 2):
                status = "Zaakceptowano" if r in self.accepted_records else "Odrzucono"
                comment = r["source"] if status == "Zaakceptowano" else r.get("reject_reason", "")

                row_vals = [r["id"], r["length"], round(r["gc_pct"], 2), round(r["n_pct"] * 100, 2), r["quality_score"],
                            status, comment]
                ws.append(row_vals)

                row_fill = fill_ok if status == "Zaakceptowano" else fill_err
                for col_idx in range(1, 8):
                    cell = ws.cell(row=row_idx, column=col_idx)
                    cell.fill = row_fill
                    cell.border = thin_border
                    if col_idx in [2, 3, 4, 5, 6]:
                        cell.alignment = center_align

            for col in ws.columns:
                max_len = max(len(str(cell.value or '')) for cell in col)
                col_letter = get_column_letter(col[0].column)
                ws.column_dimensions[col_letter].width = max(max_len + 3, 12)

            wb.save(save_path)
            QMessageBox.information(self, "Sukces", f"Zapisano pomyślnie raport Excel:\n{save_path}")
        except Exception as e:
            QMessageBox.critical(self, "Błąd", f"Nie udało się wygenerować raportu Excel: {e}")

    def export_pdf(self):
        all_data = self.accepted_records + self.rejected_records
        if not all_data:
            return

        save_path, _ = QFileDialog.getSaveFileName(self, "Zapisz raport PDF", "Raport_Phage_QC.pdf",
                                                   "Dokument PDF (*.pdf)")
        if not save_path:
            return

        try:
            from matplotlib.backends.backend_pdf import PdfPages
            import matplotlib.pyplot as plt

            with PdfPages(save_path) as pdf:
                fig, (ax_text, ax_chart) = plt.subplots(2, 1, figsize=(8.5, 11))

                fig.suptitle("Raport Kontroli Jakości Sekwencji (Phage QC)\n", fontsize=16, fontweight='bold',
                             color='#2C3E50')

                ax_text.axis('off')
                total_cnt = len(all_data)

                perc_acc = (len(self.accepted_records) / total_cnt * 100) if total_cnt > 0 else 0
                perc_rej = (len(self.rejected_records) / total_cnt * 100) if total_cnt > 0 else 0

                summary_txt = (
                    f"STATYSTYKI OGÓLNE ZALADOWANYCH SEKWENCJI:\n"
                    f"=========================================\n"
                    f"• Wszystkich odczytanych sekwencji: {total_cnt}\n"
                    f"• Zaakceptowane (Passed): {len(self.accepted_records)} ({perc_acc:.1f}%)\n"
                    f"• Odrzucone (Rejected): {len(self.rejected_records)} ({perc_rej:.1f}%)\n\n"
                    f"Kryteria progowe: Min Długość = {self.min_len_spin.value()} bp | Max N = {self.max_n_spin.value()}%"
                )
                ax_text.text(0.05, 0.3, summary_txt, fontsize=11, family='monospace',
                             bbox=dict(facecolor='#F8F9F9', alpha=0.9, boxstyle='round,pad=1', edgecolor='#BDC3C7'))

                subset = all_data[:25]
                labels = [r["id"][:12] for r in subset]
                # Dodano długość do macierzy danych PDF
                data_matrix = np.array(
                    [[r["gc_pct"], r["n_pct"] * 100, r["quality_score"], r["length"]] for r in subset])

                # Normalizacja, by kolory były spójne
                color_matrix = np.zeros_like(data_matrix, dtype=float)
                for j in range(4):
                    cmax = data_matrix[:, j].max()
                    color_matrix[:, j] = data_matrix[:, j] / cmax if cmax > 0 else 0

                cax = ax_chart.imshow(color_matrix, cmap='YlGnBu', aspect='auto')
                ax_chart.set_yticks(np.arange(len(labels)))
                ax_chart.set_yticklabels(labels, fontsize=8)

                # Oś X przeniesiona na górę w pliku PDF
                ax_chart.set_xticks(np.arange(4))
                ax_chart.set_xticklabels(['% GC', '% N', 'Quality Score', 'Długość'], fontsize=10)
                ax_chart.xaxis.tick_top()

                for i in range(len(labels)):
                    for j in range(4):
                        val = data_matrix[i, j]
                        norm_val = color_matrix[i, j]
                        text_color = 'white' if norm_val > 0.65 else 'black'
                        text_str = f"{int(val)}" if j == 3 else f"{val:.1f}"
                        ax_chart.text(j, i, text_str, ha='center', va='center', color=text_color, fontsize=8)

                ax_chart.set_title("Heatmapa profili jakości (Top 25 próbek)\n\n", fontsize=11, fontweight='bold')

                plt.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)

                chunk_size = 25
                for i in range(0, len(all_data), chunk_size):
                    chunk = all_data[i:i + chunk_size]
                    fig, ax = plt.subplots(figsize=(8.5, 11))
                    ax.axis('off')

                    table_data = []
                    cell_colors = []
                    for r in chunk:
                        status = "Zaakceptowano" if r in self.accepted_records else "Odrzucono"
                        table_data.append([r["id"][:15], r["length"], f"{r['gc_pct']:.1f}", f"{r['n_pct'] * 100:.1f}",
                                           str(r["quality_score"]), status])
                        cell_colors.append(['#E8F8F5'] * 6 if status == "Zaakceptowano" else ['#FDEDEC'] * 6)

                    pdf_table = ax.table(
                        cellText=table_data,
                        colLabels=["ID Sekwencji", "Długość", "GC (%)", "N (%)", "Quality", "Status"],
                        cellColours=cell_colors, loc='center'
                    )
                    pdf_table.auto_set_font_size(False)
                    pdf_table.set_fontsize(9)
                    pdf_table.scale(1, 1.4)

                    for col_idx in range(6):
                        cell = pdf_table[0, col_idx]
                        cell.set_text_props(weight='bold', color='white')
                        cell.set_facecolor('#2C3E50')

                    ax.set_title(f"Tabela Szczegółowa Wyników Analizy QC (Pozycje {i + 1} - {i + len(chunk)})\n",
                                 fontsize=11, fontweight='bold')
                    pdf.savefig(fig)
                    plt.close(fig)

            QMessageBox.information(self, "Sukces", f"Zapisano pomyślnie raport PDF:\n{save_path}")
        except Exception as e:
            QMessageBox.critical(self, "Błąd", f"Nie udało się wygenerować raportu PDF: {e}")

    def export_fasta(self):
        save_path, _ = QFileDialog.getSaveFileName(self, "Zapisz oczyszczony plik FASTA", "zwalidowane_fagi.fasta",
                                                   "FASTA Files (*.fasta);;All Files (*)")

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


# 4. URUCHOMIENIE APLIKACJI
if __name__ == '__main__':
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    window = PhageQCApp()
    window.show()
    sys.exit(app.exec())