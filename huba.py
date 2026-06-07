#!/usr/bin/env python3
# phage_qc.py
from __future__ import annotations

import os
import csv
import sys
import math
import statistics
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
from collections import Counter
from pathlib import Path

from Bio import SeqIO

# Optional plotting (used for small histogram pages if desired)
import matplotlib.pyplot as plt

# PDF generation
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4

# PyQt6 GUI
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QTableWidget, QTableWidgetItem,
    QFileDialog, QHeaderView, QProgressBar
)
from PyQt6.QtCore import Qt

# =========================
# 🔬 QC PARAMETERS
# =========================

MIN_LENGTH = 5000
MAX_N_PCT = 0.05
MAX_AMBIGUOUS_PCT = 0.02
MIN_PHRED = 20
MAX_HOMOPOLYMER = 30

VALID_BASES = {"A", "C", "G", "T", "N"}


# =========================
# 📦 DATA MODEL
# =========================

@dataclass
class SequenceRecord:
    id: str
    sequence: str
    length: int
    gc_pct: float
    n_pct: float
    ambiguous_pct: float
    homopolymer_max: int
    mean_phred: Optional[float] = None
    source_file: str = ""


# =========================
# 🧬 CORE FUNCTIONS
# =========================

def clean_sequence(seq: str) -> Tuple[str, int]:
    seq = seq.upper()
    cleaned = []
    corrections = 0

    for b in seq:
        if b not in VALID_BASES:
            cleaned.append("N")
            corrections += 1
        else:
            cleaned.append(b)

    return "".join(cleaned), corrections


def gc_content(seq: str) -> float:
    return (seq.count("G") + seq.count("C")) / len(seq) * 100 if seq else 0.0


def n_content(seq: str) -> float:
    return seq.count("N") / len(seq) if seq else 0.0


def ambiguous_content(seq: str) -> float:
    return sum(1 for b in seq if b not in VALID_BASES) / len(seq) if seq else 0.0


def longest_homopolymer(seq: str) -> int:
    if not seq:
        return 0

    max_run = 1
    current = 1

    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            current += 1
        else:
            max_run = max(max_run, current)
            current = 1

    return max(max_run, current)


def mean_phred_quality(record) -> Optional[float]:
    quals = getattr(record, "letter_annotations", {}).get("phred_quality", None)
    if not quals:
        return None
    return sum(quals) / len(quals)


# =========================
# 📥 FILE LOADING
# =========================

def detect_format(file_path: str) -> str:
    ext = os.path.splitext(file_path)[1].lower()

    return {
        ".fasta": "fasta",
        ".fa": "fasta",
        ".fna": "fasta",
        ".fastq": "fastq",
        ".fq": "fastq",
        ".gb": "genbank",
        ".gbk": "genbank",
        ".embl": "embl",
    }.get(ext, "fasta")


def load_sequences(file_path: str):
    fmt = detect_format(file_path)
    return list(SeqIO.parse(file_path, fmt))


def build_record(rec, source_file: str) -> Tuple[SequenceRecord, int]:
    raw_seq = str(rec.seq)
    cleaned_seq, corrections = clean_sequence(raw_seq)

    record = SequenceRecord(
        id=getattr(rec, "id", getattr(rec, "name", "unknown")),
        sequence=cleaned_seq,
        length=len(cleaned_seq),
        gc_pct=gc_content(cleaned_seq),
        n_pct=n_content(cleaned_seq),
        ambiguous_pct=ambiguous_content(cleaned_seq),
        homopolymer_max=longest_homopolymer(cleaned_seq),
        mean_phred=mean_phred_quality(rec),
        source_file=source_file
    )

    return record, corrections


def evaluate_qc(record: SequenceRecord) -> Tuple[bool, str]:
    if record.length < MIN_LENGTH:
        return False, f"TOO_SHORT ({record.length})"

    if record.n_pct > MAX_N_PCT:
        return False, f"HIGH_N ({record.n_pct:.2%})"

    if record.ambiguous_pct > MAX_AMBIGUOUS_PCT:
        return False, f"AMBIGUOUS ({record.ambiguous_pct:.2%})"

    if record.homopolymer_max > MAX_HOMOPOLYMER:
        return False, f"HOMOPOLYMER ({record.homopolymer_max})"

    if record.mean_phred is not None and record.mean_phred < MIN_PHRED:
        return False, f"LOW_PHRED ({record.mean_phred:.1f})"

    return True, "OK"


def process_file(file_path: str) -> Tuple[List[SequenceRecord], List[Tuple[SequenceRecord, str]], Dict]:
    raw_records = load_sequences(file_path)

    accepted: List[SequenceRecord] = []
    rejected: List[Tuple[SequenceRecord, str]] = []

    stats = {
        "total": 0,
        "corrections": 0,
        "format": detect_format(file_path)
    }

    for r in raw_records:
        record, corrections = build_record(r, file_path)
        stats["corrections"] += corrections
        stats["total"] += 1

        ok, reason = evaluate_qc(record)
        if ok:
            accepted.append(record)
        else:
            rejected.append((record, reason))

    return accepted, rejected, stats


# =========================
# 📊 SUMMARY & EXPORT HELPERS
# =========================

def compute_summary(accepted: List[SequenceRecord], rejected: List[Tuple[SequenceRecord, str]]) -> Dict:
    all_records = accepted + [r for r, _ in rejected]
    lengths = [r.length for r in all_records] if all_records else []
    gc = [r.gc_pct for r in all_records] if all_records else []

    summary = {
        "total": len(all_records),
        "accepted": len(accepted),
        "rejected": len(rejected),
        "length_mean": statistics.mean(lengths) if lengths else 0,
        "length_median": statistics.median(lengths) if lengths else 0,
        "length_min": min(lengths) if lengths else 0,
        "length_max": max(lengths) if lengths else 0,
        "gc_mean": statistics.mean(gc) if gc else 0,
    }

    # breakdown of rejection reasons
    reasons = Counter(reason for _, reason in rejected)
    summary["rejection_reasons"] = dict(reasons)

    return summary


def write_sequence_report(accepted: List[SequenceRecord], rejected: List[Tuple[SequenceRecord, str]], path: str):
    rows = []
    for r in accepted:
        rows.append({
            "id": r.id,
            "length": r.length,
            "gc_pct": f"{r.gc_pct:.2f}",
            "n_pct": f"{r.n_pct:.4f}",
            "ambiguous_pct": f"{r.ambiguous_pct:.4f}",
            "homopolymer": r.homopolymer_max,
            "mean_phred": f"{r.mean_phred:.1f}" if r.mean_phred is not None else "",
            "status": "ACCEPTED",
            "reason": ""
        })

    for r, reason in rejected:
        rows.append({
            "id": r.id,
            "length": r.length,
            "gc_pct": f"{r.gc_pct:.2f}",
            "n_pct": f"{r.n_pct:.4f}",
            "ambiguous_pct": f"{r.ambiguous_pct:.4f}",
            "homopolymer": r.homopolymer_max,
            "mean_phred": f"{r.mean_phred:.1f}" if r.mean_phred is not None else "",
            "status": "REJECTED",
            "reason": reason
        })

    if not rows:
        raise ValueError("No records to write to CSV.")

    fieldnames = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def write_fasta(records: List[SequenceRecord], path: str):
    # simple FASTA writer
    with open(path, "w", encoding="utf-8") as f:
        for r in records:
            f.write(f">{r.id}\n")
            # wrap sequence at 80 chars
            seq = r.sequence
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")


def write_pdf_report(summary: Dict, accepted: List[SequenceRecord], rejected: List[Tuple[SequenceRecord, str]], path: str):
    c = canvas.Canvas(path, pagesize=A4)
    w, h = A4
    margin = 50
    y = h - margin

    c.setFont("Helvetica-Bold", 16)
    c.drawString(margin, y, "Phage QC Report")
    y -= 30

    c.setFont("Helvetica", 10)
    # Basic summary
    for k in ("total", "accepted", "rejected", "length_mean", "length_median", "length_min", "length_max", "gc_mean"):
        v = summary.get(k, "")
        c.drawString(margin, y, f"{k}: {v}")
        y -= 14
        if y < margin + 100:
            c.showPage()
            y = h - margin
            c.setFont("Helvetica", 10)

    # Rejection reasons
    y -= 6
    c.setFont("Helvetica-Bold", 12)
    c.drawString(margin, y, "Rejection reasons")
    y -= 14
    c.setFont("Helvetica", 10)
    for reason, count in summary.get("rejection_reasons", {}).items():
        c.drawString(margin, y, f"{reason}: {count}")
        y -= 12
        if y < margin + 60:
            c.showPage()
            y = h - margin
            c.setFont("Helvetica", 10)

    # Optionally add small table of first N records
    y -= 10
    c.setFont("Helvetica-Bold", 12)
    c.drawString(margin, y, "Sample records")
    y -= 16
    c.setFont("Helvetica", 9)
    sample = (accepted + [r for r, _ in rejected])[:50]
    for r in sample:
        line = f"{r.id} | len={r.length} | GC={r.gc_pct:.2f}% | N={r.n_pct:.4f}"
        c.drawString(margin, y, line)
        y -= 12
        if y < margin + 40:
            c.showPage()
            y = h - margin
            c.setFont("Helvetica", 9)

    c.showPage()
    c.save()


# =========================
# 🎛 UI COMPONENT
# =========================

class StatBar(QWidget):
    def __init__(self, label, color="#50E3C2"):
        super().__init__()

        layout = QVBoxLayout(self)

        self.lbl = QLabel(f"{label}: 0")
        self.bar = QProgressBar()

        self.bar.setTextVisible(False)
        self.bar.setFixedHeight(8)

        self.bar.setStyleSheet(f"""
            QProgressBar {{
                background-color:#1a1e24;
                border-radius:4px;
            }}
            QProgressBar::chunk {{
                background-color:{color};
                border-radius:4px;
            }}
        """)

        layout.addWidget(self.lbl)
        layout.addWidget(self.bar)

    def update(self, value, total, text):
        self.bar.setMaximum(max(total, 1))
        self.bar.setValue(value)
        self.lbl.setText(text)


# =========================
# 🧬 MAIN APP
# =========================

class PhageQCApp(QMainWindow):

    def __init__(self):
        super().__init__()

        self.setWindowTitle("Phage QC Tool")
        self.resize(1100, 700)

        self.accepted: List[SequenceRecord] = []
        self.rejected: List[Tuple[SequenceRecord, str]] = []
        self.stats: Dict = {}

        central = QWidget()
        self.setCentralWidget(central)

        layout = QVBoxLayout(central)

        # STAT BARS
        stats_layout = QHBoxLayout()
        self.sb_total = StatBar("Total", "#58a6ff")
        self.sb_acc = StatBar("Accepted", "#50e3c2")
        self.sb_rej = StatBar("Rejected", "#ff7b72")
        stats_layout.addWidget(self.sb_total)
        stats_layout.addWidget(self.sb_acc)
        stats_layout.addWidget(self.sb_rej)
        layout.addLayout(stats_layout)

        # BUTTONS
        btn_layout = QHBoxLayout()
        self.btn_load = QPushButton("Load Files")
        self.btn_csv = QPushButton("Export CSV")
        self.btn_pdf = QPushButton("Export PDF")
        self.btn_fasta = QPushButton("Export FASTA")
        btn_layout.addWidget(self.btn_load)
        btn_layout.addWidget(self.btn_csv)
        btn_layout.addWidget(self.btn_pdf)
        btn_layout.addWidget(self.btn_fasta)
        layout.addLayout(btn_layout)

        # TABLE
        self.table = QTableWidget()
        self.table.setColumnCount(7)
        self.table.setHorizontalHeaderLabels([
            "ID", "Length", "GC%", "N%", "Ambiguous%", "Homopolymer", "Status"
        ])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        layout.addWidget(self.table)

        # SIGNALS
        self.btn_load.clicked.connect(self.load_files)
        self.btn_csv.clicked.connect(self.export_csv)
        self.btn_pdf.clicked.connect(self.export_pdf)
        self.btn_fasta.clicked.connect(self.export_fasta)

    # LOAD FILES
    def load_files(self):

        files, _ = QFileDialog.getOpenFileNames(
            self,
            "Select sequencing files",
            "",
            "Seq Files (*.fasta *.fa *.fna *.fastq *.fq *.gb *.gbk *.embl)"
        )

        if not files:
            return

        self.accepted.clear()
        self.rejected.clear()

        total_files = len(files)
        processed = 0

        for f in files:
            records, rej, stats = process_file(f)
            self.accepted.extend(records)
            self.rejected.extend(rej)
            processed += 1
            self.sb_total.update(processed, total_files, f"Files: {processed}/{total_files}")

        self.update_stats()
        self.update_table()

    # UPDATE STATS
    def update_stats(self):
        total = len(self.accepted) + len(self.rejected)
        acc = len(self.accepted)
        rej = len(self.rejected)

        self.sb_total.update(total, total, f"Total: {total}")
        self.sb_acc.update(acc, total, f"Accepted: {acc}")
        self.sb_rej.update(rej, total, f"Rejected: {rej}")

        self.summary = compute_summary(self.accepted, self.rejected)

    # TABLE
    def update_table(self):
        rows = self.accepted + [r for r, _ in self.rejected]
        self.table.setRowCount(len(rows))

        for i, r in enumerate(rows):
            status = "ACCEPTED" if r in self.accepted else "REJECTED"
            self.table.setItem(i, 0, QTableWidgetItem(r.id))
            self.table.setItem(i, 1, QTableWidgetItem(str(r.length)))
            self.table.setItem(i, 2, QTableWidgetItem(f"{r.gc_pct:.2f}"))
            self.table.setItem(i, 3, QTableWidgetItem(f"{r.n_pct:.4f}"))
            self.table.setItem(i, 4, QTableWidgetItem(f"{r.ambiguous_pct:.4f}"))
            self.table.setItem(i, 5, QTableWidgetItem(str(r.homopolymer_max)))
            self.table.setItem(i, 6, QTableWidgetItem(status))

    # EXPORT CSV
    def export_csv(self):
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save CSV",
            "",
            "CSV Files (*.csv)"
        )
        if not path:
            return
        write_sequence_report(self.accepted, self.rejected, path)

    # EXPORT PDF
    def export_pdf(self):
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save PDF",
            "",
            "PDF Files (*.pdf)"
        )
        if not path:
            return
        write_pdf_report(self.summary, self.accepted, self.rejected, path)

    # EXPORT FASTA (accepted / rejected)
    def export_fasta(self):
        base_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save accepted FASTA (base name)",
            "",
            "FASTA Files (*.fasta *.fa)"
        )
        if not base_path:
            return
        # write accepted and rejected with suffixes
        accepted_path = base_path if base_path.endswith(".fasta") or base_path.endswith(".fa") else base_path + "_accepted.fasta"
        rejected_path = base_path + "_rejected.fasta"
        write_fasta(self.accepted, accepted_path)
        write_fasta([r for r, _ in self.rejected], rejected_path)


# =========================
# 🚀 MAIN
# =========================

def main():
    app = QApplication(sys.argv)
    win = PhageQCApp()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
