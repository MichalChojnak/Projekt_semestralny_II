#!/usr/bin/env python3
# phage_qc.py
from __future__ import annotations

import os
import csv
import sys
import tempfile
import statistics
from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional
from collections import Counter

from Bio import SeqIO

import matplotlib.pyplot as plt
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4
from reportlab.lib.utils import ImageReader

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QTableWidget, QTableWidgetItem,
    QFileDialog, QHeaderView, QProgressBar
)

# =========================
# QC PARAMETERS
# =========================

MIN_LENGTH = 5000
MAX_N_PCT = 0.05
MAX_AMBIGUOUS_PCT = 0.02
MAX_HOMOPOLYMER = 30

VALID_BASES = {"A", "C", "G", "T", "N"}

# =========================
# DATA MODEL
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
    source_file: str = ""

# =========================
# UTILITIES
# =========================

def detect_format(path: str) -> str:
    ext = os.path.splitext(path)[1].lower()
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

def trim_ns_ends(seq: str) -> Tuple[str, int]:
    start = 0
    end = len(seq) - 1
    while start <= end and seq[start] == "N":
        start += 1
    while end >= start and seq[end] == "N":
        end -= 1
    if start > end:
        return "", len(seq)
    return seq[start:end+1], start + (len(seq) - 1 - end)

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
        if seq[i] == seq[i-1]:
            current += 1
        else:
            max_run = max(max_run, current)
            current = 1
    return max(max_run, current)

# =========================
# PROCESSING
# =========================

def load_records(path: str):
    fmt = detect_format(path)
    return list(SeqIO.parse(path, fmt))

def build_record(rec, source_file: str, do_trim: bool = True) -> Tuple[SequenceRecord, int]:
    raw_seq = str(rec.seq)
    corrections = 0
    # For FASTQ the SeqRecord still contains sequence only here; qualities ignored
    if do_trim:
        seq_trimmed, trimmed_bases = trim_ns_ends(raw_seq)
        if seq_trimmed == "":
            seq_trimmed = raw_seq  # fallback to original cleaned
    else:
        seq_trimmed = raw_seq
    cleaned, corr = clean_sequence(seq_trimmed)
    corrections += corr
    record = SequenceRecord(
        id=getattr(rec, "id", getattr(rec, "name", "unknown")),
        sequence=cleaned,
        length=len(cleaned),
        gc_pct=gc_content(cleaned),
        n_pct=n_content(cleaned),
        ambiguous_pct=ambiguous_content(cleaned),
        homopolymer_max=longest_homopolymer(cleaned),
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
    return True, "OK"

def process_file(path: str, do_trim: bool = True) -> Tuple[List[SequenceRecord], List[Tuple[SequenceRecord,str]], Dict]:
    recs = load_records(path)
    accepted = []
    rejected = []
    stats = {"total": 0, "corrections": 0, "format": detect_format(path)}
    for r in recs:
        rec, corr = build_record(r, path, do_trim)
        stats["corrections"] += corr
        stats["total"] += 1
        ok, reason = evaluate_qc(rec)
        if ok:
            accepted.append(rec)
        else:
            rejected.append((rec, reason))
    return accepted, rejected, stats

# =========================
# SUMMARY AND EXPORTS
# =========================

def compute_summary(accepted: List[SequenceRecord], rejected: List[Tuple[SequenceRecord,str]]) -> Dict:
    all_recs = accepted + [r for r,_ in rejected]
    lengths = [r.length for r in all_recs] if all_recs else []
    gc = [r.gc_pct for r in all_recs] if all_recs else []
    summary = {
        "total": len(all_recs),
        "accepted": len(accepted),
        "rejected": len(rejected),
        "length_mean": round(statistics.mean(lengths),2) if lengths else 0,
        "length_median": round(statistics.median(lengths),2) if lengths else 0,
        "length_min": min(lengths) if lengths else 0,
        "length_max": max(lengths) if lengths else 0,
        "gc_mean": round(statistics.mean(gc),2) if gc else 0
    }
    reasons = Counter(reason for _,reason in rejected)
    summary["rejection_reasons"] = dict(reasons)
    return summary

def write_sequence_csv(accepted: List[SequenceRecord], rejected: List[Tuple[SequenceRecord,str]], path: str):
    rows = []
    for r in accepted:
        rows.append({
            "id": r.id,
            "length": r.length,
            "gc_pct": f"{r.gc_pct:.2f}",
            "n_pct": f"{r.n_pct:.4f}",
            "status": "ACCEPTED",
            "reason": ""
        })
    for r, reason in rejected:
        rows.append({
            "id": r.id,
            "length": r.length,
            "gc_pct": f"{r.gc_pct:.2f}",
            "n_pct": f"{r.n_pct:.4f}",
            "status": "REJECTED",
            "reason": reason
        })
    if not rows:
        raise ValueError("No records to write")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

def write_fasta(records: List[SequenceRecord], path: str):
    with open(path, "w", encoding="utf-8") as f:
        for r in records:
            f.write(f">{r.id}\n")
            seq = r.sequence
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")

def plot_histograms(accepted: List[SequenceRecord], rejected: List[Tuple[SequenceRecord,str]], out_dir: str) -> Dict[str,str]:
    all_recs = accepted + [r for r,_ in rejected]
    lengths = [r.length for r in all_recs] if all_recs else []
    gc = [r.gc_pct for r in all_recs] if all_recs else []
    os.makedirs(out_dir, exist_ok=True)
    imgs = {}
    if lengths:
        plt.figure(figsize=(6,3.5))
        plt.hist(lengths, bins=50, color="#4C72B0", edgecolor="black")
        plt.xlabel("Length bp")
        plt.ylabel("Count")
        plt.title("Length distribution")
        p = os.path.join(out_dir, "hist_length.png")
        plt.tight_layout()
        plt.savefig(p, dpi=150)
        plt.close()
        imgs["length"] = p
    if gc:
        plt.figure(figsize=(6,3.5))
        plt.hist(gc, bins=50, color="#55A868", edgecolor="black")
        plt.xlabel("GC%")
        plt.ylabel("Count")
        plt.title("GC distribution")
        p = os.path.join(out_dir, "hist_gc.png")
        plt.tight_layout()
        plt.savefig(p, dpi=150)
        plt.close()
        imgs["gc"] = p
    return imgs

def write_pdf(summary: Dict, accepted: List[SequenceRecord], rejected: List[Tuple[SequenceRecord,str]], path: str):
    tmp = tempfile.mkdtemp(prefix="phageqc_")
    imgs = plot_histograms(accepted, rejected, tmp)
    c = canvas.Canvas(path, pagesize=A4)
    w, h = A4
    margin = 50
    y = h - margin
    c.setFont("Helvetica-Bold", 16)
    c.drawString(margin, y, "Phage QC Summary")
    y -= 28
    c.setFont("Helvetica", 10)
    c.drawString(margin, y, f"Total sequences: {summary.get('total',0)}")
    y -= 14
    c.drawString(margin, y, f"Accepted: {summary.get('accepted',0)}")
    y -= 14
    c.drawString(margin, y, f"Rejected: {summary.get('rejected',0)}")
    y -= 18
    for k in ("length_mean","length_median","length_min","length_max","gc_mean"):
        c.drawString(margin, y, f"{k}: {summary.get(k,'')}")
        y -= 14
        if y < margin + 120:
            c.showPage()
            y = h - margin
            c.setFont("Helvetica", 10)
    if imgs.get("length"):
        try:
            img = ImageReader(imgs["length"])
            c.drawImage(img, margin, y-180, width=500, height=160, preserveAspectRatio=True)
            y -= 180 + 10
        except Exception:
            pass
    if imgs.get("gc"):
        if y < margin + 180:
            c.showPage()
            y = h - margin
        try:
            img = ImageReader(imgs["gc"])
            c.drawImage(img, margin, y-180, width=500, height=160, preserveAspectRatio=True)
            y -= 180 + 10
        except Exception:
            pass
    if y < margin + 120:
        c.showPage()
        y = h - margin
    c.setFont("Helvetica-Bold", 12)
    c.drawString(margin, y, "Top rejection reasons")
    y -= 16
    c.setFont("Helvetica", 10)
    for reason, count in summary.get("rejection_reasons", {}).items():
        c.drawString(margin, y, f"{reason}: {count}")
        y -= 12
        if y < margin + 60:
            c.showPage()
            y = h - margin
            c.setFont("Helvetica", 10)
    if y < margin + 120:
        c.showPage()
        y = h - margin
    c.setFont("Helvetica-Bold", 12)
    c.drawString(margin, y, "Sample records")
    y -= 16
    c.setFont("Helvetica", 9)
    sample = (accepted + [r for r,_ in rejected])[:80]
    for r in sample:
        line = f"{r.id} | len={r.length} | GC={r.gc_pct:.2f}% | N={r.n_pct:.4f} | status={'ACCEPTED' if r in accepted else 'REJECTED'}"
        c.drawString(margin, y, line)
        y -= 12
        if y < margin + 40:
            c.showPage()
            y = h - margin
            c.setFont("Helvetica", 9)
    c.showPage()
    c.save()
    # cleanup
    try:
        for p in imgs.values():
            if os.path.exists(p):
                os.remove(p)
        os.rmdir(tmp)
    except Exception:
        pass

# =========================
# GUI
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
        self.bar.setMaximum(max(total,1))
        self.bar.setValue(value)
        self.lbl.setText(text)

class PhageQCApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Phage QC Tool - FASTA Normalizer")
        self.resize(1100,700)
        self.accepted: List[SequenceRecord] = []
        self.rejected: List[Tuple[SequenceRecord,str]] = []
        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)
        stats_layout = QHBoxLayout()
        self.sb_total = StatBar("Total", "#58a6ff")
        self.sb_acc = StatBar("Accepted", "#50e3c2")
        self.sb_rej = StatBar("Rejected", "#ff7b72")
        stats_layout.addWidget(self.sb_total)
        stats_layout.addWidget(self.sb_acc)
        stats_layout.addWidget(self.sb_rej)
        layout.addLayout(stats_layout)
        btn_layout = QHBoxLayout()
        self.btn_load = QPushButton("Load FASTA files")
        self.btn_csv = QPushButton("Export CSV")
        self.btn_pdf = QPushButton("Export PDF")
        self.btn_fasta = QPushButton("Export merged FASTA")
        btn_layout.addWidget(self.btn_load)
        btn_layout.addWidget(self.btn_csv)
        btn_layout.addWidget(self.btn_pdf)
        btn_layout.addWidget(self.btn_fasta)
        layout.addLayout(btn_layout)
        self.table = QTableWidget()
        self.table.setColumnCount(6)
        self.table.setHorizontalHeaderLabels(["ID","Length","GC%","N%","Homopolymer","Status"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        layout.addWidget(self.table)
        self.btn_load.clicked.connect(self.load_files)
        self.btn_csv.clicked.connect(self.export_csv)
        self.btn_pdf.clicked.connect(self.export_pdf)
        self.btn_fasta.clicked.connect(self.export_fasta)

    def load_files(self):
        files, _ = QFileDialog.getOpenFileNames(
            self,
            "Select FASTA or FASTQ files",
            "",
            "Seq Files (*.fasta *.fa *.fna *.fastq *.fq);;All Files (*)"
        )
        if not files:
            return
        self.accepted.clear()
        self.rejected.clear()
        total = len(files)
        processed = 0
        for f in files:
            try:
                a, r, stats = process_file(f, do_trim=True)
            except Exception as e:
                print(f"Error processing {f}: {e}")
                continue
            self.accepted.extend(a)
            self.rejected.extend(r)
            processed += 1
            self.sb_total.update(processed, total, f"Files: {processed}/{total}")
        self.update_stats()
        self.update_table()

    def update_stats(self):
        total = len(self.accepted) + len(self.rejected)
        acc = len(self.accepted)
        rej = len(self.rejected)
        self.sb_total.update(total, max(total,1), f"Total: {total}")
        self.sb_acc.update(acc, max(total,1), f"Accepted: {acc}")
        self.sb_rej.update(rej, max(total,1), f"Rejected: {rej}")
        self.summary = compute_summary(self.accepted, self.rejected)

    def update_table(self):
        rows = self.accepted + [r for r,_ in self.rejected]
        self.table.setRowCount(len(rows))
        for i, r in enumerate(rows):
            status = "ACCEPTED" if r in self.accepted else "REJECTED"
            self.table.setItem(i,0,QTableWidgetItem(r.id))
            self.table.setItem(i,1,QTableWidgetItem(str(r.length)))
            self.table.setItem(i,2,QTableWidgetItem(f"{r.gc_pct:.2f}"))
            self.table.setItem(i,3,QTableWidgetItem(f"{r.n_pct:.4f}"))
            self.table.setItem(i,4,QTableWidgetItem(str(r.homopolymer_max)))
            self.table.setItem(i,5,QTableWidgetItem(status))

    def export_csv(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save CSV", "", "CSV Files (*.csv)")
        if not path:
            return
        write_sequence_csv(self.accepted, self.rejected, path)

    def export_pdf(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save PDF", "", "PDF Files (*.pdf)")
        if not path:
            return
        write_pdf(self.summary, self.accepted, self.rejected, path)

    def export_fasta(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save merged FASTA (accepted)", "", "FASTA Files (*.fasta *.fa)")
        if not path:
            return
        accepted_path = path if path.endswith(".fasta") or path.endswith(".fa") else path + ".fasta"
        rejected_path = os.path.splitext(accepted_path)[0] + "_rejected.fasta"
        write_fasta(self.accepted, accepted_path)
        write_fasta([r for r,_ in self.rejected], rejected_path)

# =========================
# MAIN
# =========================

def main():
    app = QApplication(sys.argv)
    win = PhageQCApp()
    win.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
