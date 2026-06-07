#!/usr/bin/env python3
# phage_qc.py
from __future__ import annotations

import os
import csv
import sys
import zipfile
import tempfile
import statistics
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
from collections import Counter

from Bio import SeqIO

# plotting
import matplotlib.pyplot as plt

# PDF generation
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4
from reportlab.lib.utils import ImageReader

# PyQt6 GUI
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QTableWidget, QTableWidgetItem,
    QFileDialog, QHeaderView, QProgressBar
)

# =========================
# QC PARAMETERS (możesz zmienić)
# =========================

MIN_LENGTH = 5000
MAX_N_PCT = 0.05
MAX_AMBIGUOUS_PCT = 0.02
MIN_PHRED = 20
MAX_HOMOPOLYMER = 30

# trimming params
TRIM_QUALITY_THRESHOLD = 20   # phred threshold for end trimming
TRIM_MIN_LENGTH = 50          # minimal length after trimming to keep sequence

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
    mean_phred: Optional[float] = None
    source_file: str = ""
    fastqc_info: Optional[Dict] = None

# =========================
# CORE: cleaning, trimming, metrics
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

# trimming: for FASTQ use quality; for FASTA trim Ns at ends
def trim_by_quality(record) -> Tuple[str, int]:
    """
    Trim low-quality ends based on TRIM_QUALITY_THRESHOLD.
    Returns (trimmed_sequence, number_of_bases_trimmed)
    Works only if record has letter_annotations phred_quality.
    """
    quals = getattr(record, "letter_annotations", {}).get("phred_quality", None)
    seq = str(record.seq)
    if not quals:
        # fallback: trim Ns at ends
        return trim_ns_ends(seq)
    # trim from left
    left = 0
    right = len(quals) - 1
    while left <= right and quals[left] < TRIM_QUALITY_THRESHOLD:
        left += 1
    while right >= left and quals[right] < TRIM_QUALITY_THRESHOLD:
        right -= 1
    if left > right:
        return "", len(seq)
    trimmed = seq[left:right+1]
    trimmed, corrections = clean_sequence(trimmed)
    return trimmed, left + (len(seq) - 1 - right)

def trim_ns_ends(seq: str) -> Tuple[str, int]:
    start = 0
    end = len(seq) - 1
    while start <= end and seq[start] == "N":
        start += 1
    while end >= start and seq[end] == "N":
        end -= 1
    if start > end:
        return "", len(seq)
    trimmed = seq[start:end+1]
    return trimmed, start + (len(seq) - 1 - end)

# =========================
# FASTQC parsing (zip or folder)
# =========================

def parse_fastqc_from_zip(zip_path: str) -> Dict:
    info = {}
    try:
        with zipfile.ZipFile(zip_path, 'r') as z:
            candidates = [n for n in z.namelist() if n.endswith("fastqc_data.txt")]
            if not candidates:
                return info
            name = candidates[0]
            with z.open(name) as fh:
                text = fh.read().decode('utf-8', errors='ignore')
                info = parse_fastqc_text(text)
    except Exception:
        return {}
    return info

def parse_fastqc_text(text: str) -> Dict:
    info = {}
    lines = text.splitlines()
    section = None
    basic = {}
    module_status = {}
    for ln in lines:
        if ln.startswith(">>"):
            parts = ln.split("\t")
            header = parts[0][2:]
            status = parts[1] if len(parts) > 1 else ""
            if header == "END_MODULE":
                section = None
            else:
                section = header
                module_status[header] = status
            continue
        if section == "Basic Statistics":
            if "\t" in ln:
                k, v = ln.split("\t", 1)
                basic[k.strip()] = v.strip()
    if basic:
        try:
            info["Total Sequences"] = int(basic.get("Total Sequences", "0"))
        except ValueError:
            info["Total Sequences"] = 0
        info["Sequence length"] = basic.get("Sequence length", "")
        try:
            info["%GC"] = float(basic.get("%GC", "0"))
        except ValueError:
            info["%GC"] = 0.0
    info["modules"] = module_status
    return info

def parse_fastqc(path: str) -> Dict:
    if os.path.isdir(path):
        for root, dirs, files in os.walk(path):
            for f in files:
                if f == "fastqc_data.txt":
                    with open(os.path.join(root, f), "r", encoding="utf-8", errors="ignore") as fh:
                        return parse_fastqc_text(fh.read())
        return {}
    elif zipfile.is_zipfile(path):
        return parse_fastqc_from_zip(path)
    else:
        if os.path.basename(path) == "fastqc_data.txt" and os.path.exists(path):
            with open(path, "r", encoding="utf-8", errors="ignore") as fh:
                return parse_fastqc_text(fh.read())
    return {}

# =========================
# FILE LOADING & PROCESSING
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

def build_record(rec, source_file: str, do_trim: bool = True) -> Tuple[SequenceRecord, int]:
    """
    Build SequenceRecord and optionally trim ends.
    Returns (SequenceRecord, corrections_count)
    """
    corrections = 0
    seq_raw = str(rec.seq)
    if do_trim:
        # try quality trimming first (for FASTQ)
        trimmed_seq, trimmed_bases = trim_by_quality(rec)
        if trimmed_seq == "" and trimmed_bases > 0:
            # fully trimmed away; fallback to cleaning original
            cleaned_seq, corr = clean_sequence(seq_raw)
            corrections += corr
            seq = cleaned_seq
        else:
            seq = trimmed_seq if trimmed_seq else seq_raw
    else:
        seq = seq_raw
    cleaned_seq, corr = clean_sequence(seq)
    corrections += corr
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

def process_file(file_path: str, do_trim: bool = True) -> Tuple[List[SequenceRecord], List[Tuple[SequenceRecord, str]], Dict]:
    raw_records = load_sequences(file_path)
    accepted: List[SequenceRecord] = []
    rejected: List[Tuple[SequenceRecord, str]] = []
    stats = {"total": 0, "corrections": 0, "format": detect_format(file_path)}
    for r in raw_records:
        rec, corrections = build_record(r, file_path, do_trim)
        stats["corrections"] += corrections
        stats["total"] += 1
        ok, reason = evaluate_qc(rec)
        if ok:
            accepted.append(rec)
        else:
            rejected.append((rec, reason))
    return accepted, rejected, stats

# =========================
# SUMMARY, HISTOGRAMS, EXPORTS
# =========================

def compute_summary(accepted: List[SequenceRecord], rejected: List[Tuple[SequenceRecord, str]], fastqc_reports: Dict[str, Dict]) -> Dict:
    all_records = accepted + [r for r, _ in rejected]
    lengths = [r.length for r in all_records] if all_records else []
    gc = [r.gc_pct for r in all_records] if all_records else []
    summary = {
        "total": len(all_records),
        "accepted": len(accepted),
        "rejected": len(rejected),
        "length_mean": round(statistics.mean(lengths), 2) if lengths else 0,
        "length_median": round(statistics.median(lengths), 2) if lengths else 0,
        "length_min": min(lengths) if lengths else 0,
        "length_max": max(lengths) if lengths else 0,
        "gc_mean": round(statistics.mean(gc), 2) if gc else 0,
    }
    reasons = Counter(reason for _, reason in rejected)
    summary["rejection_reasons"] = dict(reasons)
    summary["fastqc_reports"] = fastqc_reports
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
    with open(path, "w", encoding="utf-8") as f:
        for r in records:
            f.write(f">{r.id}\n")
            seq = r.sequence
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")

def plot_histograms(accepted: List[SequenceRecord], rejected: List[Tuple[SequenceRecord, str]], out_dir: str) -> Dict[str, str]:
    """
    Create histograms for lengths and GC and save PNGs.
    Returns dict with paths.
    """
    all_records = accepted + [r for r, _ in rejected]
    lengths = [r.length for r in all_records] if all_records else []
    gc = [r.gc_pct for r in all_records] if all_records else []

    os.makedirs(out_dir, exist_ok=True)
    imgs = {}

    if lengths:
        plt.figure(figsize=(6,3.5))
        plt.hist(lengths, bins=50, color="#4C72B0", edgecolor="black")
        plt.xlabel("Sequence length (bp)")
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

def write_pdf_report(summary: Dict, accepted: List[SequenceRecord], rejected: List[Tuple[SequenceRecord, str]], path: str):
    # create histograms in temp dir
    tmpdir = tempfile.mkdtemp(prefix="phageqc_")
    imgs = plot_histograms(accepted, rejected, tmpdir)

    c = canvas.Canvas(path, pagesize=A4)
    w, h = A4
    margin = 50
    y = h - margin

    c.setFont("Helvetica-Bold", 16)
    c.drawString(margin, y, "Phage QC Report")
    y -= 28

    c.setFont("Helvetica", 10)
    c.drawString(margin, y, f"Total sequences: {summary.get('total',0)}")
    y -= 14
    c.drawString(margin, y, f"Accepted: {summary.get('accepted',0)}")
    y -= 14
    c.drawString(margin, y, f"Rejected: {summary.get('rejected',0)}")
    y -= 18

    # summary metrics
    metrics = [
        ("Length mean", summary.get("length_mean", "")),
        ("Length median", summary.get("length_median", "")),
        ("Length min", summary.get("length_min", "")),
        ("Length max", summary.get("length_max", "")),
        ("GC mean", summary.get("gc_mean", "")),
    ]
    for k, v in metrics:
        c.drawString(margin, y, f"{k}: {v}")
        y -= 14
        if y < margin + 120:
            c.showPage()
            y = h - margin
            c.setFont("Helvetica", 10)

    # insert histograms
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

    # rejection reasons
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

    # FastQC mini-section
    if y < margin + 120:
        c.showPage()
        y = h - margin
    c.setFont("Helvetica-Bold", 12)
    c.drawString(margin, y, "FastQC reports (if provided)")
    y -= 14
    c.setFont("Helvetica", 9)
    for fname, frep in summary.get("fastqc_reports", {}).items():
        c.drawString(margin, y, f"{fname}: TotalSeq={frep.get('Total Sequences', 'NA')}  %GC={frep.get('%GC', 'NA')}")
        y -= 12
        mods = frep.get("modules", {})
        if mods:
            short = ", ".join(f"{m}:{s}" for m, s in list(mods.items())[:6])
            c.drawString(margin + 10, y, short)
            y -= 12
        if y < margin + 60:
            c.showPage()
            y = h - margin
            c.setFont("Helvetica", 9)

    # sample records table (first N)
    if y < margin + 120:
        c.showPage()
        y = h - margin
    c.setFont("Helvetica-Bold", 12)
    c.drawString(margin, y, "Sample records")
    y -= 16
    c.setFont("Helvetica", 9)
    sample = (accepted + [r for r, _ in rejected])[:80]
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

    # cleanup temp images
    try:
        for p in imgs.values():
            if os.path.exists(p):
                os.remove(p)
        os.rmdir(tmpdir)
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
        self.bar.setMaximum(max(total, 1))
        self.bar.setValue(value)
        self.lbl.setText(text)

class PhageQCApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Phage QC Tool")
        self.resize(1200, 750)
        self.accepted: List[SequenceRecord] = []
        self.rejected: List[Tuple[SequenceRecord, str]] = []
        self.fastqc_reports: Dict[str, Dict] = {}
        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)
        # stat bars
        stats_layout = QHBoxLayout()
        self.sb_total = StatBar("Total", "#58a6ff")
        self.sb_acc = StatBar("Accepted", "#50e3c2")
        self.sb_rej = StatBar("Rejected", "#ff7b72")
        stats_layout.addWidget(self.sb_total)
        stats_layout.addWidget(self.sb_acc)
        stats_layout.addWidget(self.sb_rej)
        layout.addLayout(stats_layout)
        # buttons
        btn_layout = QHBoxLayout()
        self.btn_load = QPushButton("Load Files")
        self.btn_csv = QPushButton("Export CSV")
        self.btn_pdf = QPushButton("Export PDF")
        self.btn_fasta = QPushButton("Export cleaned FASTA")
        btn_layout.addWidget(self.btn_load)
        btn_layout.addWidget(self.btn_csv)
        btn_layout.addWidget(self.btn_pdf)
        btn_layout.addWidget(self.btn_fasta)
        layout.addLayout(btn_layout)
        # table
        self.table = QTableWidget()
        self.table.setColumnCount(9)
        self.table.setHorizontalHeaderLabels([
            "ID", "Length", "GC%", "N%", "Ambiguous%", "Homopolymer", "MeanPhred", "FastQC_Total", "Status"
        ])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        layout.addWidget(self.table)
        # signals
        self.btn_load.clicked.connect(self.load_files)
        self.btn_csv.clicked.connect(self.export_csv)
        self.btn_pdf.clicked.connect(self.export_pdf)
        self.btn_fasta.clicked.connect(self.export_fasta)

    def load_files(self):
        files, _ = QFileDialog.getOpenFileNames(
            self,
            "Select sequencing files and/or FastQC reports",
            "",
            "Seq Files (*.fasta *.fa *.fna *.fastq *.fq *.gb *.gbk *.embl);;FastQC (*.zip);;All Files (*)"
        )
        if not files:
            return
        self.accepted.clear()
        self.rejected.clear()
        self.fastqc_reports.clear()
        total_files = len(files)
        processed = 0
        for f in files:
            if f.lower().endswith("_fastqc.zip") or os.path.basename(f) == "fastqc_data.txt" or os.path.isdir(f):
                rep = parse_fastqc(f)
                if rep:
                    self.fastqc_reports[os.path.basename(f)] = rep
                processed += 1
                self.sb_total.update(processed, total_files, f"Files: {processed}/{total_files}")
                continue
            try:
                records, rej, stats = process_file(f, do_trim=True)
            except Exception as e:
                print(f"Error processing {f}: {e}")
                continue
            base = os.path.splitext(os.path.basename(f))[0]
            matched = None
            for k in self.fastqc_reports.keys():
                if k.startswith(base):
                    matched = self.fastqc_reports[k]
                    break
            if matched:
                for r in records:
                    r.fastqc_info = matched
                for r, _ in rej:
                    r.fastqc_info = matched
            self.accepted.extend(records)
            self.rejected.extend(rej)
            processed += 1
            self.sb_total.update(processed, total_files, f"Files: {processed}/{total_files}")
        self.update_stats()
        self.update_table()

    def update_stats(self):
        total = len(self.accepted) + len(self.rejected)
        acc = len(self.accepted)
        rej = len(self.rejected)
        self.sb_total.update(total, max(total,1), f"Total: {total}")
        self.sb_acc.update(acc, max(total,1), f"Accepted: {acc}")
        self.sb_rej.update(rej, max(total,1), f"Rejected: {rej}")
        self.summary = compute_summary(self.accepted, self.rejected, self.fastqc_reports)

    def update_table(self):
        rows = self.accepted + [r for r, _ in self.rejected]
        self.table.setRowCount(len(rows))
        for i, r in enumerate(rows):
            status = "ACCEPTED" if r in self.accepted else "REJECTED"
            fastqc_total = ""
            if r.fastqc_info:
                fastqc_total = str(r.fastqc_info.get("Total Sequences", ""))
            self.table.setItem(i, 0, QTableWidgetItem(r.id))
            self.table.setItem(i, 1, QTableWidgetItem(str(r.length)))
            self.table.setItem(i, 2, QTableWidgetItem(f"{r.gc_pct:.2f}"))
            self.table.setItem(i, 3, QTableWidgetItem(f"{r.n_pct:.4f}"))
            self.table.setItem(i, 4, QTableWidgetItem(f"{r.ambiguous_pct:.4f}"))
            self.table.setItem(i, 5, QTableWidgetItem(str(r.homopolymer_max)))
            self.table.setItem(i, 6, QTableWidgetItem(f"{r.mean_phred:.1f}" if r.mean_phred is not None else ""))
            self.table.setItem(i, 7, QTableWidgetItem(fastqc_total))
            self.table.setItem(i, 8, QTableWidgetItem(status))

    def export_csv(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save CSV", "", "CSV Files (*.csv)")
        if not path:
            return
        write_sequence_report(self.accepted, self.rejected, path)

    def export_pdf(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save PDF", "", "PDF Files (*.pdf)")
        if not path:
            return
        write_pdf_report(self.summary, self.accepted, self.rejected, path)

    def export_fasta(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save cleaned FASTA (accepted)", "", "FASTA Files (*.fasta *.fa)")
        if not path:
            return
        accepted_path = path if path.endswith(".fasta") or path.endswith(".fa") else path + ".fasta"
        rejected_path = os.path.splitext(accepted_path)[0] + "_rejected.fasta"
        write_fasta(self.accepted, accepted_path)
        write_fasta([r for r, _ in self.rejected], rejected_path)

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
