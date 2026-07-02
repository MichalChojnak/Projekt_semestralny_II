from pathlib import Path

from PyQt6.QtWidgets import (
    QWidget,
    QLabel,
    QPushButton,
    QFileDialog,
    QVBoxLayout,
    QHBoxLayout,
    QTableWidget,
    QTableWidgetItem,
    QTextEdit,
    QMessageBox
)

from PyQt6.QtCore import Qt


class ORFPage(QWidget):

    def __init__(self):
        super().__init__()

        self.input_file = ""

        layout = QVBoxLayout(self)

        #######################################################
        # TITLE
        #######################################################

        title = QLabel("ORF Finder")

        title.setStyleSheet("""
            font-size:24px;
            font-weight:bold;
        """)

        layout.addWidget(title)

        #######################################################
        # FILE
        #######################################################

        row = QHBoxLayout()

        self.file_label = QLabel("No genome selected")

        self.btn_upload = QPushButton("Upload Genome")

        row.addWidget(self.file_label)

        row.addStretch()

        row.addWidget(self.btn_upload)

        layout.addLayout(row)

        #######################################################
        # BUTTON
        #######################################################

        self.btn_find = QPushButton("Find ORFs")

        self.btn_find.setEnabled(False)

        layout.addWidget(self.btn_find)

        #######################################################
        # TABLE
        #######################################################

        self.table = QTableWidget()

        self.table.setColumnCount(5)

        self.table.setHorizontalHeaderLabels([
            "Gene",
            "Start",
            "End",
            "Length",
            "Strand"
        ])

        layout.addWidget(self.table)

        #######################################################
        # LOG
        #######################################################

        self.log = QTextEdit()

        self.log.setReadOnly(True)

        self.log.setMaximumHeight(150)

        layout.addWidget(self.log)

        #######################################################
        # SIGNALS
        #######################################################

        self.btn_upload.clicked.connect(self.select_genome)

        self.btn_find.clicked.connect(self.run_orf)

    ###########################################################

    def log_message(self, text):

        self.log.append(text)

    ###########################################################

    def select_genome(self):

        filename, _ = QFileDialog.getOpenFileName(

            self,

            "Select FASTA",

            "",

            "FASTA (*.fasta *.fa *.fna)"
        )

        if filename:

            self.input_file = filename

            self.file_label.setText(Path(filename).name)

            self.btn_find.setEnabled(True)

            self.log_message(f"Loaded genome:\n{filename}")

    ###########################################################

    def run_orf(self):

        if self.input_file == "":

            QMessageBox.warning(

                self,

                "Warning",

                "Select genome first."

            )

            return

        self.log_message("Starting ORF prediction...")

        #
        # tutaj później wywołamy PHANOTATE
        #

        self.log_message("ORF prediction finished.")

        #
        # TESTOWE DANE
        #

        fake_data = [

            ("gene_1", 120, 980, 286, "+"),

            ("gene_2", 1120, 1800, 226, "-"),

            ("gene_3", 2100, 3500, 466, "+")

        ]

        self.table.setRowCount(len(fake_data))

        for row, gene in enumerate(fake_data):

            for col, value in enumerate(gene):

                self.table.setItem(

                    row,

                    col,

                    QTableWidgetItem(str(value))

                )

        self.log_message("Displayed results.")