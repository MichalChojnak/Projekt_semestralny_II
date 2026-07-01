from PyQt6.QtWidgets import QWidget, QLabel, QVBoxLayout


class ORFPage(QWidget):

    def __init__(self):
        super().__init__()

        layout = QVBoxLayout(self)

        layout.addWidget(QLabel("ORF Finder"))

        layout.addStretch()