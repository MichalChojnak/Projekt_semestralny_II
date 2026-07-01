from PyQt6.QtWidgets import (
    QMainWindow,
    QWidget,
    QVBoxLayout,
    QTabWidget,
    QLabel,
    QStatusBar
)

from gui.home_page import HomePage
from gui.qc_page import QCPage
from gui.orf_page import ORFPage
from gui.annotation_page import AnnotationPage
from gui.host_page import HostPage
from gui.reports_page import ReportsPage


class MainWindow(QMainWindow):

    def __init__(self):
        super().__init__()

        self.setWindowTitle("Phage Host Predictor")
        self.resize(1400, 900)

        self.tabs = QTabWidget()

        self.tabs.addTab(HomePage(), "🏠 Home")
        self.tabs.addTab(QCPage(), "🧬 QC")
        self.tabs.addTab(ORFPage(), "🧬 ORF Finder")
        self.tabs.addTab(AnnotationPage(), "📖 Annotation")
        self.tabs.addTab(HostPage(), "🎯 Host Prediction")
        self.tabs.addTab(ReportsPage(), "📄 Reports")

        central = QWidget()

        layout = QVBoxLayout(central)

        title = QLabel("Phage Host Predictor")
        title.setStyleSheet("""
            QLabel{
                font-size:28px;
                font-weight:bold;
            }
        """)

        layout.addWidget(title)
        layout.addWidget(self.tabs)

        self.setCentralWidget(central)

        self.status = QStatusBar()
        self.status.showMessage("Ready")

        self.setStatusBar(self.status)