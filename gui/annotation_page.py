from PyQt6.QtWidgets import QWidget, QLabel, QVBoxLayout


class AnnotationPage(QWidget):

    def __init__(self):
        super().__init__()

        layout = QVBoxLayout(self)

        layout.addWidget(QLabel("Protein Annotation"))

        layout.addStretch()