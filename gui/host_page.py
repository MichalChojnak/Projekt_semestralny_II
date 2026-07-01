from PyQt6.QtWidgets import QWidget, QLabel, QVBoxLayout


class HostPage(QWidget):

    def __init__(self):
        super().__init__()

        layout = QVBoxLayout(self)

        title = QLabel("Host Prediction")

        layout.addWidget(title)

        layout.addStretch()