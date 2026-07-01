from PyQt6.QtWidgets import (
    QWidget,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QFrame
)


class HomePage(QWidget):

    def __init__(self):
        super().__init__()

        layout = QVBoxLayout(self)

        title = QLabel("Welcome to Phage Host Predictor")

        title.setStyleSheet("""
            QLabel{
                font-size:24px;
                font-weight:bold;
            }
        """)

        subtitle = QLabel(
            "Pipeline:\n"
            "QC → ORF Detection → Annotation → Feature Selection → Host Prediction"
        )

        subtitle.setStyleSheet("""
            QLabel{
                font-size:14px;
            }
        """)

        btn_project = QPushButton("Open Genome")

        btn_project.setMinimumHeight(45)

        pipeline = QFrame()

        pipeline.setFrameShape(QFrame.Shape.Box)

        pipeline_layout = QVBoxLayout(pipeline)

        pipeline_layout.addWidget(
            QLabel(
                "Pipeline status\n\n"
                "① QC\n"
                "② ORF Finder\n"
                "③ Annotation\n"
                "④ Feature Selection\n"
                "⑤ Host Prediction"
            )
        )

        layout.addWidget(title)
        layout.addWidget(subtitle)
        layout.addSpacing(20)
        layout.addWidget(btn_project)
        layout.addSpacing(20)
        layout.addWidget(pipeline)
        layout.addStretch()