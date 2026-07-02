from pathlib import Path
import subprocess
import shutil


class ORFManager:

    def __init__(self):

        self.phanotate = Path("external/phanotate.exe")

    def check_dependencies(self):

        if not self.phanotate.exists():

            raise FileNotFoundError(
                "Cannot find external/phanotate.exe"
            )

    def prepare_output(self, genome_path):

        genome_path = Path(genome_path)

        out_dir = Path("output") / genome_path.stem

        out_dir.mkdir(parents=True, exist_ok=True)

        return out_dir

    def run(self, genome_path):

        self.check_dependencies()

        out_dir = self.prepare_output(genome_path)

        gff_file = out_dir / "genes.gff"

        cmd = [

            str(self.phanotate),

            str(genome_path)

        ]

        with open(gff_file, "w") as f:

            subprocess.run(

                cmd,

                stdout=f,

                stderr=subprocess.PIPE,

                text=True,

                check=True

            )

        return gff_file