# clean_fasta.py
import os

input_path = "data/raw/phrogs.fasta"
output_path = "data/raw/phrogs_cleaned.fasta"

print("Czyszczenie pliku FASTA...")

with open(input_path, 'r', encoding='utf-8', errors='ignore') as infile:
    with open(output_path, 'w', encoding='utf-8') as outfile:
        for line in infile:
            # Usuwamy wszystkie znaki sterujące (ASCII < 32), zostawiając tylko znak nowej linii (10)
            cleaned_line = "".join([c for c in line if ord(c) >= 32 or ord(c) == 10])
            outfile.write(cleaned_line)

print(f"Gotowe! Czysty plik zapisano jako: {output_path}")