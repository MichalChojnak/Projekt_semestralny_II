import time

plik_wejsciowy = "phage_final_db_clean.fasta"
plik_wyjsciowy = "taxonomy.csv"

print("Rozpoczynam generowanie pliku taxonomy.csv...")
start_time = time.time()

with open(plik_wejsciowy, "r", encoding="utf-8") as fasta, open(plik_wyjsciowy, "w", encoding="utf-8") as csv:
    # Zapisz nagłówki kolumn wymagane przez DIAMOND
    csv.write("sseqid,species,genus\n")

    # Przeszukaj plik FASTA linia po linii
    for line in fasta:
        if line.startswith(">"):
            # Weź pierwszy wyraz, usuń znak '>'
            seq_id = line.split()[0][1:]

            # Zapisz ID z tymczasowym 'unknown'
            csv.write(f"{seq_id},unknown,unknown\n")

print(f"Gotowe! Plik CSV utworzony w {round(time.time() - start_time, 2)} sekund.")