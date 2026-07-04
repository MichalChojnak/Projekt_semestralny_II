import csv
import time

# --- NAZWY PLIKÓW ---
# Upewnij się, że nazwy poniżej dokładnie odpowiadają plikom w Twoim folderze!
plik_phrogs = "phrogs_table_almostfinal_plusGO_wNA_utf8.tsv"
plik_ncbi = "sequences.csv"
plik_szablon = "taxonomy.csv"  # Ten plik wygenerowaliśmy wcześniej (z 'unknown')
plik_gotowy = "taxonomy_final.csv"  # To będzie Twój ostateczny plik dla DIAMOND

slownik_taksonomii = {}
start_time = time.time()

# 1. Wczytywanie danych z bazy PHROGs
print("1/3 Wczytuję dane z bazy PHROGs...")
with open(plik_phrogs, mode='r', encoding='utf-8') as f:
    # Plik PHROGs z reguły jest rozdzielany tabulatorem
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        # Przechwytujemy nazwy kolumn (różne wersje pliku mogą mieć 'phrog' lub '#phrog')
        phrog_id = row.get('phrog') or row.get('#phrog')

        if phrog_id:
            # Usuwamy przecinki z nazw, żeby nie zepsuć formatu CSV
            species = row.get('Annotation', 'unknown').replace(',', '')
            genus = row.get('Category', 'unknown').replace(',', '')
            slownik_taksonomii[phrog_id] = (species, genus)

# 2. Wczytywanie danych z bazy NCBI
print("2/3 Wczytuję dane z bazy NCBI...")
with open(plik_ncbi, mode='r', encoding='utf-8') as f:
    # Plik z NCBI jest standardowym CSV
    reader = csv.DictReader(f, delimiter=',')
    for row in reader:
        acc = row.get('Accession')

        if acc:
            species = row.get('Species', 'unknown').replace(',', '')
            genus = row.get('Genus', 'unknown').replace(',', '')
            slownik_taksonomii[acc] = (species, genus)

# 3. Tworzenie ostatecznego pliku CSV
print("3/3 Podmieniam 'unknown' w głównym pliku. To potrwa kilka sekund...")
znaleziono = 0
nie_znaleziono = 0

with open(plik_szablon, mode='r', encoding='utf-8') as f_in, \
        open(plik_gotowy, mode='w', encoding='utf-8', newline='') as f_out:
    for line in f_in:
        # Przepisz nagłówek
        if line.startswith("sseqid"):
            f_out.write(line)
            continue

        parts = line.strip().split(',')
        seq_id = parts[0]

        # Sprawdzamy, czy nasze ID jest w słowniku z PHROGs lub NCBI
        if seq_id in slownik_taksonomii:
            species, genus = slownik_taksonomii[seq_id]
            # Zabezpieczenie przed pustymi komórkami
            if not species: species = "unknown"
            if not genus: genus = "unknown"

            f_out.write(f"{seq_id},{species},{genus}\n")
            znaleziono += 1
        else:
            # Jeśli jakimś cudem ID nie ma w obu plikach, zostawiamy 'unknown'
            f_out.write(line)
            nie_znaleziono += 1

czas_trwania = round(time.time() - start_time, 2)
print("-" * 30)
print(f"ZAKOŃCZONO w {czas_trwania} s!")
print(f"Pomyślnie dopasowano nazw: {znaleziono}")
print(f"Nie znaleziono w metadanych: {nie_znaleziono}")
print(f"Twój plik docelowy to: {plik_gotowy}")