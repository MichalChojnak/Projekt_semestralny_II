import os
import csv
import time
from Bio import Entrez, SeqIO
import urllib.error
import http.client

# ==========================================
# KONFIGURACJA NCBI
# ==========================================
# MUSISZ podać swój adres email, aby NCBI nie zablokowało Twojego IP!
Entrez.email = "chojnacki.michal84@gmail.com"

# Liczba genomów fagowych do pobrania (na początek 500 to dobry test)
# Pełna baza może wymagać ustawienia tej wartości na np. 10000
MAX_RECORDS = 10000
BATCH_SIZE = 50

OUTPUT_FILE = "host_database.csv"

# Słowa kluczowe białek, które nas interesują (żeby nie pobierać białek hipotetycznych)
INFORMATIVE_PROTEINS = [
    "tail fiber", "tail spike", "receptor binding", "rbp", "depolymerase",
    "adsorption", "baseplate", "lysin", "endolysin", "holin", "portal",
    "capsid", "coat", "head", "polymerase", "integrase", "terminase"
]


def is_informative(product_name):
    name_lower = product_name.lower()
    return any(keyword in name_lower for keyword in INFORMATIVE_PROTEINS)


def extract_host(record):
    """Próbuje wyciągnąć nazwę gospodarza z adnotacji rekordu GenBank."""
    # 1. Próba znalezienia jawnego pola 'host'
    for feat in record.features:
        if feat.type == "source":
            if "host" in feat.qualifiers:
                return feat.qualifiers["host"][0]

    # 2. Próba wywnioskowania z nazwy organizmu (np. "Escherichia phage T4" -> "Escherichia")
    organism = record.annotations.get("organism", "")
    if "phage" in organism.lower() or "virus" in organism.lower():
        parts = organism.split()
        if len(parts) > 0 and parts[0].lower() not in ["bacteriophage", "phage", "unidentified"]:
            return parts[0]

    return "Unknown"


def main():
    print("🔍 Krok 1: Wyszukiwanie kompletnych genomów fagów w NCBI...")
    try:
        # Szukamy fagów (vira), które mają kompletny genom
        handle = Entrez.esearch(db="nucleotide",
                                term="(phage[Title] OR bacteriophage[Title]) AND complete genome[Title]",
                                retmax=MAX_RECORDS)
        record = Entrez.read(handle)
        handle.close()

        id_list = record["IdList"]
        total_ids = len(id_list)
        print(f"✅ Znaleziono {total_ids} genomów do pobrania.")

    except Exception as e:
        print(f"❌ Błąd wyszukiwania: {e}")
        return

    print(f"💾 Krok 2: Pobieranie adnotacji i budowa {OUTPUT_FILE}...")

    with open(OUTPUT_FILE, mode='w', newline='', encoding='utf-8') as csv_file:
        writer = csv.writer(csv_file)
        # Nagłówki naszej tabeli
        writer.writerow(["accession", "protein_function", "host_genus", "host_species"])

        proteins_saved = 0

        for i in range(0, total_ids, BATCH_SIZE):
            batch_ids = id_list[i:i + BATCH_SIZE]
            print(f"   Pobieranie paczki {i + 1} do {min(i + BATCH_SIZE, total_ids)}...")

            try:
                fetch_handle = Entrez.efetch(db="nucleotide", id=batch_ids, rettype="gb", retmode="text")
                records = list(SeqIO.parse(fetch_handle, "genbank"))
                fetch_handle.close()

                for rec in records:
                    host = extract_host(rec)
                    if host == "Unknown":
                        continue

                    # Ekstrakcja rodzaju z nazwy gospodarza (np. "Escherichia coli" -> "Escherichia")
                    host_genus = host.split()[0]

                    # Szukanie białek (Features typu CDS)
                    for feature in rec.features:
                        if feature.type == "CDS":
                            product = feature.qualifiers.get("product", [""])[0]
                            protein_id = feature.qualifiers.get("protein_id", [""])[0]

                            # Zapisujemy tylko jeśli białko jest informatywne i ma numer akcesyjny
                            if protein_id and product and is_informative(product):
                                writer.writerow([protein_id, product, host_genus, host])
                                proteins_saved += 1

            except (urllib.error.HTTPError, http.client.IncompleteRead) as e:
                print(f"   ⚠️ Błąd sieci podczas pobierania paczki. Przeskakuję. ({e})")
            except Exception as e:
                print(f"   ⚠️ Nieoczekiwany błąd: {e}")

            # NCBI wymaga, by nie wysyłać więcej niż 3 zapytań na sekundę
            time.sleep(1)

    print("==================================================")
    print(f"🎉 Zakończono! Zapisano {proteins_saved} białek do pliku {OUTPUT_FILE}.")
    print("==================================================")


if __name__ == "__main__":
    main()