# Otwieramy plik w trybie binarnym ('rb'), żeby dokładnie odczytać wszystkie znaki
with open("phage_final_db.fasta", "rb") as f_in:
    zawartosc = f_in.read()

# Znak ASCII 26 w zapisie bajtowym to \x1a. Zastępujemy go pustością.
zawartosc_oczyszczona = zawartosc.replace(b'\x1a', b'')

# Zapisujemy nowy, czysty plik
with open("phage_final_db_clean.fasta", "wb") as f_out:
    f_out.write(zawartosc_oczyszczona)

print("Plik oczyszczony. Możesz użyć 'phage_final_db_clean.fasta'.")