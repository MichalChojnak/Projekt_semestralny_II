# protein_classifier.py
import re

# Punkt 13: Nowe wagi
PROTEIN_WEIGHTS = {
    "Receptor Binding Protein (RBP)": 20,
    "Tail Fiber": 18,
    "Tail Spike": 17,
    "Depolimerazy": 15,
    "Adsorpcja": 14,
    "Baseplate": 10,
    "Lizyny i Holiny": 4,  # Uśrednione wg Twojej tabeli (Endolysin 4, Holin 3) lub rozdzielone poniżej
    "Endolizyny": 4,
    "Holiny": 3,
    "Kapsyd": 1,
    "Białka portalu": 1,
    "Polimerazy": 1,
    "Terminazy": 1,
    "Integrazy": 1,
    "Pozostałe": 0.1
}
CATEGORY_COLORS = {
    "Receptor Binding Protein (RBP)": "#ff7f0e",
    "Tail Fiber": "#ffbb78",
    "Tail Spike": "#ff9896",
    "Depolimerazy": "#8c564b",
    "Adsorpcja": "#e377c2",
    "Baseplate": "#1f77b4",
    "Lizyny i Holiny": "#d62728",
    "Endolizyny": "#d62728",
    "Holiny": "#ff4d4d",
    "Kapsyd": "#2ca02c",
    "Białka portalu": "#17becf",
    "Polimerazy": "#9467bd",
    "Terminazy": "#c5b0d5",
    "Integrazy": "#f7b6d2",
    "Pozostałe": "#7f7f7f",
    "Brak adnotacji (ORFan)": "#c7c7c7"
}

# Punkt 1, 11 i 12: Jedna potężna funkcja klasyfikująca z wieloma synonimami
def kategoryzuj_bialko(stitle):
    if not isinstance(stitle, str):
        return "Pozostałe"

    opis = stitle.lower()

    # Rozszerzone synonimy
    if any(x in opis for x in ["rbp", "receptor binding", "receptor-binding", "host recognition"]):
        return "Receptor Binding Protein (RBP)"
    if any(x in opis for x in
           ["tail fiber", "tail fibre", "tail-fiber", "long tail fiber", "short tail fiber", "ltf", "distal tail"]):
        return "Tail Fiber"
    if any(x in opis for x in ["tail spike", "tail-spike", "tailspike"]):
        return "Tail Spike"
    if any(x in opis for x in ["depolymerase", "lyase", "hydrolase", "pectate"]):
        return "Depolimerazy"
    if any(x in opis for x in ["adsorption", "attachment"]):
        return "Adsorpcja"
    if any(x in opis for x in ["baseplate", "base plate", "wedge"]):
        return "Baseplate"
    if any(x in opis for x in ["endolysin", "lysin", "peptidoglycan", "murein", "muramidase"]):
        return "Endolizyny"
    if any(x in opis for x in ["holin", "pinholin"]):
        return "Holiny"
    if any(x in opis for x in ["capsid", "coat", "head", "major capsid", "minor capsid"]):
        return "Kapsyd"
    if any(x in opis for x in ["portal", "head-tail", "connector"]):
        return "Białka portalu"
    if any(x in opis for x in ["polymerase", "replisome", "primase"]):
        return "Polimerazy"
    if any(x in opis for x in ["terminase", "packaging"]):
        return "Terminazy"
    if any(x in opis for x in ["integrase", "recombinase", "excisionase"]):
        return "Integrazy"

    return "Pozostałe"