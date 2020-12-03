import re
import sys

import pandas
from Bio import Entrez

csv_path = sys.argv[1]
strict = None
if len(sys.argv) == 3: 
    strict = sys.argv[2]

df = pandas.read_csv(csv_path).dropna(axis=1,how="all")
try:
    for name in df["Name"]:
        new_name = re.sub(r"[\s-]","_",re.sub("[.()]","",name)).lower()
        if strict == "--strict":
            print(new_name)
        else:
            genus_species = "_".join(new_name.split("_")[:2])
            print(genus_species)
except KeyError as e:
    print(f"Unexpected data format. No column with name: {e}")
