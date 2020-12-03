#!/usr/bin/env python
import re
import sys

import pandas
from Bio import Entrez

csv_path = sys.argv[1]

df = pandas.read_csv(csv_path).dropna(axis=1,how="all")
try:
    for name in df["Name"]:
        new_name = re.sub(r"[\s-]","_",re.sub("[.()]","",name)).lower()
        print(new_name)
except KeyError as e:
    print(f"Unexpected data format. No column with name: {e}")
