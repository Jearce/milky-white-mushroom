from typing import List, Dict

import os
import glob
import subprocess
import argparse

from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser()
parser.add_argument("--proteins_dir")
parser.add_argument("--busco_results")
args = parser.parse_args()

proteins_dir = args.proteins_dir
busco_results = args.busco_results

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
MUSCLE_PATH = os.path.join(ROOT_DIR, "parallel_muscle.sh")


def make_dir(dir_name):
    if os.path.exists(dir_name):
        os.rmdir(dir_name)
    os.mkdir(dir_name)
    return dir_name


# after running muscle convert fasta alignments to phylip format and fix names for tree construction
def get_organism_name(des):
    # get index of open and close brackets looking from right to left
    first = len(des) - des[::-1].find("[")
    last = (len(des) - des[::-1].find("]")) - 1
    return des[first:last]


proteome_map = {
    os.path.basename(file): SeqIO.index(file, "fasta")
    for file in glob.glob(f"{proteins_dir}/*.faa")
}


def is_not_comment_and_is_complete(line):
    return (
        not line.startswith("#") and line.split("\t")[1].lower() == "complete"
    )


result_file = "run_basidiomycota_odb10/full_table.tsv"


busco_results_dir = [
    file
    for file in os.listdir(busco_results)
    if file.startswith("busco") and file.endswith(".faa")
]

# get sequences that match with the same busco id
grouped_by_busco_id: Dict[str, List[SeqRecord]] = {}
for root_dir in busco_results_dir:
    with open(f"{busco_results}/{root_dir}/{result_file}") as file:
        basename = root_dir.replace("busco_", "")
        for line in file:
            if is_not_comment_and_is_complete(line):
                busco_id, seq_id = line.strip().split("\t")[:3:2]
                if busco_id not in grouped_by_busco_id:
                    grouped_by_busco_id[busco_id] = []

                group: List[SeqRecord] = grouped_by_busco_id[busco_id]
                group.append(proteome_map[basename][seq_id])

results_dir = make_dir("seqs_grouped_by_busco_id")
for busco_id, protein_set in grouped_by_busco_id.items():
    if len(protein_set) == len(busco_results_dir):
        SeqIO.write(protein_set, f"{results_dir}/{busco_id}.faa", "fasta")

if not os.path.exists(MUSCLE_PATH):
    print("cannot find parellel_muscle.sh")
    exit()

alignments_dir= "alignments_by_busco_id"
subprocess.run(
    [
        "bash",
        MUSCLE_PATH,
        "seqs_grouped_by_busco_id",
        alignments_dir,
    ],
    check=True,
)

# fasta alignment to phylip format
phylip_dir = make_dir("phylip_alignments")
for file_name in glob.glob(f"{alignments_dir}/*.aln"):
    records = list(SeqIO.parse(file_name, "fasta"))
    for i, record in enumerate(records):
        des = record.description

        name = ""
        if record.id.startswith("g"):
            name = "Calocybe indica"
        else:
            name = get_organism_name(des)

        record.description = ""
        record.id = name.replace(" ", "_").replace(".", "")
        records[i] = record

    phylip_file = os.path.basename(file_name).replace("aln", "phy")
    SeqIO.write(records, f"{phylip_dir}/{phylip_file}", "phylip-relaxed")

# sort and concatenate all alignments
result_file = "combined.phy"
combined = None
for i, file_name in enumerate(glob.glob(f"{phylip_dir}/*.phy")):
    align = AlignIO.read(file_name, "phylip-relaxed")
    align.sort()
    if i == 0:
        combined = align
    else:
        combined += align

SeqIO.write(combined, result_file, "phylip-relaxed")
