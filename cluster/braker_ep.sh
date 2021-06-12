#!/bin/bash
#SBATCH -J br-ep
#SBATCH -o br-ep.o%j
#SBATCH -c 24
#SBATCH --mem=32G
#SBATCH -t 05:00:00
#SBATCH --mail-type=END,FAIL

set -e
set -u
set -o pipefail

set +eu
source ~/.bashrc

export genome=$1
export proteins=$2
export dir=$3
export t=24

conda activate braker

export ALIGNMENT_TOOL_PATH="/project/balan/bin/ProtHint-2.6.0/dependencies"
export PROTHINT_PATH="/project/balan/bin/ProtHint-2.6.0/bin"
export GENEMARK_PATH="/project/balan/bin/gmes_linux_64"

braker.pl --species="Calocybe_indi" --genome=${genome} --prot_seq=${proteins} --epmode \
 --softmasking --fungus --cores ${t} --workingdir=${dir}

conda deactivate
