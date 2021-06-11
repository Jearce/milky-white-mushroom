#!/bin/bash
#SBATCH -J br-et
#SBATCH -o br-et.o%j
#SBATCH -c 24
#SBATCH --mem=32G
#SBATCH -t 05:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joshua.zhuang7@gmail.com

set -e
set -u
set -o pipefail

set +eu
source ~/.bashrc

export genome=$1
export bam=$2
export dir=$3
export t=24

conda activate braker

export ALIGNMENT_TOOL_PATH="/project/balan/bin/ProtHint-2.6.0/dependencies"
export PROTHINT_PATH="/project/balan/bin/ProtHint-2.6.0/bin"
export GENEMARK_PATH="/project/balan/bin/gmes_linux_64"

braker.pl --species="Calocye_indi" --genome=${genome} --bam=${bam} \
 --softmasking --fungus --cores ${t} --workingdir=${dir}

conda deactivate
