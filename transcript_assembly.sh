#!/bin/bash
#SBATCH -J rna-assemble
#SBATCH -o rna-assemble.o%j
#SBATCH -c 24
#SBATCH --mem=64G
#SBATCH -t 15:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joshua.zhuang7@gmail.com
set -e
set -u
set -o pipefail

set +eu
source ~/.bashrc

export r1=$1;export r2=$2

t=24

conda activate trinity
export trinity_out="../../../trinity/trinity_corrected"

Trinity --seqType fq --max_memory 64G --left ${1} --right ${2} --CPU 24 --output ${trinity_out}

conda deactivate
