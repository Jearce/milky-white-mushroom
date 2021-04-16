#!/bin/bash
#SBATCH -J kmergenie
#SBATCH -o kmergenie.o%j
#SBATCH -c 16
#SBATCH --mem=18G
#SBATCH -t 06:00:00
#SBATCH --mail-type=END,FAIL

set -e
set -u
set -o pipefail

set +eu
source ~/.bashrc


if [ $# -ne 2 ]
then
  echo "Incorrect number of arguments"
  echo "${0} <read #1 fasta file> <read #2 fasta file>"
  exit
fi

export fasta=$1
export fasta2=$2

conda activate repeat-env

kmergenie ${fasta} -o read1_histogram
kmergenie ${fasta2} -o read2_histogram

conda deactivate


