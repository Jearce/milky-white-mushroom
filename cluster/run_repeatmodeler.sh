#!/bin/bash
#SBATCH -J repeat
#SBATCH -o repeat.o%j
#SBATCH -c 16
#SBATCH --mem=18G
#SBATCH -t 06:00:00 
#SBATCH --mail-type=END,FAIL

set -e
set -u 
set -o pipefail

set +eu
source ~/.bashrc 


if [ $# -ne 1 ]
then
  echo "Incorrect number of arguments"
  echo "${0} <contigs fasta file>"
  exit
fi

export fasta=$1

conda activate repeatmodeler-env

BuildDatabase -name mushroom ${fasta} 
RepeatModeler -database mushroom -pa 15

conda deactivate
