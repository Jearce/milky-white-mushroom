#!/bin/bash
#SBATCH -J blastp
#SBATCH -o blastp.o%j
#SBATCH -c 20
#SBATCH --mem=18G
#SBATCH -t 10:00:00
#SBATCH --mail-type=END,FAIL

set -e
set -u
set -o pipefail

#set local anaconda3 env located in project dir
set +eu
source ~/.bashrc 


if [ $# -ne 1 ]
then
  echo "Incorrect number of arguments"
  echo "${0} <protein fasta file>"
  exit
fi

export fasta=$1

conda activate java11-env

interproscan.sh \
  -d interpro_out\
  -b mwm \
  -f tsv, gff3, html, json\
  -goterms\
  -cpu 20\
  -i ${fasta}

conda deactivate 
