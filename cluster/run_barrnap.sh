#!/bin/bash
#SBATCH -J barrnap
#SBATCH -o barrnap.o%j
#SBATCH -c 10
#SBATCH --mem=18G
#SBATCH -t 00:30:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jearce2@uh.edu
set -e
set -u
set -o pipefail

if [ $# -ne 1 ]
then
  echo "Incorrect number of arguments"
  echo "run_barrnap.sh <contigs>"
  exit
fi


#set local anaconda3 env located in project dir
set +eu
source ~/.bashrc 

conda activate barrnap

export contigs=$1
barrnap --threads 10 --kingdom euk --evalue 1e-16 -o rrna.fa < ${contigs} > rrna.gff

conda deactivate


