#!/bin/bash
#SBATCH -J EDTA
#SBATCH -o EDTA.o%j
#SBATCH -c 16
#SBATCH --mem=24G
#SBATCH -t 08:00:00 
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

conda activate EDTA-env

#cds file in genomic sequence format (.fna)
export cds="prot.fasta.codingseq"

perl /project/balan/anaconda3/envs/EDTA-env/bin/EDTA.pl --genome ${fasta} --cds ${cds} --sensitive 1 --anno 1 --evaluate 1 --species others -t 20

conda deactivate 
