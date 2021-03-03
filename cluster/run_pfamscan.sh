#!/bin/bash
#SBATCH -J pfam 
#SBATCH -o pfam.o%j
#SBATCH -c 16
#SBATCH --mem=10G
#SBATCH -t 02:00:00
#SBATCH --mail-type=END,FAIL
set -e
set -u
set -o pipefail

set +eu
source ~/.bashrc

if [ $# -ne 1 ]
then
   echo "Incorrect number of arguments"
   echo "${0} <fasta file>"
   exit
fi

export fasta=$1
export db="/project/balan/milky-white-mushoom/db/pfamdb"

conda activate pfam-env

pfam_scan.pl -fasta ${fasta} -dir ${db} -cpu 16 -outfile pfam_hits.txt -json pretty

conda deactivate 
