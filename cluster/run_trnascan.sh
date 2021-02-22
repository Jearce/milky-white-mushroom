#!/bin/bash
#SBATCH -J trna
#SBATCH -o trna.o%j
#SBATCH -c 12
#SBATCH --mem=10G
#SBATCH -t 02:00:00
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
	echo "$0 <contigs>"
	exit
fi

export contigs=$1

conda activate trnascan-se-env

tRNAscan-SE -H -Q -E -o trnas.txt -f trnas_stuctures.txt -m trna.models -a trans.fasta -s isospecific.txt ${contigs}

conda deactivate
