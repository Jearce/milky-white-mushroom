#!/bin/bash
#SBATCH -J brakerP
#SBATCH -o brakerP.o%j
#SBATCH -c 22
#SBATCH --mem=20G
#SBATCH -t 05:00:00
#SBATCH --mail-type=END,FAIL

set -e
set -u
set -o pipefail

#set local anaconda3 env located in project dir
set +eu
source ~/.bashrc 

if [ $# -ne 2 ]
then
	echo "Incorrect number of arguments"
	echo "run_braker.sh <softmasked_contigs> <proteins>"
	exit
fi

export contigs=$1; export proteins=$2;export t=22;

if ! egrep -e "[atgc]" ${contigs}
then
	echo "${contigs} must be softmasked"
	exit
fi
conda activate braker-env

prothint.py ${contigs} ${proteins} --fungus --threads ${t}

export hints="prothint_augustus.gff"

if [ -f ${hints} ]
then
	braker.pl --species="calocybe_indica" --genome=${contigs} --hints=${hints} --epmode --softmasking --fungus --cores ${t}
else
	"${hints} does not exist"
fi

conda deactivate
