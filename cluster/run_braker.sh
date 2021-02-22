#!/bin/bash
#SBATCH -J braker
#SBATCH -o braker.o%j
#SBATCH -c 20
#SBATCH --mem=15G
#SBATCH -t 05:00:00
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
	echo "run_braker.sh <softmasked_contigs>"
	exit
fi

export contigs=$1

conda activate braker-env

#pipeline: contigs --> genemark es --> augutus training --> augustus --> resutls
braker.pl --species="calocybe_indica" --genome=${contigs} --softmasking --esmode --fungus --cores 20

conda deactivate
