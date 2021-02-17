#!/bin/bash
#SBATCH -J blastp
#SBATCH -o blastp.o%j
#SBATCH -c 20
#SBATCH --cpus-per-task=5
#SBATCH --mem=18G
#SBATCH -t 10:00:00

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jearce2@uh.edu
set -e
set -u
set -o pipefail

#set local anaconda3 env located in project dir
set +eu
source ~/.bashrc 

if [ $# -ne 2 ]
then
	echo "Incorrect number of arguments"
	echo "run_blastp.sh <query> <db_name> <out>"
	exit
fi

export query=$1
export db_name=$2
export out=$3

conda activate blast-env

blastp -query ${query} -db ${db_name} -outfmt 7 -out ${out} -num_threads 20 -evalue 1e-6

conda deactivate
