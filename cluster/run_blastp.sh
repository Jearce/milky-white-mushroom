#!/bin/bash
#SBATCH -J blastp
#SBATCH -o blastp.o%j
#SBATCH -c 20
#SBATCH --cpus-per-task=5
#SBATCH --mem=18G
#SBATCH -t 05:00:00
#SBATCH --mail-type=END,FAIL

set -e
set -u
set -o pipefail

#set local anaconda3 env located in project dir
set +eu
source ~/.bashrc 

if [ $# -ne 3 ]
then
	echo "Incorrect number of arguments"
	echo "${0} <query> <db_name> <out>"
	exit
fi

export query=$1; export db_name=$2; export out=$3;

conda activate blast-env

blastp -query ${query} -db ${db_name}\
	-outfmt "6 qacc sacc qstart qend sstart send evalue pident qlen slen bitscore qcovs qcovhsp"\
	-out ${out}\
	-num_threads 20\
	-evalue 1e-6

conda deactivate
