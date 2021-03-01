#!/bin/bash
#SBATCH -J dbcan
#SBATCH -o dbcan.o%j
#SBATCH -c 12
#SBATCH --mem=12G
#SBATCH -t 03:00:00
#SBATCH --mail-type=END,FAIL
set -e
set -u
set -o pipefail

export t=12

if [ $# -eq 0 ]
then
  echo "Incorrect number of arguments"
  echo "${0} [fasta files]"
  exit
fi

#set local anaconda3 env located in project dir
set +eu
source ~/.bashrc 

conda activate dbcan-env

export db_dir="/project/balan/milky-white-mushoom/db/dbcan"
for file in ${@}
do
  run_dbcan.py ${file} protein \
    --out_dir $(basename ${file})\
    --tools all\
    --hmm_cpu ${t}\
    --dia_cpu ${t}\
    --hotpep_cpu ${t}\
    --db_dir ${db_dir}
done

conda deactivate
