#!/bin/bash
#SBATCH -J flyeSub
#SBATCH -o flyeSub.o%j
#SBATCH -c 20
#SBATCH --mem=18G
#SBATCH -t 05:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jearce2@uh.edu
set -e
set -u
set -o pipefail

if [ $# -ne 3 ]
then
  echo "Incorrect number of arguments"
  echo "run_flye-subassembly.sh <assembly 1> <assembly 2> <outdir>"
  exit
fi

export a1=$1; export a2=$2; export outdir=$3

#set local anaconda3 env located in project dir
set +eu
source ~/.bashrc 

conda activate flye-env

flye --subassemblies ${a1} ${a2} --out-dir ${outdir} -t 8 -g 28m -i 3

conda deactivate 
