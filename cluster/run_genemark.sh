#!/bin/bash
#SBATCH -J gmes
#SBATCH -o gmes.o%j
#SBATCH -c 20
#SBATCH --mem=12G
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
  echo "run_genemark.sh <genome>"
  exit
fi

export genome=$1

conda activate genemark-env

gmes_petap.pl --sequence ${genome} --ES --fungus --cores 20

conda deactivate
