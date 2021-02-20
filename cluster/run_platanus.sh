#!/bin/bash
#SBATCH -J mushroom
#SBATCH -o mushroom.o%j
#SBATCH -c 16
#SBATCH --mem=18G
#SBATCH -t 03:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<your_email>
set -e
set -u
set -o pipefail

if [ $# -ne 2 ]
then
   echo "Incorrect number of arguments"
   echo "run_platanus.sh <r1> <r2>"
   exit
fi

export r1=$1; export r2=$2;

platanus assemble -o Mwm -t $SLURM_JOB_CPUS_PER_NODE -m 16 \
-f ${r1}  ${r2} \
>mwm_log.txt

