#!/bin/bash
#SBATCH -J repeat
#SBATCH -o repeat.o%j
#SBATCH -c 16
#SBATCH --mem=18G
#SBATCH -t 06:00:00 
#SBATCH --mail-type=END,FAIL

set -e
set -u 
set -o pipefail

set +eu
source ~/.bashrc 

#Prerequisites: Updated Dfam library added to ${ENV}/share/RepeatMasker/Libraries/ and the NINJA tool

if [ $# -ne 1 ]
then
  echo "Incorrect number of arguments"
  echo "${0} <contigs fasta file>"
  exit
fi

export fasta=$1

#find segmental duplications using default settings, followed by reversed-complemented
asgart -k 101 --out mushroom ${fasta}
asgart -k 101 -RCv --out mushroom_RC ${fasta}

#plotting duplications in chord format 
asgart-plot mushroom.json --colorize by-position --out=mushroom_chord.svg chord
asgart-plot mushroom_RC.json --colorize by-position --out=mushroom_RC_chord.svg chord

conda activate repeat-env

export NINJA_DIR="../../../../anaconda3/envs/repeat-env/bin/NINJA-0.97-cluster_only/NINJA"

BuildDatabase -name mushroom -engine ncbi ${fasta} 

RepeatModeler -database mushroom -engine ncbi -pa 15 -LTRStruct -ninja_dir ${NINJA_DIR}  

conda deactivate
