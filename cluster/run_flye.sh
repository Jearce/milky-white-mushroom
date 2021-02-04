#!/bin/bash
#SBATCH -J flye
#SBATCH -o flye.o%j
#SBATCH -c 20
#SBATCH --mem=18G
#SBATCH -t 15:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jearce2@uh.edu
set -e
set -u
set -o pipefail

set +eu
source ~/.bashrc

if [ $# -ne 4 ]
then
  echo "Incorrect number of arguments"
  echo "run_flye.sh <nanoreads> <r1> <2> <outdir>"
  exit
fi

conda activate flye-env

export nanoreads=$1
export r1=$2
export r2=$3
export prefix=$4

export assembly_folder=${prefix}_flye
flye --nano-raw ${nanoreads} --out-dir ${assembly_folder} --threads 20 -i 3

cd ${assembly_folder}

#set up polishing environment 
conda activate polish-env
module load SAMtools/1.9-intel-2017b

#three rounds of polishing
assembly="assembly.fasta"
for i in {1..3}
  polish_dir="polish_round_${i}"
  mkdir polish_folder

  aln="${polish_dir}/aln-sorted.bam"
  pilon_out="${polish_dir}/pilon_out"
  minimap3 -ax -t 18 sr ${assembly} ${r1} ${r2} | samtools view -u | samtools sort -@ 20 > ${aln}
  pilon --genome ${assembly} --frags ${aln} --outdir ${pilon_out}
  assembly="${pilon_out}/pilon.fasta"

do


conda deactivate
