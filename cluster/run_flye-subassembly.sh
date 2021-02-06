#!/bin/bash
#SBATCH -J flyeSub
#SBATCH -o flyeSub.o%j
#SBATCH -c 20
#SBATCH --mem=32G
#SBATCH -t 15:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jearce2@uh.edu
set -e
set -u
set -o pipefail

if [ $# -ne 3 ]
then
  echo "Incorrect number of arguments"
  echo "run_flye-subassembly.sh <assembly 1> <assembly 2> <r1> <r2> <outdir>"
  exit
fi

export a1=$1; export a2=$2; export r1=$3; export r2=$4; export out_dir=$5;

t=20

#set local anaconda3 env located in project dir
set +eu
source ~/.bashrc 

conda activate flye-env

flye --subassemblies ${a1} ${a2} --out-dir ${out_dir} -t ${t} -g 28m -i 3

conda deactivate 

#set up polishing environment 
conda activate polish-env
module load SAMtools/1.9-intel-2017b

#three rounds of polishing
assembly="${out_dir}/assembly.fasta"
for i in {1..3}
do
  polish_dir="${out_dir}/polish_round_${i}"
  mkdir ${polish_dir}
  sorted_aln="${polish_dir}/aln-sorted.bam"
  pilon_out="${polish_dir}/pilon_out"
  minimap2 -ax sr ${assembly} ${r1} ${r2} | samtools view -u | samtools sort -@ ${t} > ${sorted_aln}
  samtools index ${sorted_aln}
  pilon --genome ${assembly} --frags ${sorted_aln} --outdir ${pilon_out}
  assembly="${pilon_out}/pilon.fasta"
done

conda deactivate
