#!/bin/bash
#SBATCH -J wengan
#SBATCH -o wengan.o%j
#SBATCH -c 20
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH -t 08:00:00

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jearce2@uh.edu
set -e
set -u
set -o pipefail

#set local anaconda3 env located in project dir
set +eu
source ~/.bashrc

if [ $# -ne 3 ]
then
  echo "Incorrect number of arguments"
  echo "run_wengan.sh <nanoreads> <r1> <r2>"
  exit
fi

export nano=$1; export r1=$2; export r2=$3;

conda activate wengan-runtime-env

export LD_LIBRARY_PATH=/project/balan/anaconda3/envs/wengan-runtime-env/lib
prefix="white";
wengan.pl -x ontraw -a M -s ${r1},${r2} -l ${nano} -p ${prefix} -t 20 -g 28

conda deactivate

#set up polishing environment 
conda activate polish-env
module load SAMtools/1.9-intel-2017b

#three rounds of polishing
assembly="${prefix}.SPolished.asm.wengan.fasta"
for i in {1..3}
do
  polish_dir="polish_round_${i}"
  mkdir ${polish_dir}
  sorted_aln="${polish_dir}/aln-sorted.bam"
  pilon_out="${polish_dir}/pilon_out"
  minimap2 -ax sr ${assembly} ${r1} ${r2} | samtools view -u | samtools sort -@ ${t} > ${sorted_aln}
  samtools index ${sorted_aln}
  pilon --genome ${assembly} --frags ${sorted_aln} --outdir ${pilon_out}
  assembly="${pilon_out}/pilon.fasta"
done

conda deactivate
