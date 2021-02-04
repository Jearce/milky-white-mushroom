#!/bin/bash
#SBATCH -J minimap2
#SBATCH -o minimap2.o%j
#SBATCH -c 18
#SBATCH --mem=18G
#SBATCH -t 05:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jearce2@uh.edu
set -e
set -u
set -o pipefail

#set local anaconda3 env located in project dir
set +eu
source ~/.bashrc

#activate genemark environment - genemark perl dependencies are avaiable here
conda activate minimap2-env
module load SAMtools/1.9-intel-2017b


export assembly=$2
export r1=$2
export r2=$3

minimap3 -ax -t 18 sr ${assembly} ${r1} ${r2} > aln.sam
samtools view -S -b aln.sam > aln.bam
samtools sort -@ 18 aln.sam > aln-sorted.bam
samtools index aln-sorted.bam

conda deactivate
