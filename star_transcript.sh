#!/bin/bash
#SBATCH -J star
#SBATCH -o star.o%j
#SBATCH -c 24
#SBATCH --mem=128G
#SBATCH -t 05:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joshua.zhuang7@gmail.com

set -e
set -u
set -o pipefail

set +eu
source ~/.bashrc

t=24
export index_dir=$1
export fasta=$2
export rna1=$3
export rna2=$4
export rna3=$5
export rna4=$6

conda activate star

STAR --runThreadN ${t} --runMode genomeGenerate --genomeDir ${index_dir} --genomeFastaFiles ${fasta} \
 --genomeSAindexNbases 11

STAR --runThreadN ${t} --genomeDir ${index_dir} --readFilesIn ${rna1},${rna3} ${rna2},${rna4} \
 --quantTranscriptomeBan Singleend --outSAMtype BAM Unsorted SortedByCoordinate --outFileNamePrefix mwm --twopassMode Basic \
 --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --outFilterMismatchNmax 999 \
 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 20000 --alignMatesGapMax 20000 \
 --outFilterIntronMotifs RemoveNoncanonical --outReadsUnmapped Fastx

conda deactivate 
