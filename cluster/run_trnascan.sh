#!/bin/bash
#SBATCH -J trna
#SBATCH -o trna.o%j
#SBATCH -c 12
#SBATCH --mem=10G
#SBATCH -t 02:00:00
#SBATCH --mail-type=END,FAIL

set -e
set -u
set -o pipefail

#set local anaconda3 env located in project dir
set +eu
source ~/.bashrc 


if [ $# -eq 0 ]
then
	echo "Incorrect number of arguments"
	echo "$0 [contigs files]"
	exit
fi

conda activate trnascan-se-env

for contig in ${@}
do
  export out_dir=$(basename -- "${file%.*}")
  mkdir ${out_dir}
  tRNAscan-SE -H -Q -E\
    -o ${out_dir}/trnas.txt\
    -f ${out_dir}/trnas_stuctures.txt\
    -m ${out_dir}/trna.models\
    -a ${out_dir}/trans.fasta\
    -s ${out_dir}/isospecific.txt\
    ${contigs}
done

conda deactivate
