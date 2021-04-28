#!/bin/bash
set -u
set -e
set -o pipefail

if [ $# -ne 2 ]
then
 echo "Incorrect number of args"
 echo "${0} <in dir> <out dir>"
 exit
fi

export in_dir=$1; export out_dir=$2;

if [ -d ${out_dir} ]
then
 echo "${out_dir} already exists"
 exit
fi

mkdir ${out_dir}

for file in $(ls ${in_dir}/*.{faa,fasta});
do
 basename=$(basename ${file})
 muscle -in $file -out ${out_dir}/${basename%.*}.aln &
done
wait
