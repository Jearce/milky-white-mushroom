#!/bin/bash
set -e
set -u
set -o pipefail


if [ $# -eq 0 ]
then
  echo "usage: ./download_fungal_genomes.sh <csv file>"
  echo " Any failed downloads are logged to fungal.downloads.log file."
  exit
fi

export filename=$1
export ftp_dir="ftp://ftp.ensemblgenomes.org/pub/fungi/release-49/fasta/fungi_basidiomycota1_collection"

for name in $(extract_and_format_fungi_names.py ${filename} --strict)
do
  echo "------------------------------------"
  echo "Downloading proteins for ${name}"
  echo "-------------------------------------"
  if ! wget "${ftp_dir}/${name}/cds/*.cds.all.fa.gz"
  then
    echo "${name} was not downloaded" >> fungal.downloads.log
  fi
done
exit
