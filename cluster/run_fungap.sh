#!/bin/bash
#SBATCH -J fungap
#SBATCH -o fungap.o%j
#SBATCH -c 30
#SBATCH --mem=32G
#SBATCH -t 6-10:00:00

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<your email>
set -e
set -u
set -o pipefail

#set local anaconda3 env located in project dir
set +eu
source ~/.bashrc

#set paths to prepare fungap correctly 
export PROJECT_DIR="/project/balan"
export FUNGAP_DIR="${PROJECT_DIR}/FunGAP"
export SISTER_DIR="${PROJECT_DIR}/milky-white-mushroom/annotations/sister_orgs"

if [ $# -ne 3 ]
then
  echo "  Incorrect number of arguments  "
  echo "run_fungap.sh <a> <r1> <r2>"
  echo "Arguments:"
  echo "  a: genome assembly in fasta format"
  echo "  r1: forward read of rna-seq data in fastq format"
  echo "  r2: reverse read of rna-seq data in fastq format"
  exit
fi


#get genome assembly and  rna reads
export assembly=$1; export r1=$2; export r2=$3;

#prepare fungap
conda activate maker
export MAKER_DIR=$(dirname $(which maker))
conda activate fungap
cd $FUNGAP_DIR
./set_dependencies.py --pfam_db_path db/pfam/ --genemark_path external/gmes_linux_64/ --maker_path ${MAKER_DIR}


cd /project/balan/milky-white-mushroom/annotations

#annotate with fungap
${FUNGAP_DIR}/fungap.py --output_dir fungap_out\
 --trans_read_1 ${r1} --trans_read_2 ${r2}\
 --genome_assembly ${ASSEMBLY}\
 --augustus_species coprinus_cinereus\
 --sister_proteome ${SISTER_DIR}/prot_db.faa\
 --num_cores 20 --busco_dataset basidoiomycota_odb10

#get out of fungap env
conda deactivate
#get out of maker env
conda deactivate
#get out of base env
conda deactivate
