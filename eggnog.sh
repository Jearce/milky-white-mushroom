#!/bin/bash
#SBATCH -J eggn
#SBATCH -o eggn.o%j
#SBATCH -c 20
#SBATCH --mem=64G
#SBATCH -t 03:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joshua.zhuang7@gmail.com

set -e
set -u
set -o pipefail

set +eu
source ~/.bashrc

conda activate

t=20
export LD_LIBRARY_PATH="/project/balan/mzhuang2/miniconda3/lib"
export data_dir="/project/balan/bin/eggnog-mapper-2.1.3/data"
export mmseqs_dir="/project/balan/bin/eggnog-mapper-2.1.3/data/fungi.mmseqs/fungi.mmseqs"
export proteins=$1

emapper.py -i ${proteins} -m mmseqs --mmseqs_db ${mmseqs_dir} -o mwm --dbmem --data_dir ${data_dir} \
 --tax_scope 5338,155619,5204,4751,33154,2759,1 --cpu ${t}

conda deactivate
