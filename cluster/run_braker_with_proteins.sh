#SBATCH -J braker
#SBATCH -o braker.o%j
#SBATCH -c 20
#SBATCH --mem=15G
#SBATCH -t 05:00:00
#SBATCH --mail-type=END,FAIL

set -e
set -u
set -o pipefail

#set local anaconda3 env located in project dir
set +eu
source ~/.bashrc 

if [ $# -ne 2 ]
then
	echo "Incorrect number of arguments"
	echo "run_braker.sh <softmasked_contigs> <proteins>"
	exit
fi

export contigs=$1; export proteins=$2;

if ! egrep -e "[atgc]" ${contigs}
then
	echo "${contigs} must be softmasked"
	exit
fi
conda activate braker-env

prothint.py ${contigs} ${proteins} --threads 20

export hints="prothint_augustus.gff"

if [ -f ${hints} ]
then
	braker.pl --genome=${contigs} --hints=${hints} --softmasking --fungus --cores 20
else
	"${hints} does not exist"
fi

conda deactivate
