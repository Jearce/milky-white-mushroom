# milkywhite-mushroom

Scripts to assist with the genome assembly and annotation of Calocybe indica

## Scripts

### download_fungal_genomes.sh

#### Python dependencies

- biopython
- pandas

#### Set up

Clone repo:

```
git clone https://github.com/Jearce/milkywhite-mushroom.git
cd milkywhite-mushroom
```

After repo is cloned. Add the repo directory to your PATH.

Then make download_fungal_genomes.sh and extract_and_format_fungi_names.py executeble:

```
chmod +x download_fungal_genomes.sh
chmod +x extract_and_format_fungi_names.py
```

#### Usage

download_fungal_genome.sh taks a single csv file from [EnsemblFungi](http://fungi.ensembl.org/species.html)
and will use the name of the organisms to download all CDS in fasta format.

Below is the help message:
```
usage: ./download_fungal_genomes.sh <csv file>
 Any failed downloads are logged to fungal.downloads.log file.
```

