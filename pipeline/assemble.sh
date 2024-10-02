#!/bin/bash

# Parameterized script for filtering Nanopore long reads, assembling with Flye, and polishing with Illumina short reads
# Usage: ./filter_assemble.sh <ACCESSION_LONG> <ACCESSION_SHORT> <THREADS> <OUTPUT_DIR>

# Check if correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <ACCESSION_LONG> <ACCESSION_SHORT> <THREADS> <OUTPUT_DIR>"
    exit 1
fi

# Variables from command-line arguments
ACCESSION_LONG=$1      # Nanopore long reads SRA accession number
ACCESSION_SHORT=$2     # Illumina short reads SRA accession number
THREADS=$3             # Number of threads for filtering and assembly
OUTPUT_DIR=$4          # Output directory for the assembly

# Step 0: Check dependencies
echo "Checking dependencies..."

REQUIRED_TOOLS=("NanoFilt" "flye" "fastp" "bwa" "samtools" "pilon" "java")

for TOOL in "${REQUIRED_TOOLS[@]}"
do
    if ! command -v $TOOL &> /dev/null
    then
        echo "$TOOL could not be found. Please install $TOOL."
        exit 1
    fi
done

echo "All required tools are installed."

# Step 1: Perform QC filtering on Nanopore long reads
echo "Filtering Nanopore long reads with NanoFilt..."
zcat "${ACCESSION_LONG}.fastq.gz" | NanoFilt -q 7 -l 1000 | gzip > "${ACCESSION_LONG}_filtered.fastq.gz"

# Step 2: Perform QC filtering on Illumina short reads
echo "Filtering Illumina short reads with fastp..."
fastp -i "${ACCESSION_SHORT}_1.fastq.gz" -I "${ACCESSION_SHORT}_2.fastq.gz" \
      -o "${ACCESSION_SHORT}_1_filtered.fastq.gz" -O "${ACCESSION_SHORT}_2_filtered.fastq.gz" \
      --thread $THREADS --qualified_quality_phred 20 --length_required 50

# Step 3: Assemble Nanopore long reads with Flye
echo "Running Flye assembly..."
flye --nano-raw "${ACCESSION_LONG}_filtered.fastq.gz" --out-dir $OUTPUT_DIR --threads $THREADS

# Step 4: Polishing assembly with Illumina short reads

# Step 4.1: Index the assembly
echo "Indexing assembly..."
ASSEMBLY="${OUTPUT_DIR}/assembly.fasta"
bwa index $ASSEMBLY

# Step 4.2: Align Illumina short reads to the assembly
echo "Aligning short reads to the assembly..."
bwa mem -t $THREADS $ASSEMBLY "${ACCESSION_SHORT}_1_filtered.fastq.gz" "${ACCESSION_SHORT}_2_filtered.fastq.gz" | \
samtools view -@ $THREADS -bS - | samtools sort -@ $THREADS -o alignment_sorted.bam

# Step 4.3: Index the BAM file
echo "Indexing BAM file..."
samtools index alignment_sorted.bam

# Step 4.4: Run Pilon for polishing
echo "Running Pilon for polishing..."
pilon --genome $ASSEMBLY --frags alignment_sorted.bam --output pilon_polished --threads $THREADS

echo "Assembly and polishing completed. Results are available in the directory: $OUTPUT_DIR"

