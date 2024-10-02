#!/bin/bash

# Parameterized script to download Nanopore long reads and Illumina short reads
# Usage: ./download_data.sh <ACCESSION_LONG> <ACCESSION_SHORT> <THREADS>

# Check if correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <ACCESSION_LONG> <ACCESSION_SHORT> <THREADS>"
    exit 1
fi

# Variables from command-line arguments
ACCESSION_LONG=$1     # Nanopore long reads SRA accession number
ACCESSION_SHORT=$2    # Illumina short reads SRA accession number
THREADS=$3            # Number of threads for downloading

# Step 1: Download Nanopore long reads
echo "Downloading Nanopore long read data for accession $ACCESSION_LONG..."
fasterq-dump --threads $THREADS $ACCESSION_LONG

# Compress the long read FASTQ file to save space
echo "Compressing long read FASTQ file..."
gzip "${ACCESSION_LONG}.fastq"

# Step 2: Download Illumina short reads
echo "Downloading Illumina short read data for accession $ACCESSION_SHORT..."
fasterq-dump --split-files --threads $THREADS $ACCESSION_SHORT

# Compress the short read FASTQ files to save space
echo "Compressing short read FASTQ files..."
gzip "${ACCESSION_SHORT}_1.fastq"
gzip "${ACCESSION_SHORT}_2.fastq"

echo "Download and compression completed."

