#!/bin/bash

# Download the sample FASTQ file
echo "Downloading FASTQ file..."
wget -O sample_reads.fastq.gz https://github.com/rieseberglab/fastq-examples/raw/refs/heads/master/data/HI.4019.002.index_7.ANN0831_R1.fastq.gz

# Install FastQC if not already installed
if ! command -v fastqc &> /dev/null; then
    echo "FastQC not found. Installing..."
    sudo apt-get update
    sudo apt-get install fastqc -y
fi

# Run FastQC on the downloaded file
echo "Running FastQC..."
fastqc sample_reads.fastq.gz

# Locate and print report file names
echo "FastQC report files generated:"
ls -1 sample_reads_fastqc.*

# User guidance to view the report
echo ""
echo "To view the FastQC report, download or open the file:"
echo "  sample_reads_fastqc.html"
echo "from your current directory in a web browser."
