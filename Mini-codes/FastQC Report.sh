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




# 1. Open the FastQC .html report from Task 1
scp -r a_adegite@135.181.163.242:~/Eadencre8ives/FastQC/sample_reads_fastqc.html ./

# 2. Base 1-109 has high quality and are mostly in the green zone
     Base 110-149 has declined in quality but in the warning zone
     Base 150 Above has a poor quality and are in the red zone
     No warnings about adapter sequences in "Adapter Content"
     Yes Overrepresented sequences were flagged

The sequencing data in sample_reads.fastq.gz is of generally high quality, with no reads flagged as poor and a consistent read length of 150 bp. The average per-read Phred score is above 30, indicating reliable base calling across most of the dataset. However, the "Per base sequence quality" module reveals a noticeable decline in quality toward the final 15–20 bases, where scores drop below 20, suggesting increased error rates at the read ends.
Additionally, the presence of overrepresented adapter sequences—such as TruSeq and Illumina multiplexing primers—indicates that adapter contamination is present and may interfere with downstream analyses like alignment or variant calling. These issues can lead to misalignments, inflated duplication rates, and biased base composition metrics. To mitigate these effects, trimming low-quality ends and removing adapter sequences using tools like Trimmomatic, Cutadapt, or fastp is strongly recommended before proceeding with further bioinformatics workflows.







Suggested Improvements: The script currently assumes a Linux environment with 'apt-get' as the package manager. To improve portability, consider checking or providing installation support for other operating systems or package managers.
