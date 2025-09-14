#!/bin/bash

# ============================
# AMR Detection Pipeline
# Author: Opeoluwa
# Date: September 2025
# ============================

# Step 1: Download raw paired-end reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/059/SRR13554759/SRR13554759_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/059/SRR13554759/SRR13554759_2.fastq.gz

# Step 2: Run FastQC for quality control
mkdir -p qc_reports
fastqc SRR13554759_1.fastq.gz SRR13554759_2.fastq.gz -o qc_reports

# Step 3: Trim reads using fastp
mkdir -p trim
fastp \
  -i SRR13554759_1.fastq.gz \
  -I SRR13554759_2.fastq.gz \
  -o trim/SRR13554759_1.trim.fastq.gz \
  -O trim/SRR13554759_2.trim.fastq.gz \
  -h trim/fastp_report.html \
  -j trim/fastp_report.json \
  --detect_adapter_for_pe \
  --thread 4

# Step 4: Assemble genome using SPAdes with corrected odd k-mer values
mkdir -p assembly
spades.py \
  -1 trim/SRR13554759_1.trim.fastq.gz \
  -2 trim/SRR13554759_2.trim.fastq.gz \
  -o assembly \
  -k 41,85

# Step 5: Detect AMR genes using Abricate
mkdir -p abricate_results
abricate --db resfinder assembly/contigs.fasta > abricate_results/resfinder_output.txt

# Step 6: Generate a clean resistance report
echo "Detected Antimicrobial Resistance Genes:" > resistance_report.txt
grep -v "^#" abricate_results/resfinder_output.txt | awk '{print $1, $5, $6, $7, $8}' >> resistance_report.txt

# Optional: Summarize results
abricate --summary abricate_results/resfinder_output.txt >> resistance_report.txt

echo 







Suggested Improvements: The script is mostly comprehensive and correctly follows the outlined steps. However, it would be beneficial to include error handling for each step, such as checking for the existence of required tools and successful completion of each command. Additionally, it might be helpful to add comments or documentation for users who may not be familiar with each bioinformatics tool.
