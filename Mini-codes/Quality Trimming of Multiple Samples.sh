#!/bin/bash

# Create output directory
mkdir -p qc_reads

# Define array of sample names
samples=(ACBarrie Alsen Baxter Chara Drysdale)

# Loop through each sample and run fastp
for SAMPLE in "${samples[@]}"; do
    echo "Processing $SAMPLE..."

    fastp \
        -i ${SAMPLE}_R1.fastq.gz \
        -I ${SAMPLE}_R2.fastq.gz \
        -o qc_reads/${SAMPLE}_R1.fastq.gz \
        -O qc_reads/${SAMPLE}_R2.fastq.gz \
        -h qc_reads/${SAMPLE}_fastp.html \
        -j qc_reads/${SAMPLE}_fastp.json \
        --detect_adapter_for_pe \
        --thread 4

    echo "$SAMPLE done."
done

















#!/bin/bash
mkdir qc_reads
SAMPLES=(
  "ACBarrie"
  "Alsen"
  "Baxter"
  "Chara"
  "Drysdale"
)

for SAMPLE in "${SAMPLES[@]}"; do

  fastp \
    -i "$PWD/${SAMPLE}_R1.fastq.gz" \
    -I "$PWD/${SAMPLE}_R2.fastq.gz" \
    -o "qc_reads/${SAMPLE}_R1.fastq.gz" \
    -O "qc_reads/${SAMPLE}_R2.fastq.gz" \
    --html "qc_reads/${SAMPLE}_fastp.html" 
done

echo "All samples processed."
