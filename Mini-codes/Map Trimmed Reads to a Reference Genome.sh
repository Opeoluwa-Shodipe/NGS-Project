#!/bin/bash

# Create necessary directories
mkdir -p repaired alignment_map

# Index the reference genome
bwa index references/reference.fasta

# Define array of sample names
samples=(ACBarrie Alsene Baxter Chara Drysdale)

# Loop through each sample
for SAMPLE in "${samples[@]}"; do
    echo "Repairing reads for $SAMPLE..."

    # Run repair.sh to fix paired-end reads and handle singletons
    repair.sh \
        in1=qc_reads/${SAMPLE}_R1.fastq.gz \
        in2=qc_reads/${SAMPLE}_R2.fastq.gz \
        out1=repaired/${SAMPLE}_R1.fastq.gz \
        out2=repaired/${SAMPLE}_R2.fastq.gz \
        outs=repaired/${SAMPLE}_singletons.fastq.gz \
        overwrite=true

    echo "Mapping $SAMPLE to reference genome..."

    # Map repaired reads using bwa mem
    bwa mem references/reference.fasta \
        repaired/${SAMPLE}_R1.fastq.gz \
        repaired/${SAMPLE}_R2.fastq.gz \
        > alignment_map/${SAMPLE}.sam

    echo "Converting SAM to BAM for $SAMPLE..."

    # Convert SAM to BAM
    samtools view -b alignment_map/${SAMPLE}.sam > alignment_map/${SAMPLE}.bam

    # Optional: remove SAM to save space
    rm alignment_map/${SAMPLE}.sam

    echo "$SAMPLE processing complete."
done

echo "All samples mapped and converted to BAM format."






















#!/bin/bash

SAMPLES=(
  "ACBarrie"
  "Alsene"
  "Baxter"
  "Chara"
  "Drysdale"
)

bwa index references/reference.fasta
mkdir repaired
mkdir alignment_map

for SAMPLE in "${SAMPLES[@]}"; do

    repair.sh in1="trimmed_reads/${SAMPLE}_R1.fastq.gz" in2="trimmed_reads/${SAMPLE}_R2.fastq.gz" out1="repaired/${SAMPLE}_R1_rep.fastq.gz" out2="repaired/${SAMPLE}_R2_rep.fastq.gz" outsingle="repaired/${SAMPLE}_single.fq"
    echo $PWD
    bwa mem -t 1 \
    references/reference.fasta \
    "repaired/${SAMPLE}_R1_rep.fastq.gz" "repaired/${SAMPLE}_R2_rep.fastq.gz" \
  | samtools view -b \
  > "alignment_map/${SAMPLE}.bam"
done
