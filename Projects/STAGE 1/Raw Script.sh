#!/bin/bash

# =====================================================
# Project Script: Bash Basics + Bioinformatics Software
# Author: Opeoluwa
# =====================================================

# run_cftr_pipeline.sh
set -euo pipefail

# --- Global Config ---
THREADS=8
OUTDIR=results
mkdir -p "$OUTDIR"

# --- User Variables ---
REF=/path/to/GRCh38.primary_assembly.genome.fa
CHILD_1=child_1.fastq.gz
CHILD_2=child_2.fastq.gz
FATHER_1=father_1.fastq.gz
FATHER_2=father_2.fastq.gz
CFTR_REGION="7:117120016-117308718"  # GRCh38 coordinates for CFTR

# --- Tool Check ---
check_tools() {
  echo "Checking required tools..."
  for t in fastp bwa samtools picard bcftools bgzip tabix vep; do
    command -v "$t" >/dev/null 2>&1 || { echo "Error: $t is required but not installed."; exit 1; }
  done
}

# --- Quality Control & Trimming ---
run_fastp() {
  echo "Running fastp for child..."
  fastp -i "$CHILD_1" -I "$CHILD_2" \
        -o "$OUTDIR/child_1.trim.fastq.gz" -O "$OUTDIR/child_2.trim.fastq.gz" \
        -w "$THREADS" -h "$OUTDIR/child_fastp.html" -j "$OUTDIR/child_fastp.json"

  echo "Running fastp for father..."
  fastp -i "$FATHER_1" -I "$FATHER_2" \
        -o "$OUTDIR/father_1.trim.fastq.gz" -O "$OUTDIR/father_2.trim.fastq.gz" \
        -w "$THREADS" -h "$OUTDIR/father_fastp.html" -j "$OUTDIR/father_fastp.json"
}

# --- Reference Indexing ---
index_reference() {
  echo "Indexing reference genome if needed..."
  [ -f "${REF}.bwt" ] || bwa index "$REF"
  [ -f "${REF}.fai" ] || samtools faidx "$REF"
  [ -f "${REF%.fa}.dict" ] || java -jar "$(which picard)" CreateSequenceDictionary R="$REF" O="${REF%.fa}.dict"
}

#1 To repair trimmed reads before aligning to reference genome
repair.sh in1-trimmed_reads/child_1.fastq.gz in2-trimmed_reads/child_2.fastq.gz out1-child_1_rep.fastq.gz  out2-child_2_rep.fastq.gz outsingle-single.fq
#2 To index
bwa index reference/hg38/hg38.fasta

# --- Alignment ---
align_reads() {
  echo "Aligning child reads..."
  bwa mem -t "$THREADS" "$REF" "$OUTDIR/child_1.trim.fastq.gz" "$OUTDIR/child_2.trim.fastq.gz" \
    | samtools view -b -@ "$THREADS" -o "$OUTDIR/child.unsorted.bam" -

  echo "Aligning father reads..."
  bwa mem -t "$THREADS" "$REF" "$OUTDIR/father_1.trim.fastq.gz" "$OUTDIR/father_2.trim.fastq.gz" \
    | samtools view -b -@ "$THREADS" -o "$OUTDIR/father.unsorted.bam" -
}

# Alignment code
bwa mem reference/hg38/hg38.fasta child_1_rep.fastq.gz child_2_rep.fastq.gz > alignment/child.sam

# --- Sorting, Deduplication, Indexing ---
process_bams() {
  echo "Sorting and marking duplicates..."
  samtools sort -@ "$THREADS" -o "$OUTDIR/child.sorted.bam" "$OUTDIR/child.unsorted.bam"
  samtools sort -@ "$THREADS" -o "$OUTDIR/father.sorted.bam" "$OUTDIR/father.unsorted.bam"

  picard MarkDuplicates I="$OUTDIR/child.sorted.bam" O="$OUTDIR/child.dedup.bam" M="$OUTDIR/child.dup.metrics.txt"
  picard MarkDuplicates I="$OUTDIR/father.sorted.bam" O="$OUTDIR/father.dedup.bam" M="$OUTDIR/father.dup.metrics.txt"

  samtools index "$OUTDIR/child.dedup.bam"
  samtools index "$OUTDIR/father.dedup.bam"
}

# --- Coverage Check ---
check_coverage() {
  echo "Checking coverage over CFTR region..."
  samtools depth -r "$CFTR_REGION" "$OUTDIR/child.dedup.bam" | awk '{sum+=$3; cnt++} END{print "child_mean_dp", sum/cnt}'
  samtools depth -r "$CFTR_REGION" "$OUTDIR/father.dedup.bam" | awk '{sum+=$3; cnt++} END{print "father_mean_dp", sum/cnt}'
}

# --- Variant Calling ---
call_variants() {
  echo "Calling variants in CFTR region..."
  bcftools mpileup -f "$REF" -r "$CFTR_REGION" "$OUTDIR/child.dedup.bam" "$OUTDIR/father.dedup.bam" -Ou \
    | bcftools call -mv -Oz -o "$OUTDIR/joint.cftr.vcf.gz"
  bcftools index "$OUTDIR/joint.cftr.vcf.gz"
}

# --- Variant Filtering ---
filter_variants() {
  echo "Filtering variants..."
  bcftools view -i 'QUAL>=20 && INFO/DP>=10' "$OUTDIR/joint.cftr.vcf.gz" -Oz -o "$OUTDIR/joint.cftr.filtered.vcf.gz"
  bcftools index "$OUTDIR/joint.cftr.filtered.vcf.gz"
}

# --- Annotation ---
annotate_variants() {
  echo "Annotating variants with VEP..."
  vep -i "$OUTDIR/joint.cftr.filtered.vcf.gz" \
      --cache --assembly GRCh38 --offline \
      --fasta "$REF" \
      --everything \
      --vcf -o "$OUTDIR/joint.cftr.vep.vcf.gz" --compress_output bgzip --thread "$THREADS"
  tabix -p vcf "$OUTDIR/joint.cftr.vep.vcf.gz"

  # Optional: Add ClinVar plugin if configured
  # --plugin ClinVar,"/path/to/clinvar.vcf.gz"
}

# --- Summary Tables ---
generate_tables() {
  echo "Generating summary tables..."
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ\t%INFO/CLIN_SIG\t%INFO/AF\t[%SAMPLE=%GT:%DP]\n' \
    "$OUTDIR/joint.cftr.vep.vcf.gz" > "$OUTDIR/joint_cftr_annotation.tsv"

  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE=%GT\t]\n' \
    "$OUTDIR/joint.cftr.vep.vcf.gz" > "$OUTDIR/compare_genotypes.tsv"
}

# --- Run Pipeline ---
main() {
  check_tools
  run_fastp
  index_reference
  align_reads
  process_bams
  check_coverage
  call_variants
  filter_variants
  annotate_variants
  generate_tables
  echo "Pipeline complete. Outputs saved in $OUTDIR"
}

main



# -------- Professional Profile --------
echo "Professional Profile"
echo "1. GitHub: https://github.com/Opeoluwa-Shodipe/NGS-Project-HACKBIO-/blob/main/Projects/STAGE%201/CFTR_WGS_CaseStudy.md"
echo "2. LinkedIn: Ensure video is uploaded and public."





