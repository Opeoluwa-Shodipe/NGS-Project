#! /bin/bash
# This script is used to perform quality control for our RNA-seq data for C. elegans
# create directories for output
mkdir -p qc

for filename in  clean_files/*.fq; do
	fastqc -o qc/ $filename
done

#!/bin/bash
#trimming files
mkdir -p trim

for filename in clean_files/*.fq; do
    base=$(basename "$filename" .fq)
    echo fastp -i "$filename" -o "trim/${base}_trim.fq" -h "trim/${base}_report.html"
done

#advisable to do qc on the trimmed files
