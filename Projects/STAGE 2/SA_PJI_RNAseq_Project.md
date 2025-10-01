# Transcriptomic Profiling of *Staphylococcus aureus* — Acute vs Chronic PJI (RNA‑seq)

**Deliverable:** reproducible pipeline (Bash + R) and a concise analysis/report template suitable for GitHub submission.  
**Dataset:** PRJNA867318 (8 samples: SRR20959676–SRR20959683).

---

## Overview & Goals
This project is comparing transcriptomes of *S. aureus* isolates from **acute** versus **chronic** periprosthetic joint infections (PJI) to identify expression programs linked to virulence (acute) and persistence/biofilm (chronic). Goals:
- QC, trimming, alignment to *S. aureus* reference, count matrix generation.
- Differential expression (DESeq2) acute vs chronic.
- Functional enrichment (GO/KEGG) and pathway mapping.
- Visualization: PCA, volcano, heatmap, pathway tables.
- Short clinical microbiology report summarizing findings and implications.

---

## Requirements (recommended)
Install via conda (bioconda / conda-forge):
```bash
conda create -n sa_rnaseq -c conda-forge -c bioconda \
  fastp sra-tools bowtie2 samtools subread multiqc \
  python pandas r-base r-essentials bioconductor-deseq2 r-ggplot2 r-pheatmap r-clusterProfiler r-enrichplot r-EnhancedVolcano
```
(Alternatively adjust to use `hisat2` instead of `bowtie2`.)

---

## Directory layout (suggested)
```
sa_rnaseq_project/
  data/raw/         # raw FASTQ from SRA
  data/trimmed/
  results/bam/
  results/counts/
  results/figures/
  results/reports/
  scripts/
```
Create it:
```bash
mkdir -p data/raw data/trimmed results/{bam,counts,figures,reports} scripts
```

---

## Sample sheet (samples.tsv)
Create `samples.tsv` in project root with columns: `sample,SRR,condition`

```
sample,SRR,condition
chronic_1,SRR20959676,chronic
chronic_2,SRR20959677,chronic
chronic_3,SRR20959678,chronic
chronic_4,SRR20959679,chronic
acute_1,SRR20959680,acute
acute_2,SRR20959681,acute
acute_3,SRR20959682,acute
acute_4,SRR20959683,acute
```

---

## 1) Download raw FASTQ (SRA -> FASTQ)
Use `prefetch` + `fasterq-dump`:
```bash
while IFS=, read -r sample srr condition; do
  [[ "$sample" == "sample" ]] && continue
  echo "Fetching $srr -> data/raw/${sample}.fastq.gz"
  prefetch $srr -O sra_raw
  fasterq-dump sra_raw/$srr/$srr.sra -o data/raw/${sample}.fastq --split-files -e 4
  gzip -f data/raw/${sample}_1.fastq data/raw/${sample}_2.fastq || true
done < samples.tsv
```

---

## 2) Reference genome & annotation
Download a *Staphylococcus aureus* reference (RefSeq). Example: **NCTC 8325 (NC_007795.1)**. Place `ref.fasta` and annotation GFF/GTF `ref.gff` in `data/reference/`.

Index for bowtie2 and samtools:
```bash
bowtie2-build data/reference/ref.fasta data/reference/ref_index
samtools faidx data/reference/ref.fasta
```

---

## 3) Processing pipeline (Bash): `scripts/run_sa_rnaseq.sh`
```bash
#!/usr/bin/env bash
set -euo pipefail
THREADS=8
SAMPFILE="samples.tsv"
REF="data/reference/ref.fasta"
INDEX="data/reference/ref_index"
RAW="data/raw"
TRIM="data/trimmed"
BAMDIR="results/bam"
CTSDIR="results/counts"
FIGDIR="results/figures"
mkdir -p "$TRIM" "$BAMDIR" "$CTSDIR" "$FIGDIR"

while IFS=, read -r sample srr cond; do
  [[ "$sample" == "sample" ]] && continue
  r1="${RAW}/${sample}_1.fastq.gz"
  r2="${RAW}/${sample}_2.fastq.gz"
  out1="${TRIM}/${sample}_1.trim.fastq.gz"
  out2="${TRIM}/${sample}_2.trim.fastq.gz"
  fastp -i "$r1" -I "$r2" -o "$out1" -O "$out2" -w $THREADS
  bowtie2 -x "$INDEX" -1 "$out1" -2 "$out2" -p $THREADS | samtools view -b -@ $THREADS -o "${BAMDIR}/${sample}.bam" -
  samtools sort -@ $THREADS -o "${BAMDIR}/${sample}.sorted.bam" "${BAMDIR}/${sample}.bam"
  samtools index "${BAMDIR}/${sample}.sorted.bam"
done < "$SAMPFILE"

featureCounts -T $THREADS -a data/reference/ref.gff -o "${CTSDIR}/featureCounts.txt" ${BAMDIR}/*.sorted.bam
```

---

## 4) Differential Expression (R): `scripts/deseq2_sa.R`
```r
library(DESeq2); library(pheatmap); library(EnhancedVolcano)
samples <- read.csv("samples.tsv", stringsAsFactors = FALSE)
fc <- read.table("results/counts/featureCounts.txt", header = TRUE, sep = "\t", comment.char = "#")
counts <- as.matrix(fc[ , 7:ncol(fc)])
rownames(counts) <- fc$Geneid
samples <- samples[samples$sample %in% colnames(counts), ]
dds <- DESeqDataSetFromMatrix(countData = counts, colData = samples, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","acute","chronic"))
write.csv(as.data.frame(res), "results/counts/deseq2_results.csv")
```

---

## Outputs
- **Counts:** `results/counts/featureCounts.txt`
- **DEGs:** `results/counts/deseq2_results.csv`
- **Figures:** PCA, volcano, heatmap
- **Report:** `results/reports/clinical_micro_report.md`

---

## References
- DESeq2 manual: https://bioconductor.org/packages/release/bioc/html/DESeq2.html  
- featureCounts: Liao et al., Bioinformatics 2014  
- clusterProfiler: Yu et al., OMICS 2012  
- SRA PRJNA867318

# -------- Professional Profile --------
- GitHub: https://github.com/Opeoluwa-Shodipe/NGS-Project-HACKBIO-/blob/main/Projects/STAGE%202/SA_PJI_RNAseq_Project.md
- LinkedIn: Ensure video is uploaded and viewed public."

