# Clinical Case Study: Whole-Genome Sequencing (CFTR Analysis)

## Background

A 6-year-old boy presents with chronic cough, recurrent lung infections, and poor weight gain despite adequate nutrition. His physician suspects **cystic fibrosis (CF)**. A sweat chloride test yields borderline results (45 mmol/L). To confirm diagnosis and identify the causative genetic variant(s), the clinical team orders **whole-genome sequencing (WGS)**.  

The father’s genome has also been sequenced for inheritance analysis, while the mother’s genomic data are not available.

---

## Objectives

- Perform quality control and alignment of raw WGS reads.  
- Call and filter variants in both father and child.  
- Focus on the **CFTR locus (7q31.2)**.  
- Compare inheritance patterns between child and father.  
- Annotate variants with VEP and ClinVar.  
- Interpret whether pathogenic CFTR mutations (ΔF508, G542X, W1282X, etc.) are present.  

---

## Pipeline Overview

This pipeline was implemented in **Bash** with standard bioinformatics tools. Below is the step-by-step workflow.  

### 1. Tool Check

Ensure required tools are installed: `fastp`, `bwa`, `samtools`, `picard`, `bcftools`, `bgzip`, `tabix`, `vep`.

```bash
check_tools() {
  for t in fastp bwa samtools picard bcftools bgzip tabix vep; do
    command -v "$t" >/dev/null 2>&1 || { echo "Error: $t not installed"; exit 1; }
  done
}
```

### 2. Quality Control and Trimming

Run `fastp` to trim adapters and generate QC reports.

```bash
fastp -i child_1.fastq.gz -I child_2.fastq.gz       -o results/child_1.trim.fastq.gz -O results/child_2.trim.fastq.gz       -w 8 -h results/child_fastp.html -j results/child_fastp.json
```

### 3. Reference Indexing

Index the reference genome for alignment.

```bash
bwa index GRCh38.primary_assembly.genome.fa
samtools faidx GRCh38.primary_assembly.genome.fa
picard CreateSequenceDictionary R=GRCh38.primary_assembly.genome.fa O=GRCh38.dict
```

### 4. Alignment

Align trimmed reads to the reference.

```bash
bwa mem -t 8 GRCh38.fa results/child_1.trim.fastq.gz results/child_2.trim.fastq.gz | samtools view -b -o results/child.unsorted.bam -
```

### 5. Sorting, Deduplication, and Indexing

```bash
samtools sort -o results/child.sorted.bam results/child.unsorted.bam
picard MarkDuplicates I=results/child.sorted.bam O=results/child.dedup.bam M=results/child.dup.metrics.txt
samtools index results/child.dedup.bam
```

### 6. Coverage Check

Assess mean depth across the **CFTR locus**.

```bash
samtools depth -r 7:117120016-117308718 results/child.dedup.bam | awk '{sum+=$3; cnt++} END{print "child_mean_dp", sum/cnt}'
```

### 7. Variant Calling

Joint variant calling for child and father.

```bash
bcftools mpileup -f GRCh38.fa -r 7:117120016-117308718 results/child.dedup.bam results/father.dedup.bam -Ou | bcftools call -mv -Oz -o results/joint.cftr.vcf.gz
```

### 8. Variant Filtering

Apply basic filters (QUAL ≥ 20, DP ≥ 10).

```bash
bcftools view -i 'QUAL>=20 && INFO/DP>=10' results/joint.cftr.vcf.gz -Oz -o results/joint.cftr.filtered.vcf.gz
```

### 9. Annotation

Use **VEP** for variant annotation.

```bash
vep -i results/joint.cftr.filtered.vcf.gz --cache --assembly GRCh38 --offline     --fasta GRCh38.fa --everything --vcf -o results/joint.cftr.vep.vcf.gz --compress_output bgzip --thread 8
```

### 10. Summary Tables

Extract key fields for comparison.

```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ\t%INFO/CLIN_SIG\t%INFO/AF\t[%SAMPLE=%GT:%DP]\n' results/joint.cftr.vep.vcf.gz > results/joint_cftr_annotation.tsv
```

---

## Expected Results

- **QC Reports**: Adapter content, quality trimming statistics.  
- **Coverage**: Adequate mean depth across CFTR (>30x recommended).  
- **Variants**: Candidate CFTR mutations extracted in annotated VCF.  
- **Inheritance**: Comparison of father vs child genotypes for inheritance pattern.  
- **Annotation**: Pathogenic mutations identified via ClinVar.  

---

## Clinical Interpretation

The presence of **two pathogenic CFTR variants in trans** would confirm cystic fibrosis. If only one variant is identified, further maternal testing is required to confirm inheritance. Coverage analysis is critical to ensure no missed variants.

---

## Professional Profile

- **GitHub**: [Project Repository](https://github.com/Opeoluwa-Shodipe/NGS-Project-HACKBIO-/blob/main/Projects/STAGE%201/CFTR_WGS_CaseStudy.md)  
- **LinkedIn**: A short video/post will be uploaded explaining the project in lay terms.  

---

## Conclusion

This case study demonstrates the workflow for analyzing WGS data to identify pathogenic variants in **CFTR**. The pipeline integrates QC, alignment, variant calling, annotation, and inheritance analysis, providing clinically relevant insights into suspected cystic fibrosis cases.
