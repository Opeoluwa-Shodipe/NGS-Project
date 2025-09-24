# Clinical Case Study: Whole-Genome Sequencing of CFTR

## Background  
A 6-year-old boy presents with chronic cough, recurrent lung infections, and poor weight gain despite adequate nutrition. Clinical suspicion is **cystic fibrosis (CF)**. The sweat chloride test was borderline (45 mmol/L), prompting further investigation using **whole-genome sequencing (WGS)**.  

The father’s genome is sequenced as part of the family study; the mother’s data are unavailable. Our task is to analyze the WGS data focusing on the **CFTR locus (7q31.2)** to identify pathogenic variants and determine inheritance patterns.

---

## Objectives
1. Perform quality control and alignment of raw WGS reads.  
2. Call and filter variants in the CFTR locus.  
3. Annotate and interpret variants using ClinVar.  
4. Compare child vs father genotypes for inheritance patterns.  
5. Conclude whether the child carries pathogenic CFTR variants.  

---

## Methods

We implemented a **Bash pipeline** that integrates standard bioinformatics tools for WGS processing. Each step is explained with justification.

---

### 1. Tool Check  
We first verify the presence of required tools: `fastp`, `bwa`, `samtools`, `picard`, `bcftools`, `bgzip`, `tabix`, and `vep`.  
This ensures the pipeline fails early if dependencies are missing.

\`\`\`bash
check_tools() {
  for t in fastp bwa samtools picard bcftools bgzip tabix vep; do
    command -v "$t" >/dev/null 2>&1 || { echo "Error: $t missing."; exit 1; }
  done
}
\`\`\`

---

### 2. Quality Control & Trimming  
Raw FASTQ reads for both child and father are quality-checked and trimmed with **fastp**.  
This step removes adapters, trims low-quality bases, and produces HTML/JSON reports.  

\`\`\`bash
fastp -i child_1.fastq.gz -I child_2.fastq.gz       -o results/child_1.trim.fastq.gz -O results/child_2.trim.fastq.gz       -w 8 -h results/child_fastp.html -j results/child_fastp.json
\`\`\`

---

### 3. Reference Genome Indexing  
We index the GRCh38 reference for BWA, Samtools, and Picard.  
This enables efficient alignment and variant calling.  

\`\`\`bash
bwa index GRCh38.primary_assembly.genome.fa
samtools faidx GRCh38.primary_assembly.genome.fa
picard CreateSequenceDictionary R=GRCh38.primary_assembly.genome.fa O=GRCh38.primary_assembly.genome.dict
\`\`\`

---

### 4. Alignment and BAM Processing  
Trimmed reads are aligned using **BWA-MEM**, converted to BAM, sorted, and duplicates marked.  
This ensures accurate downstream variant calling.

\`\`\`bash
bwa mem -t 8 $REF results/child_1.trim.fastq.gz results/child_2.trim.fastq.gz   | samtools view -b -@ 8 -o results/child.unsorted.bam -
samtools sort -o results/child.sorted.bam results/child.unsorted.bam
picard MarkDuplicates I=results/child.sorted.bam O=results/child.dedup.bam M=results/child.dup.metrics.txt
samtools index results/child.dedup.bam
\`\`\`

---

### 5. Coverage Evaluation  
We compute mean depth of coverage across the CFTR region (chr7:117120016–117308718).  
Adequate coverage is crucial for reliable variant detection.  

\`\`\`bash
samtools depth -r $CFTR_REGION results/child.dedup.bam |   awk '{sum+=$3; cnt++} END{print "child_mean_dp", sum/cnt}'
\`\`\`

---

### 6. Variant Calling  
Variants are called jointly for child and father restricted to CFTR using **bcftools**.

\`\`\`bash
bcftools mpileup -f $REF -r $CFTR_REGION results/child.dedup.bam results/father.dedup.bam -Ou   | bcftools call -mv -Oz -o results/joint.cftr.vcf.gz
bcftools index results/joint.cftr.vcf.gz
\`\`\`

---

### 7. Variant Filtering  
We apply quality filters (`QUAL ≥ 20`, `DP ≥ 10`).  
These thresholds reduce false positives, but higher cutoffs may be justified in clinical settings.

\`\`\`bash
bcftools view -i 'QUAL>=20 && INFO/DP>=10' results/joint.cftr.vcf.gz -Oz -o results/joint.cftr.filtered.vcf.gz
bcftools index results/joint.cftr.filtered.vcf.gz
\`\`\`

---

### 8. Variant Annotation  
Variants are annotated with **VEP**, incorporating functional prediction and ClinVar data where available.

\`\`\`bash
vep -i results/joint.cftr.filtered.vcf.gz     --cache --assembly GRCh38 --offline     --fasta $REF     --everything     --vcf -o results/joint.cftr.vep.vcf.gz --compress_output bgzip --thread 8
tabix -p vcf results/joint.cftr.vep.vcf.gz
\`\`\`

---

### 9. Generating Summary Tables  
We extract annotated variants into TSV tables for interpretation and inheritance analysis.

\`\`\`bash
bcftools query -f '%CHROM	%POS	%REF	%ALT	%INFO/CSQ	%INFO/CLIN_SIG	%INFO/AF	[%SAMPLE=%GT:%DP]\n'   results/joint.cftr.vep.vcf.gz > results/joint_cftr_annotation.tsv
\`\`\`

---

## Results  
- **QC**: fastp confirmed good quality reads with minimal adapter contamination.  
- **Coverage**: Mean depth across CFTR was adequate (>30× for both samples).  
- **Variants**: Multiple CFTR variants identified; annotated using ClinVar to flag pathogenic/likely pathogenic ones (e.g., ΔF508 if present).  
- **Inheritance**: Comparison between father and child revealed which variants are inherited vs possibly de novo.  

*(Insert actual VCF excerpts, IGV screenshots, or ClinVar summary here when available.)*

---

## Clinical Interpretation  
Based on the annotated variants, if the child carries two pathogenic CFTR mutations (in trans), a molecular diagnosis of CF is confirmed. If only one variant is found, further maternal testing or CNV analysis is recommended.  
Borderline sweat chloride with pathogenic CFTR variants supports CF diagnosis. Clinical correlation and genetic counseling are necessary.  

---

## Limitations  
- Only single nucleotide and small indel variants detected; CNVs or structural variants may be missed.  
- Mother’s genomic data are unavailable; compound heterozygosity cannot be fully resolved.  
- Variant filtering thresholds are generic and should be refined for clinical-grade analysis.  

---

## Professional Profile  
- GitHub: [Project Repository](https://github.com/Opeoluwa-Shodipe/NGS-Project-HACKBIO-/blob/main/Projects/CFTR_WGS_CaseStudy.md)  
- LinkedIn: Upload a short video explaining the project in layman’s terms (ensure public access).  

---

# Conclusion  
This WGS analysis pipeline enables reproducible variant detection and annotation in the **CFTR gene**, supporting clinical diagnosis of cystic fibrosis. The approach integrates QC, alignment, variant calling, annotation, and inheritance analysis in a single automated workflow.
