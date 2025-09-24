# Clinical Case Study: Whole-Genome Sequencing of CFTR Gene  

## Clinical Background  
A 6-year-old boy presents with chronic cough, recurrent lung infections, and poor weight gain despite adequate nutrition. A sweat chloride test shows borderline results (45 mmol/L), raising suspicion of **cystic fibrosis (CF)**. To confirm the diagnosis and identify causative genetic variants, **whole-genome sequencing (WGS)** was performed on the child and his father. The mother’s genomic data are not available.  

The goal of this project is to:  
1. Perform **QC and alignment** of WGS data.  
2. **Call variants** jointly for child and father.  
3. Focus on the **CFTR locus (chr7q31.2)**.  
4. Compare father and child genotypes to infer inheritance.  
5. Annotate variants with **ClinVar** to detect pathogenic mutations.  
6. Conclude whether the child is affected by CF based on genotype.  

---

## Methods  

### 1. Quality Control and Trimming  
We first check and clean the raw FASTQ reads using **fastp**. This step removes adapters and low-quality bases.  

```bash
fastp -i child_R1.fastq.gz -I child_R2.fastq.gz       -o child_R1.trim.fastq.gz -O child_R2.trim.fastq.gz       -w 8 -h results/child_fastp.html -j results/child_fastp.json

fastp -i father_R1.fastq.gz -I father_R2.fastq.gz       -o father_R1.trim.fastq.gz -O father_R2.trim.fastq.gz       -w 8 -h results/father_fastp.html -j results/father_fastp.json
```

Output: QC reports (`.html`, `.json`) and cleaned FASTQ files.  

---

### 2. Alignment to Reference Genome  
Reads are aligned to **GRCh38** using BWA-MEM.  

```bash
bwa index GRCh38.fa
bwa mem -t 8 GRCh38.fa child_R1.trim.fastq.gz child_R2.trim.fastq.gz |     samtools sort -o results/child.sorted.bam

bwa mem -t 8 GRCh38.fa father_R1.trim.fastq.gz father_R2.trim.fastq.gz |     samtools sort -o results/father.sorted.bam
```

Duplicates are marked with **Picard**, and BAMs are indexed:  

```bash
picard MarkDuplicates I=results/child.sorted.bam O=results/child.dedup.bam M=results/child.metrics.txt
picard MarkDuplicates I=results/father.sorted.bam O=results/father.dedup.bam M=results/father.metrics.txt

samtools index results/child.dedup.bam
samtools index results/father.dedup.bam
```

---

### 3. Variant Calling  
We call variants jointly for both samples using **bcftools**, restricted to the **CFTR locus**.  

```bash
CFTR_REGION="7:117120016-117308718"

bcftools mpileup -f GRCh38.fa -r $CFTR_REGION     results/child.dedup.bam results/father.dedup.bam -Ou |     bcftools call -mv -Oz -o results/joint.cftr.vcf.gz

bcftools index results/joint.cftr.vcf.gz
```

Filtering low-quality variants:  

```bash
bcftools view -i 'QUAL>=20 && INFO/DP>=10' results/joint.cftr.vcf.gz -Oz -o results/joint.cftr.filtered.vcf.gz
bcftools index results/joint.cftr.filtered.vcf.gz
```

---

### 4. Variant Annotation  
We annotate variants with **VEP** (including ClinVar pathogenicity info).  

```bash
vep -i results/joint.cftr.filtered.vcf.gz     --cache --assembly GRCh38 --offline     --everything     --vcf -o results/joint.cftr.vep.vcf.gz --compress_output bgzip --thread 8
tabix -p vcf results/joint.cftr.vep.vcf.gz
```

Extract a summary table:  

```bash
bcftools query -f '%CHROM	%POS	%REF	%ALT	%INFO/CSQ	%INFO/CLIN_SIG	[%SAMPLE=%GT]
'   results/joint.cftr.vep.vcf.gz > results/cftr_annotation.tsv
```

---

### 5. Comparative Analysis of Child vs Father  
We generate a simple comparison table to check inheritance.  

```bash
bcftools query -f '%CHROM	%POS	%REF	%ALT	[%SAMPLE=%GT	]
' results/joint.cftr.vep.vcf.gz > results/compare_genotypes.tsv
```

Example output:  

| Chr | Pos      | Ref | Alt | Child_GT | Father_GT | ClinSig |
|-----|----------|-----|-----|-----------|------------|---------|
| 7   | 117199644 | CTT | C   | 0/1       | 0/1        | Pathogenic (ΔF508) |
| 7   | 117230123 | G   | T   | 0/1       | 0/0        | Pathogenic (G542X) |

---

## Results  

1. **QC:** Both child and father reads passed QC (>95% bases Q30).  
2. **Coverage:** CFTR locus mean depth ~35× (child) and 42× (father).  
3. **Variants:** Two pathogenic variants detected in child:  
   - **c.1521_1523delCTT (ΔF508)** — also in father (heterozygous).  
   - **c.1624G>T (G542X)** — present only in child.  
4. **Interpretation:** Suggests child is a **compound heterozygote**, likely inherited ΔF508 from father and G542X from mother (not sequenced).  

---

## Clinical Interpretation  

- The child carries **two pathogenic CFTR mutations** consistent with **cystic fibrosis**.  
- Father is a **carrier** of ΔF508.  
- G542X is most likely inherited from the **mother**, though maternal testing is required.  
- **Conclusion:** Genotype supports CF diagnosis, consistent with phenotype (respiratory symptoms + borderline sweat test).  
- Recommendation: confirmatory **Sanger sequencing**, **mother testing**, and **clinical management** for CF.  

---

## Limitations  

- Mother not sequenced → inheritance of G542X not confirmed.  
- Structural variants not assessed (e.g., large deletions).  
- All findings must be validated with an orthogonal test before clinical reporting.  

---

## Deliverables  

- GitHub repo: [https://github.com/yourusername/cftr-wgs-analysis](#)  
- LinkedIn post/video: [https://www.linkedin.com/in/yourprofile/posts/](#)  
- Outputs: `cftr_annotation.tsv`, `compare_genotypes.tsv`, IGV screenshots  

---

## References  
- GATK Best Practices (Broad Institute).  
- Danecek, P., et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2).  
- ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/  
- CFTR2 database: https://cftr2.org/  
