#!/bin/bash
# run_cftr_pipeline.sh
set -euo pipefail
THREADS=8
OUTDIR=results
mkdir -p $OUTDIR

# --- User variables (edit) ---
REF=/path/to/GRCh38.primary_assembly.genome.fa
CHILD_1=child_1.fastq.gz
CHILD_2=child_2.fastq.gz
FATHER_1=father_1.fastq.gz
FATHER_2=father_2.fastq.gz
CFTR_REGION="7:117120016-117308718"    # GRCh38 coordinates for CFTR

# --- 0. Check tools ---
for t in fastp bwa samtools picard bcftools bgzip tabix vep; do
  command -v $t >/dev/null 2>&1 || { echo "$t is required"; exit 1; }
done

# --- 1. QC & trimming (fastp) ---
fastp -i $CHILD_1 -I $CHILD_2 \
      -o $OUTDIR/child_1.trim.fastq.gz -O $OUTDIR/child_2.trim.fastq.gz \
      -w $THREADS -h $OUTDIR/child_fastp.html -j $OUTDIR/child_fastp.json

fastp -i $FATHER_R1 -I $FATHER_R2 \
      -o $OUTDIR/father_1.trim.fastq.gz -O $OUTDIR/father_2.trim.fastq.gz \
      -w $THREADS -h $OUTDIR/father_fastp.html -j $OUTDIR/father_fastp.json

# --- 1b. Optional: fastqc reporting ---
# fastqc -t $THREADS $OUTDIR/*.trim.fastq.gz -o $OUTDIR/fastqc

# --- 2. Reference indexing (if not already done) ---
[ -f ${REF}.bwt ] || bwa index $REF
[ -f ${REF}.fai ] || samtools faidx $REF
[ -f ${REF%.fa}.dict ] || java -jar $(which picard) CreateSequenceDictionary R=$REF O=${REF%.fa}.dict

# --- 3. Alignment (bwa mem) + convert to BAM ---
bwa mem -t $THREADS $REF $OUTDIR/child_1.trim.fastq.gz $OUTDIR/child_2.trim.fastq.gz \
  | samtools view -b -@ $THREADS -o $OUTDIR/child.unsorted.bam -

bwa mem -t $THREADS $REF $OUTDIR/father_1.trim.fastq.gz $OUTDIR/father_2.trim.fastq.gz \
  | samtools view -b -@ $THREADS -o $OUTDIR/father.unsorted.bam -

# --- 4. Sort, mark duplicates, index ---
samtools sort -@ $THREADS -o $OUTDIR/child.sorted.bam $OUTDIR/child.unsorted.bam
samtools sort -@ $THREADS -o $OUTDIR/father.sorted.bam $OUTDIR/father.unsorted.bam

picard MarkDuplicates I=$OUTDIR/child.sorted.bam O=$OUTDIR/child.dedup.bam M=$OUTDIR/child.dup.metrics.txt
picard MarkDuplicates I=$OUTDIR/father.sorted.bam O=$OUTDIR/father.dedup.bam M=$OUTDIR/father.dup.metrics.txt

samtools index $OUTDIR/child.dedup.bam
samtools index $OUTDIR/father.dedup.bam

# --- 5. Coverage check across CFTR (quick) ---
samtools depth -r $CFTR_REGION $OUTDIR/child.dedup.bam | awk '{sum+=$3; cnt++} END{print "child_mean_dp", sum/cnt}'
samtools depth -r $CFTR_REGION $OUTDIR/father.dedup.bam | awk '{sum+=$3; cnt++} END{print "father_mean_dp", sum/cnt}'

# --- 6. Joint variant calling restricted to CFTR (bcftools) ---
bcftools mpileup -f $REF -r $CFTR_REGION $OUTDIR/child.dedup.bam $OUTDIR/father.dedup.bam -Ou \
  | bcftools call -mv -Oz -o $OUTDIR/joint.cftr.vcf.gz
bcftools index $OUTDIR/joint.cftr.vcf.gz

# --- 7. Filter variants (example) ---
# Keep variants with QUAL>=20 and DP>=10 (tweak thresholds to data)
bcftools view -i 'QUAL>=20 && INFO/DP>=10' $OUTDIR/joint.cftr.vcf.gz -Oz -o $OUTDIR/joint.cftr.filtered.vcf.gz
bcftools index $OUTDIR/joint.cftr.filtered.vcf.gz

# --- 8. Annotate (VEP) ---
vep -i $OUTDIR/joint.cftr.filtered.vcf.gz \
    --cache --assembly GRCh38 --offline \
    --fasta $REF \
    --everything \
    --vcf -o $OUTDIR/joint.cftr.vep.vcf.gz --compress_output bgzip --thread $THREADS
tabix -p vcf $OUTDIR/joint.cftr.vep.vcf.gz

# --- 9. Create summary tables ---
# Basic annotation table (includes CLIN_SIG if present)
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ\t%INFO/CLIN_SIG\t%INFO/AF\t[%SAMPLE=%GT:%DP]\n' \
  $OUTDIR/joint.cftr.vep.vcf.gz > $OUTDIR/joint_cftr_annotation.tsv

# Compare genotypes per sample
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE=%GT\t]\n' $OUTDIR/joint.cftr.vep.vcf.gz > $OUTDIR/compare_genotypes.tsv

echo "Done. Check $OUTDIR for outputs."


# -------- Professional Profile --------
echo "Professional Profile"
echo "1. GitHub: https://github.com/Opeoluwa-Shodipe/NGS-Project-HACKBIO-/blob/main/Projects/Stage%20one.md"
echo "2. LinkedIn: Ensure video is uploaded and public."





