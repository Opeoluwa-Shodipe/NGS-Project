#download genome
wget http://ftp.ensemblgenomes.org/pub/metazoa/release-61/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz
gunzip Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz 

# rename the file to c_elegans.fa
mv Caenorhabditis_elegans.WBcel235.dna.toplevel.fa c_elegans.fa
# create genome index directory
STAR --runMode genomeGenerate --genomeDir genomeIndex --genomeFastaFiles c_elegans.fa

# create directories for output
mkdir -p mapped
#! /bin/bash
for infile in trim/*.fq ; do
	outfile=$(basename "$infile".fq)
	STAR --genomeDir genome/genomeIndex --readFilesIn $infile --outFileNamePrefix mapped/$outfile --outSAMtype BAM SortedByCoordinate --outSAMattributes All
done

#play with the bam output
samtools view mapped/fm_1.fq.fq_trim.fq.fqAligned.sortedByCoord.out.bam | head
samtools view -c mapped/fm_1.fq.fq_trim.fq.fqAligned.sortedByCoord.out.bam | head
samtools  flagstat mapped/fm_1.fq.fq_trim.fq.fqAligned.sortedByCoord.out.bam | head

# create a directory for IGV viewing
mkdir IGV
cp mapped/*.bam IGV/
