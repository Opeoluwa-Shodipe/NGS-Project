# Project 1: Bash Basic

```bash
#!/bin/bash

# 1. Print your name
echo "Opeoluwa"

# 2. Create a folder titled your name
mkdir Opeoluwa

# 3. Create 'biocomputing' directory and change into it in one line
mkdir biocomputing && cd biocomputing

# 4. Download the 3 files (optimized with wget)
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna 
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk 
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk 

# 5. Move the .fna file to the folder titled your name
mv wildtype.fna ../Opeoluwa/

# 6. Delete the duplicate .gbk file (keep one)
rm -r wildtype.gbk.1  

# 7. Confirm if .fna file is mutant or wild type
cd ../Opeoluwa
if grep -q "tatatata" wildtype.fna; then
    echo "Mutant"
else
    echo "Wild type"
fi

# 8. If mutant, print matching lines into a new file
grep "tatatata" wildtype.fna > mutant_lines.txt

# 9. Count number of lines (excluding header) in .gbk file
cd ../biocomputing
tail -n +2 wildtype.gbk | wc -l

# 10. Print sequence length from LOCUS tag
grep "LOCUS" wildtype.gbk | awk '{print $2}'

# 11. Print source organism from SOURCE tag
grep "SOURCE" wildtype.gbk | head -n 1 | awk '{$1=""; print $0}'

# 12. List all gene names
grep "/gene=" wildtype.gbk

# 13. Clear terminal and print all commands used today
clear
history

# 14. List files in both folders
echo "Files in Opeoluwa folder:"
ls ../Opeoluwa
echo "Files in biocomputing folder:"
ls


# Project 2: Installing Bioinformatics Software on the Terminal


```bash
# 1. Download the installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh 

# 2. Run the installer
bash Miniconda3-latest-Linux-x86_64.sh 

# 3. Activate your base conda environment (after a quick restart)
conda --version

# 4. Create a conda environment named funtools
conda create -n funtools

# 5. Activate the funtools environment
conda activate funtools

# 6. Install Figlet using conda (or apt-get alternative)
conda install -c conda-forge figlet

# 7. Run figlet with your name
figlet Opeoluwa

# 8. Install bwa through the bioconda channel
conda install -c bioconda bwa

# 9. Install blast through the bioconda channel
conda install -c bioconda blast

# 10. Install samtools through the bioconda channel
conda install -c bioconda samtools

# 11. Install bedtools through the bioconda channel
conda install -c bioconda bedtools

# 12. Install spades.py through the bioconda channel
conda install -c bioconda spades

# 13. Install bcftools through the bioconda channel
conda install -c bioconda bcftools

# 14. Install fastp through the bioconda channel
conda install -c bioconda fastp

# 15. Install multiqc through the bioconda channel
conda install -c bioconda multiqc


