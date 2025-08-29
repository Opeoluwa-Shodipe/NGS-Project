#!/bin/bash
# =====================================================
# Project Script: Bash Basics + Bioinformatics Software
# Author: Opeoluwa
# =====================================================

# -------- Project 1: Bash Basics --------

# 1. Print your name
echo "Opeoluwa"

# 2. Create a folder titled your name
mkdir -p Opeoluwa

# 3. Create 'biocomputing' directory and change into it in one line
mkdir -p biocomputing && cd biocomputing || exit 1

# 4. Download the 2 required files (optimized wget)
wget -O wildtype.fna https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna \
     -O wildtype.gbk https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk

# 5. Move the .fna file to the folder titled your name
mv wildtype.fna ../Opeoluwa/

# 6. Remove duplicate .gbk file if it exists
rm -r wildtype.gbk.1  

# 7. Confirm if .fna file is mutant or wild type
cd ../Opeoluwa || exit 1
if grep -q "tatatata" wildtype.fna; then
    echo "Mutant sequence detected"
elif grep -q "tata" wildtype.fna; then
    echo "Wild type sequence detected"
else
    echo "No known sequence pattern found"
fi

# 8. If mutant, print matching lines into a new file
grep "tatatata" wildtype.fna > mutant_lines.txt

# 9. Count number of lines (excluding header) in .gbk file
cd ../biocomputing || exit 1
grep -v "^LOCUS" wildtype.gbk | wc -l

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


# -------- Project 2: Installing Bioinformatics Software --------

# 1. Download the Miniconda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh 

# 2. Run the installer
bash Miniconda3-latest-Linux-x86_64.sh 

# 3. Confirm conda installation
conda --version

# 4. Create a conda environment named funtools
conda create -n funtools -y

# 5. Activate the funtools environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate funtools

# 6. Install Figlet using conda (or apt-get alternative)
conda install -c conda-forge figlet -y

# 7. Run figlet with your name
figlet Opeoluwa

# 8. Install bioinformatics tools in a loop
for tool in bwa blast samtools bedtools spades bcftools fastp multiqc; do
    echo "Installing $tool..."
    conda install -c bioconda -y "$tool"
done

echo "All bioinformatics tools installed successfully!"

# -------- Professional Profile --------
echo "Professional Profile"
echo "1. GitHub: https://github.com/Opeoluwa-Shodipe/NGS-Project-HACKBIO-/blob/main/Projects/Stage%20zero.sh"
echo "2. LinkedIn: Video link (to be added)"
