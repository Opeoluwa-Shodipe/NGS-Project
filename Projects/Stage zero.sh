#!/bin/bash
# ======================================================
# Project Script: Bash Basics + Bioinformatics Setup
# Author: Opeoluwa
# ======================================================

# ---------------------------
# Project 1: Bash Basics
# ---------------------------

# 1. Print your name
echo "Opeoluwa"

# 2. Create a folder titled your name
mkdir -p Opeoluwa

# 3. Create 'biocomputing' directory and change into it
mkdir -p biocomputing && cd biocomputing || exit

# 4. Download the 2 required files (optimized with wget)
wget -O wildtype.fna https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna
wget -O wildtype.gbk https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk

# 5. Move the .fna file to the folder titled your name
mv wildtype.fna ../Opeoluwa/

# 6. Delete duplicate .gbk file if it exists
[ -f wildtype.gbk.1 ] && rm -r wildtype.gbk.1

# 7. Confirm if .fna file is mutant or wild type
cd ../Opeoluwa || exit
if grep -q "tatatata" wildtype.fna; then
    echo "Mutant"
    # 8. If mutant, print matching lines into a new file
    grep "tatatata" wildtype.fna > mutant_lines.txt
else
    echo "Wild type"
fi

# 9. Count number of lines (excluding header) in .gbk file
cd ../biocomputing || exit
echo "Number of lines (excluding header):"
tail -n +2 wildtype.gbk | wc -l

# 10. Print sequence length from LOCUS tag
echo "Sequence length from LOCUS:"
grep "LOCUS" wildtype.gbk | awk '{print $2}'

# 11. Print source organism from SOURCE tag
echo "Source organism:"
grep "SOURCE" wildtype.gbk | head -n 1 | awk '{$1=""; print $0}'

# 12. List all gene names
echo "Gene names:"
grep "/gene=" wildtype.gbk

# 13. Print all commands used today (instead of clearing terminal)
echo "Command history (today):"
history

# 14. List files in both folders
echo "Files in Opeoluwa folder:"
ls ../Opeoluwa
echo "Files in biocomputing folder:"
ls


# ---------------------------
# Project 2: Bioinformatics Software Installation
# ---------------------------

# 1. Download the installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# 2. Run the installer
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda

# 3. Activate your base conda environment
export PATH="$HOME/miniconda/bin:$PATH"
conda --version

# 4. Create a conda environment named funtools
conda create -n funtools -y

# 5. Activate the funtools environment
source activate funtools

# 6. Install Figlet using conda
conda install -c conda-forge -y figlet

# 7. Run figlet with your name
figlet Opeoluwa

# 8. Install bioinformatics tools via loop
tools=(bwa blast samtools bedtools spades bcftools fastp multiqc)
for tool in "${tools[@]}"; do
    echo "Installing $tool..."
    conda install -c bioconda -y "$tool"
done

echo "All tools installed successfully."


# ---------------------------
# Professional Profile
# ---------------------------
echo "Professional Profile:"
echo "1. GitHub: https://github.com/Opeoluwa-Shodipe/NGS-Project-HACKBIO-/blob/main/Projects/Stage%200.md"
echo "2. LinkedIn: (Add your LinkedIn video link here)"
