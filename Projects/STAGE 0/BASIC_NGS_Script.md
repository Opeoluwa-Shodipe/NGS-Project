#!/bin/bash
# =====================================================
# Project Script: Bash Basics + Bioinformatics Software
# Author: Opeoluwa
# =====================================================

# -------- Pre-checks --------
echo "Checking dependencies..."
command -v wget >/dev/null 2>&1 || { echo "wget not found. Please install it."; exit 1; }
command -v conda >/dev/null 2>&1 || { echo "conda not found. Please install Miniconda/Anaconda."; exit 1; }

# -------- Project 1: Bash Basics --------
echo "Starting Project 1: Bash Basics"

# 1. Print your name
echo "Name: Opeoluwa"

# 2. Create a folder titled your name (check if it exists first)
if [ ! -d "../Opeoluwa" ]; then
    mkdir -p ../Opeoluwa
    echo "Created directory: Opeoluwa"
else
    echo "Directory 'Opeoluwa' already exists."
fi

# 3. Create 'biocomputing' directory and change into it
if [ ! -d "../biocomputing" ]; then
    mkdir -p ../biocomputing
    echo "Created directory: biocomputing"
fi
cd ../biocomputing || exit 1

# 4. Download the 2 required files separately
for file in wildtype.fna wildtype.gbk; do
    if wget -q https://raw.githubusercontent.com/josoga2/dataset-repos/main/$file; then
        echo "Downloaded $file"
    else
        echo "Failed to download $file"
        exit 1
    fi
done

# 5. Move the .fna file to the folder titled your name
mv -f wildtype.fna ../Opeoluwa/
echo "Moved wildtype.fna to Opeoluwa folder."

# 6. Handle duplicate .gbk file if present
duplicates=$(ls wildtype.gbk.* 2>/dev/null | wc -l)
if [ "$duplicates" -gt 0 ]; then
    echo "Found duplicate .gbk files. Removing..."
    rm -f wildtype.gbk.*
else
    echo "No duplicate .gbk files found."
fi

# 7. Confirm if .fna file is mutant or wild type
cd ../Opeoluwa || exit 1
if grep -q "tatatata" wildtype.fna; then
    echo "Mutant sequence detected"
    grep -n "tatatata" wildtype.fna
elif grep -q "tata" wildtype.fna; then
    echo "Wild type sequence detected"
else
    echo "No known sequence pattern found"
fi

# 8. Save mutant lines to a new file (if any)
grep "tatatata" wildtype.fna > mutant_lines.txt
echo "Mutant lines saved to mutant_lines.txt"

# 9. Count number of lines (excluding LOCUS header) in .gbk file
cd ../biocomputing || exit 1
line_count=$(grep -v "^LOCUS" wildtype.gbk | wc -l)
echo "Number of lines in wildtype.gbk (excluding LOCUS): $line_count"

# 10. Print sequence length from LOCUS tag
seq_length=$(grep "LOCUS" wildtype.gbk | awk '{print $2}')
echo "Sequence length: $seq_length"

# 11. Print source organism from SOURCE tag
organism=$(grep "SOURCE" wildtype.gbk | head -n 1 | awk '{$1=""; print $0}')
echo "Source organism: $organism"

# 12. List all gene names
echo "Gene names found:"
grep "/gene=" wildtype.gbk

# 13. List files in both folders
echo "Files in Opeoluwa folder:"
ls ../Opeoluwa
echo "Files in biocomputing folder:"
ls


# -------- Project 2: Installing Bioinformatics Software --------
echo "Starting Project 2: Installing Bioinformatics Software"

# 1. Download the Miniconda installer
wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
echo "Miniconda installer downloaded"

# 2. Run the installer (non-interactive mode)
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda

# 3. Confirm conda installation
source "$HOME/miniconda/etc/profile.d/conda.sh"
conda --version

# 4. Create a conda environment named funtools
conda create -n funtools -y
conda activate funtools

# 5. Install Figlet using conda
conda install -c conda-forge figlet -y
figlet Opeoluwa

# 6. Install bioinformatics tools in a loop
tools=(bwa blast samtools bedtools spades bcftools fastp multiqc)
for tool in "${tools[@]}"; do
    echo "Installing $tool..."
    conda install -c bioconda -y "$tool"
done

echo "All bioinformatics tools installed successfully."


# -------- Professional Profile --------
echo "Professional Profile"
echo "1. GitHub: https://github.com/Opeoluwa-Shodipe/NGS-Project-HACKBIO-/blob/main/Projects/Stage%20zero.sh"
echo "2. LinkedIn: Ensure video is uploaded and public."

