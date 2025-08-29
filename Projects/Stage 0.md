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
