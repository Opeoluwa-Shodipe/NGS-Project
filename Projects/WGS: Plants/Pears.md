# Pear Genomics: Whole Genome Sequencing Case Study  

In the lush orchards of Asia and Europe, pear trees have thrived for centuries, their fruits cherished for their sweetness, texture, and diversity. Yet, beneath their juicy exteriors lies a complex story of domestication, adaptation, and genetic variation.  

The Food and Agriculture Organization (FAO) has been interested in understanding how pears can contribute to global food security, most especially through improving the nutritional quality of domesticated fruits. You are presented with WGS data from a study revealing telomere-to-telomere (T2T) haplotype-resolved genomes of two pear cultivars, **‘Dangshansuli’ (Asian)** and **‘Max Red Bartlett’ (European)**.  

Your job is to uncover structural variations (SVs), Single Nucleotide Polymorphisms (SNPs), and domestication signatures that shape fruit traits like **size, sweetness, and stone cell content**. These insights open a window into pear evolution, but many questions remain:  

- How do SVs differ across pear populations?  
- Can we pinpoint genetic markers for superior nutritional quality?  

Your mission as bioinformatics students is to tackle this agricultural puzzle, using WGS to explore pear genomic diversity and inform breeding for better, tastier pears—a challenge with global implications for food security and sustainability.  

---

## Goals  

1. Identify the different species of pears collected in the project  
2. Perform high quality variant detection in the samples  
3. Identify SNP and SVs that are likely to shape nutritional quality and taste of pears between Europe and Asia  

---

## Dataset  

Use **SRA-Explorer** to download the datasets from the following project:  

- **Asian:** `SRR32764071`, `SRR32764072`  
- **European:** `SRR7135493`, `SRR7135521`  

---

## Expected Deliverables  

We expect a **code report submitted as a markdown**, i.e. a mix of code and text explanation.  
See example: [Scripting Example](https://github.com/josoga2/bash-course/blob/main/bash/module_7/scripting_1.md)  

---

## Grading Rubric (Total = 10 points)  

### **1. Data Acquisition & QC (1.0 pt)**  
- Raw data correctly downloaded from SRR IDs  
- QC with **FastQC/MultiQC** and trimming documented  
- Coverage depth reported  

### **2. Species Identification (2 pts)**  
- Correct assignment of Asian vs European samples  

### **3. SNP/Small Variant Calling (3 pts)**  
- High-quality SNP/indel detection with modern caller (e.g., GATK)  
- Filtering justified  
- Variant QC metrics shown  

### **4. Reporting (3 pts)**  
- Clear and concise interpretation in a written report  
- If possible, link to GitHub for images from the genome viewer  

### **5. Critical Thinking & Limitations (1 pt)**  
- Acknowledges small sample size, reference bias, and potential false positives  
- Suggests next experiments or validation steps  

---

## Outreach Component  

Provide the link to a **LinkedIn Post** containing either:  

- A **video recording** discussing the project (problems, results, and importance) in a way a layman will understand  
**OR**  
- A **long-form post** discussing the project (problems, results, and importance) in a way a layman will understand  

---
