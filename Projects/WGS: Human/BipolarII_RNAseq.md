# Unmasking Molecular Signatures of Bipolar II Disorder

Bipolar disorder (BD) is a chronic psychiatric illness characterized by mood instability, high relapse rates, and significant social and economic burden. While genetic heritability is strong, the underlying molecular pathways remain poorly defined, especially in Bipolar II disorder (BD-II), where depressive episodes predominate and hypomania is often underdiagnosed.  

Microglia—the brain’s resident immune cells—are increasingly implicated in neuropsychiatric disorders. They regulate synaptic pruning, inflammatory signaling, and neurodevelopmental processes. Patient-derived induced pluripotent stem cells (iPSCs) differentiated into microglia-like cells provide a human-relevant model to study cell-intrinsic molecular differences in BD-II without relying solely on post-mortem tissue.  

Both familial BD-II (FBD) and sporadic BD-II (SBD) are clinically recognized, yet whether they share molecular signatures or diverge in their immune-cell transcriptomes is largely unexplored. Establishing transcriptional differences between these subgroups, compared to matched controls (FHC and SHC), could illuminate disease-specific pathways and point toward biomarkers or therapeutic targets.  

---

## Why RNA-seq?

RNA sequencing captures the global transcriptional landscape, enabling unbiased discovery of:

- **Differentially expressed genes (DEGs)** between disease and control groups.  
- **Pathways** enriched in immune signaling, synaptic regulation, or metabolic function.  
- **Subtle distinctions** between familial and sporadic forms of BD-II.  

Given the hypothesis that BD involves dysregulated immune and synaptic biology, RNA-seq in microglia-like cells is the logical first step to map the transcriptomic perturbations at genome-wide resolution.  

---

## Project Goals

### 1. Quality Control and Preprocessing  
- Perform adapter trimming, quality filtering, and alignment to the human reference genome.  
- Evaluate read quality (FastQC, MultiQC) and alignment rates.  

### 2. Quantification of Gene Expression  
- Generate count matrices for all samples (SBD, FBD, SHC, FHC).  
- Normalize expression values to account for sequencing depth and library composition.  

### 3. Differential Expression Analysis  
Compare gene expression between:  
- SBD vs SHC  
- FBD vs FHC  
- SBD vs FBD (to detect familial/sporadic distinctions; e.g., PCA to distinguish between the two)  

Identify significantly up- or downregulated genes in each contrast.  

### 4. Pathway and Functional Enrichment  
- Conduct GO and KEGG enrichment analyses.  
- Highlight immune-related, neurodevelopmental, or synaptic signaling pathways dysregulated in BD-II.  

### 5. Integrative Analysis  
- Compare overlap in DEGs between SBD and FBD.  
- Explore unique vs shared transcriptional programs, focusing on immune activation vs synaptic regulation.  

### 6. Reporting  
- Deliver visualizations: volcano plots, heatmaps of DEGs, pathway enrichment diagrams.  
- Write a concise interpretative report to psychiatric genomics stakeholders, outlining:  
  - Common vs distinct molecular alterations in SBD and FBD.  
  - Potential biomarkers and therapeutic pathway candidates.  

---

## Dataset

| S/N | SRA Accession Number | Age | State             |
|-----|-----------------------|-----|------------------|
| 1   | SRR33243164           | 28  | Bipolar disorder |
| 2   | SRR33243165           | 33  | Bipolar disorder |
| 3   | SRR33243166           | 38  | Bipolar disorder |
| 4   | SRR33243167           | 36  | Bipolar disorder |
| 5   | SRR33243168           | 25  | Healthy control  |
| 6   | SRR33243169           | 24  | Healthy control  |
| 7   | SRR33243170           | 45  | Healthy control  |
| 8   | SRR33243171           | 67  | Healthy control  |

---
