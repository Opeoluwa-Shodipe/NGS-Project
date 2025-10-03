# Transcriptomic Profiling of *Staphylococcus aureus* During Acute vs Chronic Phases of Periprosthetic Joint Infection (PJI)

## Background and Rationale

Periprosthetic joint infections (PJIs) are among the most devastating complications of orthopedic implants. They increase morbidity, prolong hospital stays, and often require costly revision surgeries. *Staphylococcus aureus*—particularly methicillin-resistant strains (MRSA)—is a leading cause of PJIs.  

One critical feature of *S. aureus* is its ability to switch phenotypes between acute and chronic infection phases:

- **Acute phase:** Bacteria adopt an aggressive, planktonic growth mode, expressing virulence factors such as toxins, adhesins, and immune evasion genes.  
- **Chronic phase:** Bacteria adapt to a biofilm-like state, downregulating overt virulence and upregulating persistence pathways (stress response, metabolic rewiring, antibiotic tolerance).  

This adaptive flexibility makes chronic PJIs notoriously difficult to eradicate. Antibiotic regimens often fail, and host immune responses are blunted by biofilm shielding.  

---

## Why RNA-seq?

RNA sequencing provides a window into the global transcriptional programs that underpin this acute-to-chronic transition. By capturing gene expression profiles directly from *S. aureus* isolates in different clinical phases of PJI, RNA-seq can:

- Identify virulence genes uniquely expressed in acute infection.  
- Detect metabolic and stress-response pathways that dominate during chronic infection.  
- Reveal regulatory RNAs and transcriptional signatures linked to biofilm persistence.  

Such insights are crucial for designing diagnostics that distinguish acute from chronic PJIs, and for developing therapies that specifically disrupt persistence mechanisms.  

---

## Project Goals

### Preprocessing and Quality Control
- Perform read trimming, alignment to the *S. aureus* reference genome, and assessment of sequencing quality.  
- Generate count matrices of gene expression for acute and chronic isolates.  

### Differential Gene Expression Analysis
- Compare transcriptional profiles of *S. aureus* from acute vs chronic PJIs.  
- Identify virulence, stress response, and metabolic genes that are significantly up- or downregulated.  

### Functional Enrichment and Pathway Mapping
- Conduct GO/KEGG pathway enrichment analyses.  
- Highlight pathways associated with biofilm formation, immune evasion, and antibiotic resistance.  

### Regulatory and Adaptation Insights
- Explore small RNAs and regulatory elements potentially shaping the acute/chronic shift.  

### Visualization and Communication
- Generate PCA plots to separate acute vs chronic isolates.  
- Create volcano plots and heatmaps of key DEGs.  
- Write a “clinical microbiology report” summarizing the molecular strategies *S. aureus* deploys in acute vs chronic infection.  

---

## Dataset

Use **SRA-Explorer** to download the 8 samples from **PRJNA867318**.  

| S/N | Accession Number | State                                   |
|-----|------------------|-----------------------------------------|
| 1   | SRR20959676      | chronic periprosthetic joint infection  |
| 2   | SRR20959677      | chronic periprosthetic joint infection  |
| 3   | SRR20959678      | chronic periprosthetic joint infection  |
| 4   | SRR20959679      | chronic periprosthetic joint infection  |
| 5   | SRR20959680      | acute periprosthetic joint infection    |
| 6   | SRR20959681      | acute periprosthetic joint infection    |
| 7   | SRR20959682      | acute periprosthetic joint infection    |
| 8   | SRR20959683      | acute periprosthetic joint infection    |

---
