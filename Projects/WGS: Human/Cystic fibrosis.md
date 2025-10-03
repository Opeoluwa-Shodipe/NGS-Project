# Clinical Case Study: Whole-Genome Sequencing

## Clinical Background
A 6-year-old boy presents with chronic cough, recurrent lung infections, and poor weight gain despite adequate nutrition. His physician suspects cystic fibrosis (CF). A sweat chloride test has been ordered but yields borderline results (45 mmol/L). To confirm diagnosis and identify the causative genetic variant(s), the clinical team orders whole-genome sequencing (WGS).

The boy’s father has also been sequenced as part of the family study. The mother’s genomic data are not available.

---

## Your Task
You are part of the bioinformatics team tasked with analyzing the WGS data to determine whether the child carries pathogenic variants in the CFTR gene.

Specifically, you should:

- Perform quality control and alignment of raw WGS reads to the human reference genome.  
- Call variants in both father and child.  
- Focus analysis on the CFTR locus (chromosome 7q31.2).  
- Compare the child’s variants against the father’s to determine inheritance patterns.  
- Annotate and interpret variants to identify known pathogenic CFTR mutations (e.g., ΔF508, G542X, W1282X, etc.) using variant databases such as ClinVar.  
- Conclude whether the child is affected by CF based on genotype, and comment on whether the observed mutations are inherited or potentially de novo.  

---

## Deliverables
- A concise report (2–3 pages) including methods, results, and clinical interpretation.  
- Visualizations or tables of identified variants (VCF excerpts, IGV screenshots, or annotation summaries).  
- Final conclusion on the likely genetic basis of the child’s phenotype.  

---

## Grading Rubric (Total = 10 points)
We expect a code report submitted as a markdown. i.e. a mix of code and text explanation.  
See example: [Bash Course Example](https://github.com/josoga2/bash-course/blob/main/bash/module_7/scripting_1.md)

| Criterion                        | Description                                                                                   | Points |
|----------------------------------|-----------------------------------------------------------------------------------------------|--------|
| Data processing & QC             | Correct use of quality control (FastQC, trimming) and alignment (BWA/other) with justification | 2      |
| Variant calling & filtering      | Accurate execution of variant calling, filtering                                               | 2      |
| Comparative analysis             | Clear comparison of father vs child variants, highlighting inheritance                        | 2      |
| Variant annotation & interpretation | Correct use of annotation tools/databases; identification of pathogenic variants             | 2      |
| Report quality & conclusions     | Clear, logical report with accurate clinical interpretation and discussion of limitations      | 2      |

---

## Submission Requirements
- Provide the link to a **GitHub repo** containing your script.  
- Provide the link to a **LinkedIn Post** containing either:  
  - A video recording discussing the project (problems, results, and importance) in a way a layman will understand.  
  - Or a long-form post discussing the project (problems, results, and importance) in a way a layman will understand.  
