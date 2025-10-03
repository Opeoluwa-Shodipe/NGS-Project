# South African Polony Outbreak – Whole-Genome Sequencing (WGS) Project

## Background  

In early 2017, South Africa faced a **public health crisis** that would become the world's largest recorded outbreak for a particular infection. Doctors at Chris Hani and Steve Biko academic hospitals noticed an alarming spike in neonatal infections, signaling that something sinister was spreading through the food supply.  

By March 2018, the **National Institute for Communicable Diseases (NICD)** had reported:  
- **978 laboratory-confirmed cases**  
- **183 deaths (27% case fatality rate)**  

This was devastatingly high for a foodborne illness. Most cases were concentrated in:  
- Gauteng (**59%**)  
- Western Cape (**12%**)  
- KwaZulu-Natal (**7%**)  

Vulnerable groups—neonates, pregnant women, the elderly, and immunocompromised patients (especially those living with HIV)—were hardest hit.  

### Tracing the Source  

Interviews with patients pointed to processed cold meats, particularly **polony**, as the likely culprit. Further investigation traced many cases to the **Enterprise Foods facility in Polokwane**.  

This outbreak was not limited to South Africa—polony products had been exported to **15 African countries**, raising fears of a regional spread.  

### Why Whole-Genome Sequencing (WGS)?  

To:  
- Confirm the pathogen identity  
- Uncover antimicrobial resistance (AMR) profiles  
- Detect toxins driving mortality  
- Suggest **treatment strategies**  

WGS provides the tools to unlock these genetic secrets.  

---

## Project Objectives  

As a bioinformatics team, your goals are to:  

1. **Confirm organism identity** from WGS data using tools like BLAST.  
2. **Determine antimicrobial resistance (AMR) profiles** of the isolates.  
3. **Detect toxin genes** that may accelerate fatal outcomes (e.g., *hly, plcA, plcB*).  
4. **Suggest evidence-based antibiotics/treatments** to manage the outbreak.  

This is more than a lab exercise—your findings could **guide public health responses** and prevent future tragedies.  

---

## Dataset  

You are provided with **100 bacterial isolate samples** collected during the **2017–2018 South African outbreak**.  
- You may downsample to **50 isolates** if needed.  

Download the dataset with:  

```bash
wget https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/SA_Polony_100_download.sh
bash SA_Polony_100_download.sh

```

---
## Project Rubric  

We expect a **code report submitted as a Markdown file** (a mix of code and text explanation).  
For reference, see example: [Scripting Example](https://github.com/josoga2/bash-course/blob/main/bash/module_7/scripting_1.md)  

### Evaluation Criteria  

| Criterion | Description | Points |
|-----------|-------------|--------|
| Organism Identification | Correctly identifies the organism(s) using BLAST | 1 |
| AMR Gene Detection | Accurately identifies antimicrobial resistance (AMR) genes using **abricate** or similar tools | 1 |
| AMR Profile Summary | Summarizes AMR profiles across isolates, including prevalence and resistance implications | 2 |
| Antibiotic Proposals | Proposes evidence-based antibiotics (e.g., **ampicillin + gentamicin**) or alternatives based on AMR profile | 1 |
| Report Quality | Provides a clear, well-organized report with methods, results, and public health discussion | 2 |
| Scripts & Documentation | Submits functional scripts (Bash/Python/R) in a GitHub repository with clear documentation | 2 |
| **Bonus**: Toxin Gene Detection | Identifies toxin genes (e.g., *hly, plcA, plcB*) | +2 |

### Communication Requirement  

In addition to the GitHub submission, you must provide a **LinkedIn Post** containing either:  
- A **video recording** discussing the project (problems, results, importance) in a way a layperson will understand, **or**  
- A **long-form written post** explaining the same.  
