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

