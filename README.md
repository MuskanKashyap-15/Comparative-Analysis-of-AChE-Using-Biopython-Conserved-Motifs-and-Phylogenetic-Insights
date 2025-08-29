# Comparative-Analysis-of-Acetylcholinesterase-(AChE)-Using-Biopython-Conserved-Motifs-and-Phylogenetic-Insights

## 📌 Project Overview  
This project demonstrates the use of **Biopython** for analyzing **Acetylcholinesterase (AChE) protein sequences** from multiple species.  
The workflow covers **sequence retrieval, multiple sequence alignment, conserved motif identification, and phylogenetic tree construction** to study the evolutionary relationships of AChE.  

---

## 🎯 Objectives  
- Retrieve AChE protein sequences from **NCBI**.  
- Perform **Multiple Sequence Alignment (MSA)** using Clustal Omega.  
- Visualize **conserved motifs** using Jalview.  
- Construct and visualize a **phylogenetic tree** with Biopython and Matplotlib.  
- Interpret **evolutionary relationships** among selected species.  

---

## 🧬 Methods & Workflow  
1. **Sequence Retrieval**  
   - Protein sequences of AChE from different species (e.g., human, mouse, rat, fruit fly) retrieved from NCBI in FASTA format.  

2. **Multiple Sequence Alignment (MSA)**  
   - Alignment performed using **Clustal Omega** via Biopython.  
   - Output saved in `.aln` (Clustal) and `.phy` (PHYLIP) formats.  

3. **Motif Visualization**  
   - Alignment results visualized in **Jalview** to highlight conserved regions.  

4. **Phylogenetic Tree Construction**  
   - Used Biopython’s **Phylo module** with the `.phy` alignment file.  
   - Tree plotted and saved using **Matplotlib**.  

---

## 📊 Results  
- **Conserved motifs** detected across all species, indicating essential functional sites of AChE.  
- **Phylogenetic analysis**:  
  - Mammals (human, mouse, rat) clustered closely, showing high similarity.  
  - Fruit fly appeared as the most evolutionary distant species.  
- Demonstrates **conservation and divergence** of AChE sequences across evolution.  

---

## 🛠️ Tools & Technologies  
- **Biopython** – sequence handling, alignment, and phylogenetics.  
- **Clustal Omega** – multiple sequence alignment.  
- **Jalview** – alignment and motif visualization.  
- **Matplotlib** – phylogenetic tree plotting.  
- **NCBI** – sequence data source.  

---

## 🔬 Applications  
- Understanding **evolutionary conservation** of AChE.  
- Identifying conserved motifs as potential **drug target regions**.  
- Demonstrates a reproducible **bioinformatics pipeline** using Python.  

---

## 📷 Visuals  
- Multiple Sequence Alignment (Jalview screenshot).  
- Phylogenetic Tree (Matplotlib output).  

---

## 📖 References  
- Cock PJA et al. (2009). *Biopython: freely available Python tools for computational molecular biology and bioinformatics*. **Bioinformatics, 25(11):1422–1423.**  
- Sievers F & Higgins DG. (2018). *Clustal Omega for making accurate alignments of many protein sequences*. **Protein Science, 27(1):135–145.**  
- Waterhouse AM et al. (2009). *Jalview Version 2—a multiple sequence alignment editor and analysis workbench*. **Bioinformatics, 25(9):1189–1191.**  
- Hunter JD. (2007). *Matplotlib: A 2D graphics environment*. **Computing in Science & Engineering, 9(3):90–95.**  

---

## 🙌 Acknowledgements  
- **Biopython developers** for open-source tools.  
- **Clustal Omega & Jalview** for alignment visualization.  
- **NCBI** for sequence database access.  
- **Matplotlib** community for visualization support.  

