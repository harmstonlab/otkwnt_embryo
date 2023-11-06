<!-- ⚠️ This README has been generated from the file(s) "blueprint.md" ⚠️--><p align="center">
  <img src="readme_images/ccser1_title.png" alt="Logo" width="auto" height="auto" />
</p>
<h1 align="center">otkwnt_embryo</h1>
<p align="center">
  <b>RNA-seq analysis of Otk1, Wnt4, Otk2 Drosophila embryos</b></br>
  <sub><sub>
</p>

<br />



<p align="center">
		<a href="undefined"><img alt="R 4.2.2" src="https://img.shields.io/badge/R-4.2.2-blue" height="20"/></a>
<a href="undefined"><img alt="Quarto v1.3.340" src="https://img.shields.io/badge/Quarto-1.3.340-orange" height="20"/></a>
<a href="undefined"><img alt="Status:in-progress" src="https://img.shields.io/badge/status-in--progress-yellow" height="20"/></a>
	</p>


[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/vintage.png)](#table-of-contents)

## ➤ Table of Contents

* [➤ 1. About The Project](#-1-about-the-project)
* [➤ 2. Introduction](#-2-introduction)
* [➤ 3. Dataset](#-3-dataset)
* [➤ 4. Preprocessing](#-4-preprocessing)
* [➤ 5. Quality Control](#-5-quality-control)
* [➤ 6. Differential Expression Analysis](#-6-differential-expression-analysis)
* [➤ 7. Clustering](#-7-clustering)
* [➤ 8. KEGG pathview plots](#-8-kegg-pathview-plots)
* [➤ 8. Results](#-8-results)
* [➤ 9. MD5 checksums](#-9-md5-checksums)
* [➤ 10. Version info ](#-10-version-info-)
* [➤ 11. Acknowledgements](#-11-acknowledgements)



[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/vintage.png)](#1-about-the-project)

## ➤ 1. About The Project
[*Back to Top*](#-table-of-contents)

This project analyzes bulk RNA-seq data from various **Otk1**, **Otk2** and **Wnt4** mutants in Drosophila embryos.  

Specifically, we are interested in the following questions: 

1. Does Wnt4 bind to Otk1 and activate specific pathways? 
2. What is the relationship between Otk1 and Otk2? 


[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/vintage.png)](#2-introduction)

## ➤ 2. Introduction
[*Back to Top*](#-table-of-contents)

**Otk1** (Off-track 1) is a transmembrane protein that is associated with colon cancer, but its mechanism of action is still unknown. **Wnt4** has been shown to bind to **Otk1**, but whether this binding is functional and activates any particular pathway has not been studied. 

In this study, we overexpress and downregulate various combinations of these three factors (**Otk1**, **Otk2** and **Wnt4**) and analyze the resulting bulk RNA-seq data. The UAS system was used to upregulate genes, while RNAi was used to downregulate target gene expression. In particular, we are looking for specific pathways that are either activated/ repressed by the Wnt4-Otk1-Otk2 signaling axis. 


[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/vintage.png)](#3-dataset)

## ➤ 3. Dataset
[*Back to Top*](#-table-of-contents)

**Mutation details:**

Label   | Condition |  Description 
---     | -------   | -------------------------------------
Ctrl    | Control              | Wild-type Drosophila embryos
W4_OE   | UAS-wnt4             | Drosophila embryos overexpressing Wnt4
O1_OE   |  UAS-otk1             | Drosophila embryos that overexpress Otk1
O1W4_OE | UAS-otk1, UAS-wnt4 | Drosophila embryos overexpressing both Otk1 and Wnt4
O2_KO   | otk2RNAi             | Drosophila embryos downregulating Otk2

________________________________________________________________________________________________


[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/vintage.png)](#4-preprocessing)

## ➤ 4. Preprocessing

Initial quality control was performed with [fastqc](scripts/fastqc_star_rsem/00_fastqc.sh). Reads were then aligned with [STAR](scripts/fastqc_star_rsem/01_STAR.sh), and counted with [RSEM](scripts/fastqc_star_rsem/02_RSEM.sh). 

Data was of high quality overall, and no issues were flagged during the QC process. 

MultiQC report: 
- View on Github:  [data/02_aligned/03_multiqc/multiqc_report.pdf](data/02_aligned/03_multiqc/multiqc_report.pdf)

- Download interactive copy: [data/02_aligned/03_multiqc/multiqc_report_interactive.html](data/02_aligned/03_multiqc/multiqc_report_interactive.html)


[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/vintage.png)](#5-quality-control)

## ➤ 5. Quality Control
Folder: [analysis/01_QC/](analysis/01_QC/)

*Contents: Read counts per sample, gene counts per sample, correlation heatmaps, PCA plots.* 


[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/vintage.png)](#6-differential-expression-analysis)

## ➤ 6. Differential Expression Analysis
Folder: [analysis/02_DE/](analysis/02_DE/)

*Contents: Volcano plots, MA plots, log2FoldChanges compared to control. Also contains GO and KEGG for upregulated and downregulated genes.* 


[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/vintage.png)](#7-clustering)

## ➤ 7. Clustering
Folder: [analysis/03_clust/](analysis/03_clust/)

*Contents: K-means clustering, GO and KEGG enrichments. Also contains silhouette plots and WSS plots.* 


[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/vintage.png)](#8-kegg-pathview-plots)

## ➤ 8. KEGG pathview plots
Folder: [analysis/03b_keggplot/](analysis/03b_keggplot/)

*Contents: KEGG pathview plots for selected pathways. * 


[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/vintage.png)](#8-results)

## ➤ 8. Results
[*Back to Top*](#-table-of-contents)

- Master DEG list (**LRT p-values**, **log2FoldChange**, **padj**, **cluster**) [analysis_all/output/04_compile/master_de_gene_list.xlsx](analysis_all/output/04_compile/master_de_gene_list.xlsx)



[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/vintage.png)](#9-md5-checksums)

## ➤ 9. MD5 checksums
[*Back to Top*](#-table-of-contents)

- `Drosophila_melanogaster.BDGP6.22.97.chr.gtf`: cba22b0ac8e925a53f3e60f450d4888a
- `dm6.fa`: 5aadf7ccab5a6b674e76516bf75eaa09



[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/vintage.png)](#10-version-info-)

## ➤ 10. Version info 
[*Back to Top*](#-table-of-contents)

- **FastQC:** v0.11.8
- **STAR:** 2.7.1a
- **RSEM:** v1.3.1
- **MultiQC:** version 1.14 (931327d)
- **Drosophila genome:** `dm6.fa`, `Drosophila_melanogaster.BDGP6.22.97.chr.gtf`


[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/vintage.png)](#11-acknowledgements)

## ➤ 11. Acknowledgements
This table of contents was built with [https://github.com/andreasbm/readme](https://github.com/andreasbm/readme)




