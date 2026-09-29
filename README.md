Overview

IL-15 and IL-2 share receptor subunits (IL2RB, IL2RG) and signal through JAK/STAT5, but they have different roles in immune regulation. 
This project examines how IL15 and related genes are expressed across breast tissue samples, and how expression varies by tissue type (normal vs. tumor) and metastasis status.

What the script does
Loads data: the FPKM expression matrix and the GEO series matrix metadata (GEOquery).
Cleans metadata: keeps the tissue and metastasis annotations and strips the label prefixes.
Reshapes: converts the expression matrix from wide to long format and joins it to the sample metadata.
Summarizes: calculates mean and median IL15 FPKM by tissue type.
Visualizes (ggplot2):
Bar plot of IL15 FPKM per sample
Density plot of IL15 by tissue type
Violin plot of IL15 by metastasis status
Scatter plot of IL15 vs. IL2 with linear fits per tissue type
Heatmap of 27 pathway genes across all samples (heatmap_save.pdf)
Genes analyzed
Cytokines and receptors: IL2, IL15, IL2RA, IL2RB, IL2RG, IL15RA
Signaling: JAK1, JAK3, STAT3, STAT5A, STAT5B
Downstream and functional genes: FOXP3, BCL2, BCL2L1, MCL1, MYC, CCND2, CCND3, PRF1, GZMB, IFNG, TBX21, EOMES, SOCS1, SOCS3, CISH, PRDM1
Key observations

Based on the heatmap, IL2 and most immune effector genes (IFNG, PRF1, EOMES, TBX21) show very low expression in these bulk tissue samples.
IL15 and IL15RA are low but detectable, and IL15 expression is higher in a small subset of samples. Survival and signaling genes such as MCL1, STAT3, and SOCS3 show the highest and most variable expression.
A few samples stand out with high MCL1 and SOCS3.

Requirements

R with readr, tidyverse, dbplyr, GEOquery, and ggplot2.

Usage

Download GSE183947_fpkm.csv.gz and GSE183947_series_matrix.txt.gz from GEO into the working directory, then run IL15_expression_in_BC.R.

Notes and possible next steps
The heatmap uses raw FPKM, so a few high-expression genes dominate the color scale. Using log2(FPKM + 1) or per-gene scaling would show the low-expressed genes better.
Statistical testing is not yet included. Paired tests for tumor vs. normal and correlation tests (IL15 vs. IL2 and receptors) would be logical additions.
Bulk RNA-seq reflects a mix of cell types, so IL15 signal partly reflects immune infiltration.
