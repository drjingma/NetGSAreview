# NetGSAreview

We provide here the supplementary materials to the review paper `A comparative study of network-based pathway enrichment analysis methods', including 

 * the online supplementary material summarizing main findings of the comparative study,
 * numerical results for all settings considered in the review paper, and 
 * R code for reproducing the numerical results. 
   
**Datasets**

This folder consists of the data used in the comparative study, including

 * `BreastCancer.rda`: data for the TCGA breast cancer study.
 * `ProstateCancer.rda`: data for the TCGA prostate cancer study. 
 * `metabolomic_KEGGgraph.rdata`: KEGG metabolic interactions used in the metabolomics data example.
 * `metabolomic_data.csv`: data for the metabolomic study. 

Detailed description of the first two R objects are available in the **RCode** folder. 

**RCode** 

This folder consists of the R code used in the comparative study. Specifically, we provide three documentation files 

 * TCGA_BCa.html
 * TCGA_PCa.html
 * metabolomics.html

which, respectively, corresponds to data analysis for the TCGA breast cancer study, TCGA prostate cancer and the metabolomic study. 

The file `preprocess_lib.r` consists of all the R functions needed in the analysis. The file `hsa_KEGG_signalingSPIA.RData` was prepared for running the Signaling Pathway Impact Analysis (SPIA). 

**Results**

This folder consists of three subfolders:

 * BCa
 * PCa
 * Metabolomics

which, respectively, provides detailed numerical results from our analysis of the three data examples. For example, in the folder BCa, one can find the file `bca_powers_noperm_btw_t5_1.csv`, which are empirical powers with the betweenness pathway dysregulation design, without sample permutation and with mean signal 0.1 for the breast cancer study. 



