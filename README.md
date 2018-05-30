# NetGSAreview

We provide here the supplementary materials to the review paper `A comparative study of network-based pathway enrichment analysis methods', including 

 * **Supplement**: the online supplementary material summarizing main findings of the comparative study,
 * **Results**: numerical results for all settings considered in the review paper, and 
 * **RCode**: R code for reproducing the numerical results. 
   
**Datasets**

This folder consists of the data used in the comparative study, including

 * `BreastCancer.rda`: data for the TCGA breast cancer study.
 * `ProstateCancer.rda`: data for the TCGA prostate cancer study. 
 * `metabolomic_KEGGgraph.rdata`: KEGG metabolic interactions used in the metabolomic study.
 * `metabolomic_data.csv`: data for the metabolomic study. 

Detailed descriptions of the first two R objects are available in the **RCode** folder.

**RCode** 

This folder consists of the R code used in the comparative study. Specifically, we provide three documentation files 

 * TCGA_BCa.html
 * TCGA_PCa.html
 * metabolomics.html

which, respectively, corresponds to data analysis for the TCGA breast cancer, TCGA prostate cancer and the metabolomic study. One can also find the corresponding R markdown files where the R code is available. All R functions needed are available in the file `preprocess_lib.r`. Finally, the file `hsa_KEGG_signalingSPIA.RData` was prepared for running the Signaling Pathway Impact Analysis (SPIA).

**Results**

This folder consists of three subfolders:

 * /BCa
 * /PCa
 * /Metabolomics

which, respectively, provides detailed numerical results from our analysis of the three data examples. For example, in the folder /BCa, one can find the file `bca_powers_noperm_btw_t5_1.csv`, which has the empirical powers with the betweenness pathway dysregulation design, without sample permutation and with mean signal 0.1 for the breast cancer study.

Note that SPIA and Pathway-Express only work when the number of DE genes is greater than zero. Threfore at mean signal 0.1, 0.2 and 0.3, the outputs of SPIA and Pathway-Express are all zero, meaning that they didn't return any values. In addition, at mean signal 0.4 and 0.5, when the number of DE genes is indeed greater than zero, SPIA and Pathway-Express only work for a subset of pathways due to topological constraints. Pathways that were analyzed by SPIA and Pathway-Express will have a numeric significance value, whereas the remaining ones not analyzed have `NA` values.

**Supplement**

This folder consists of the online supplementary material for the main paper, which summarizes the main findings from our comparative study across different pathway dysregulation designs and studies. It has provides a synthetic data example to illustrate the use of NetGSA in simultaneously detecting changes in means and networks. 




