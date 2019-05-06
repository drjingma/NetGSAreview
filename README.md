# NetGSAreview

We provide here the supplementary materials to the review paper `A comparative study of network-based pathway enrichment analysis methods', including 

 * **Datasets**: data used in the study, 
 * **RCode**: R code and documentation for reproducing the numerical results, and
 * **Results**: numerical results for all settings considered in the review paper.  
   
**Datasets**

This folder consists of the data used in the comparative study, including

 * `breastcancer2012_ready.rda`: data for the TCGA breast cancer study.
 * `prostatecancer2015_ready.rda`: data for the TCGA prostate cancer study. 
 * `metabolomic_KEGGgraph.rdata`: KEGG metabolic interactions used in the metabolomic study.
 * `metabolomic_ready.csv`: data for the metabolomic study. 

Detailed descriptions of the first two R objects are available in the **RCode** folder.

**RCode** 

This folder consists of the R code used in the comparative study. Specifically, we provide three documentation files 

 * TCGA_BCa.html
 * TCGA_PCa.html
 * metabolomics.html

which, respectively, corresponds to data analysis for the TCGA breast cancer, prostate cancer and the metabolomic study. One can also find the corresponding R markdown files where the R code is available. All R functions needed are available in the file `preprocess_lib.r`. Finally, the two files `BCA_KEGGSPIA.RData` and `PCA_KEGGSPIA.RData` were prepared for running the Signaling Pathway Impact Analysis (SPIA) in each of the two cancer studies.

**Results**

This folder consists of the following files:
 
 * `BCA_t1error.csv`
 * `BCA_power_betweenness.csv`
 * `BCA_power_community.csv`
 * `BCA_power_neighborhood.csv`

which, respectively, provides detailed type I error and power results under different designs from our analysis of the breast cancer data. For example, the file `BCA_power_betweenness.csv` provides the empirical powers with the betweenness pathway dysregulation design for the breast cancer study, where the reader can find the empirical power of each pathway from different methods (`method`), with different levels of mean signals (`mu`) and permuted or original sample labels (`perm`). Several methods such as SPIA, Pathway-Express, CePa, PRS, and topologyGSA only work in selected pathways and at selected mean signals. In cases where these methods do not apply, their results are coded as `NA`.  

Similarly, one can find the type I error and power results for the prostate cancer study in the following files:

 * `PCA_t1error.csv`
 * `PCA_power_betweenness.csv`
 * `PCA_power_community.csv`
 * `PCA_power_neighborhood.csv`

Lastly, the type I error and power results for the metabolomics data example are available in:

 * `metabolomics_t1error.csv`
 * `metabolomics_power.csv`







