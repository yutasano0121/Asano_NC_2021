# Codes for single-cell RNA-seq analysis

## Folders
* `input`
    Input .csv files, such as read counts and annotation.
* `image`
    .png files of output plots.
* `RData`
    Analysis results stored as .RData files.
* `DEtest_out`
    .csv files of differential gene expression test results.

## Tested OS
* Ubuntu 16.04
* Mac Catalina

## Necessary CRAN R Packages
* ggplot2
* ggrepel
* gridExtra
* RColorBrewer
* viridis
* Rtsne
* pheatmap
* igraph
* dplyr

## Necessary Bioconductor Packages
* edgeR 3.26.8
* RUVSeq 1.16.1
* sva 3.32.1
* clusterProfiler 3.12.0
* org.Hs.eg.db 3.8.2
* Seurat 3.1.1
* RCy3 2.4.6

## External Data
* Immgen Microarray Data Phase 1 (GSE15907)
* Sarwal et al., NEJM 2003 (GSE343)  
*Not the entire data. Only 1340 transcripts identified as DEGs in the original paper.<br />
*Transcripts of the same gene name are aggregated, and outdated gene names are corrected.
* Ligand-receptor list from FANTOM5 (https://fantom.gsc.riken.jp/5/suppl/Ramilowski_et_al_2015/)
* Human-Mouse ortholog list from Ensembl BioMart (http://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/)

## Order of analyses
### First cohort
1. firstCohort_QC-RUVg.R...apply QC and batch normalization on the first cohort data
    * firstCohort_maketSNE.R...make t-SNE plots on the QC'ed cells  

    * firstCohort_makeViolin.R...make violin plots to visualize expression of single genes
2. firstCohort_edgeR.R...identify differentially expressed genes (DEGs) in the first cohort by edgeR
3. firstCohort_makeDEGcluster.R...cluster DEGs by their expression patterns
    * firstCohort_DEGcluster_pathwayScore.R...analyze expression level of gene signatures of interest e.g. innate immune genes
    * firstCohort_DEGcluster_GO-KEGG.R...apply GO/KEGG enrichment analysis on DEG clusters
4. correlation_analysis.R...test correlation of gene expression between cell populations; also applied to the integrated data

##### - Immgen and AHNAK analysis
1. immgen_fetchBcellData.py...fetch B cell expression data from Immgen data
2. immgen_AHNAKcorrelation.R...identify AHNAK-covariant genes
3. immgen_DEGclusterScore.R...test enrichment of AHNAK-covariant genes in DEG clusters

##### - Ligand-receptor connectome analysus
1. sarwal2003_tTest.R...identify rejection-upregulated genes in Sarwal et al. (2003) data  
2. firstCohort_connectomeAnalysis.R...identify co-upregulated ligand-receptor pairs

### Second cohort and integration
1. secondCohort_QC-RUVg.R...apply QC on the second cohort data
2. integration_SCTcombat.R...integrate the first and second cohort using Seurat/ComBat
    * integration_maketSNE.R...make t-SNE plots
3. integration_tTest.R...apply t-tests on integrated data
    * integration_makeVolcano.R...make volacno plots of t-test results
