# packages
library(RUVSeq)
library(edgeR)
library(viridis)
library(RColorBrewer)
library(ggplot2)

saveResult <- TRUE  # store QC/RUVg'ed data

# for t-SNE
genenum4tsne <- 1000  # number of genes used for normalization, default 2000
seed <- 0

# directories
workDir <- getwd()  # Or enter a path to your directory of choice.
inputDir <- file.path(workDir, "input")
imageDir <- file.path(workDir, "image/scatter")
RDataDir <- file.path(workDir, "RData")


# load functions
source(file.path(workDir, "integration_functions.R"))

# load the integrated data
load(file.path(RDataDir, "integrated_SCT-ComBat.RData"))

# load the Seurat-assigned clusters and append it to the annotation data frame
#load(paste0(RDataDir, "integrated_seuratClusters.RData"))
#anno$seuratCluster <- seuratClusters

# make t-SNE plots with subset of genes
counts.norm.cb.sub <- subset_genes(counts.norm.cb, genenum4tsne)

# run t-SNE
timeStamp("Run t-SNE.")
set.seed(seed)

tsne <- Rtsne::Rtsne(t(counts.norm.cb.sub))
d.tsne <- as.data.frame(tsne$Y)

p.patient <- plot.tsne(d.tsne, anno)
p.ig <- plot.tsne(d.tsne, anno, COLOUR="classSwitch", COLOUR.LAB="Ig class switch")
p.cluster <- plot.tsne(d.tsne, anno, COLOUR="seuratCluster", COLOUR.LAB="Seurat cluster")

if (saveResult){
    png(
        pasta0(
            imageDir, "integrated_tSNE_SCTall_scaleData.png"),
            15, 5, units="in", res=300
        )
    )
    gridExtra::grid.arrange(p.patient, p.ig, p.cluster, nrow = 1)
    dev.off()
} else {gridExtra::grid.arrange(p.patient, p.ig, p.cluster, nrow = 1)}
