# packages
library(RUVSeq)
library(edgeR)
library(viridis)
library(RColorBrewer)
library(ggplot2)

# parameters for the analysis
saveImage <- TRUE  # if save QC a plot
saveResult <- TRUE  # if store QC/RUVg'ed data as an RData
RUVgParameter = 2  # RUVG k parameter, default 2
out_name <- paste0("firstCohort_afterQC-RUVg", RUVgParameter, ".RData")


# directories
workDir <- getwd()  # Or enter a path to your directory of choice.
inputDir <- file.path(workDir, "input")
imageDir <- file.path(workDir, "image")
RDataDir <- file.path(workDir, "RData")


# Load functions.
source(file.path(workDir, "firstCohort_functions.R"))


# Load gene expression data.
timeStamp("Data loading started.")

counts <- read.csv(file.path(inputDir, "kallisto_rawCounts_allCells_geneSymbol.csv"), row.names = 1)  # gene x cell
anno <- read.csv(file.path(inputDir, "annotation.csv"))


# Subset the first cohort.
counts <- counts[, anno$Cohort == "First"]
anno <- anno[anno$Cohort == "First", ]


# Remove cells with zero counts for endogenous genes.
notZeroCounts <- colSums(counts[grep("ERCC-", rownames(counts), invert = TRUE), ]) != 0
counts <- counts[, notZeroCounts]
anno <- anno[notZeroCounts, ]


timeStamp("Data loaded.")


# QC the data
# Count the number of mapped reads.
librarySize <- colSums(counts)
anno$librarySize <- librarySize


# Count the number of expressed genes.
geneNum <- colSums(counts != 0)
anno$gene.number <- geneNum


# Count the sum cpm of expressed constant region genes.
HC <- grep("IGH[GA]\\d|IGH[DME]$",
           rownames(counts),
           perl=TRUE, value=TRUE)

# Fetch most highly expressed Ig isotype.
anno$Ig.class <- apply(counts[HC, ], 2, function(x){return(HC[which.max(x)])})

anno$classSwitch <- !anno$Ig.class %in% c("IGHM", "IGHD")

# Calculate the sum of IGHC gene counts (cpm).
HCsum.cpm <- colSums(cpm(counts)[rownames(counts) %in% HC, ])
anno$HCsum <- log(HCsum.cpm + 1)

# Plot a QC histograms.
p.libSize <- qplot(x = anno$Patient, y = log10(anno$librarySize))
p.geneNum <- qplot(geneNum, bins=60)
p.HC <- qplot(log(HCsum.cpm + 1), bins=60)
gridExtra::grid.arrange(p.libSize, p.geneNum, p.HC, nrow = 1)


# filter the data by qc criteria gene number cutoff is defined by the node of
# distribution (lower threshold) and 99.9% percentile (higher threshold) log2HCsum
# is defined by (roughly) boxplot's outlier definition
cutoff.geneNum.low <- 3000
cutoff.geneNum.high <- 15000
cutoff.HC <- 5
qc <- ((log(HCsum.cpm + 1) > cutoff.HC) & (geneNum > cutoff.geneNum.low)) & (geneNum < cutoff.geneNum.high)

# plot QC in a single plot
p <- ggplot(anno, aes(x=gene.number, y=HCsum, colour=qc)) + geom_point()
p <- p + scale_colour_manual(values = c("darkgrey", "red"))
p <- p + theme_classic() + xlab("Gene count") + ylab("Ig Heavy chain expression") + labs(colour=NULL)
p <- p + theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12, colour = "black"))
p.qc_scatter <- p + theme(legend.position="none")

if (saveImage){
    png(file.path(imageDir, 'bar_histo/firstCohort_qc_histo.png'), width=8, height=4, units='in', res=600)
    gridExtra::grid.arrange(p.libSize, p.geneNum, p.HC, nrow = 1)
    dev.off()
    pngStore(p.qc_scatter, file.path(imageDir, 'scatter/firstCohort_qc.png'), WIDTH=4, HEIGHT=4)
    print("QC plot saved.")
}


# subset by QC criteria
counts <- counts[, qc]
anno <- anno[qc, ]

# remove unexpressed genes
counts <- counts[rowSums(counts) > 0, ]

timeStamp(paste("QC done.", dim(counts)[2], "cells with", dim(counts)[1], "genes passed. Annotation summary \n"))
print(summary(anno))


# BATCH CORRECTION fetch ERCC genes
ercc <- grep("ERCC-", rownames(counts), value = TRUE)
erccCounts.beforeRUVg <- counts[ercc, ]
endoGenes <- grep("ERCC-", rownames(counts), value = TRUE, invert = TRUE)

timeStamp(paste("run RUVg. k:", toString(RUVgParameter)))
ruvg <- RUVg(round(as.matrix(counts)), ercc, k = RUVgParameter)
timeStamp("RUVg done.")


if (saveResult == TRUE) {
    # list of QC cutoffs
    qc.list <- list(cutoff.geneNum.low = cutoff.geneNum.low,
                    cutoff.geneNum.high = cutoff.geneNum.high,
                    cutoff.HC = cutoff.HC
                )
    save(counts, anno, ruvg, RUVgParameter, qc.list, file = file.path(RDataDir, out_name))
}
