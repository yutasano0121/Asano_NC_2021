# packages
library(edgeR)
library(ggplot2)

saveResult <- TRUE

# directories
workDir <- getwd()  # Or enter a path to your directory of choice.
inputDir <- file.path(workDir, "input")
imageDir <- file.path(workDir, "image")
RDataDir <- file.path(workDir, "RData")
outDir <- file.path(workDir, "DEtest_out")


# load functions
source(paste0(workDir, "firstCohort_functions.R"))


### import and process old data ###
timeStamp("Data loading started.")
#load the gene expression data
counts <- read.csv(paste0(inputDir, "kallisto_rawCounts_allCells_geneSymbol.csv"), row.names = 1)  # gene x cell
anno <- read.csv(paste0(inputDir, "annotation.csv"))


# subset first cohort
counts <- counts[, anno$Cohort == "Second"]
anno <- anno[anno$Cohort == "Second", ]

# remove ERCC genes and cells with zero counts
endoGenes <- grep("ERCC-", rownames(counts), invert=TRUE)
counts <- counts[endoGenes, ]
nonZeroCounts <- colSums(counts) != 0
counts <- counts[, nonZeroCounts]
anno <- anno[nonZeroCounts, ]


timeStamp("Data loaded.")


# QC the data
# count the number of mapped reads
librarySize <- colSums(counts)
anno$library.size <- librarySize


# count the number of expressed genes
geneNum <- colSums(counts != 0)
anno$gene.number <- geneNum


# count the sum cpm of expressed constant region genes
HC <- grep(
    "IGH[GA]\\d|IGH[DME]$",
    rownames(counts),
    perl=TRUE,
    value=TRUE
)

# fetch most highly expressed Ig isotyoe
anno$Ig.class <- apply(
    counts[HC, ],
    2,
    function(x){return(HC[which.max(x)])}
)

anno$classSwitch <- !anno$Ig.class %in% c("IGHM", "IGHD")


# calculate the sum of IGHC gene counts (cpm)
HCsum.cpm <- colSums(cpm(counts)[rownames(counts) %in% HC, ])
anno$HCsum <- log(HCsum.cpm + 1)

# plot QC histograms
p.libSize <- qplot(x = anno$Patient, y = librarySize)
p.geneNum <- qplot(geneNum, bins=60)
p.HC <- qplot(log2(HCsum.cpm + 1), bins=60)
gridExtra::grid.arrange(p.libSize, p.geneNum, p.HC, nrow = 1)




# filter the data by qc criteria
cutoff.geneNum.low <- 1000
cutoff.geneNum.high <- Inf
cutoff.HC <- 5
qc <- (((log(HCsum.cpm + 1) > cutoff.HC) &
   	    (geneNum >= cutoff.geneNum.low)) &
	    (geneNum <= cutoff.geneNum.high))


# plot QC in a single plot
#p <- ggplot(anno, aes(x=gene.number, y=HCsum, size=library.size, colour=qc)) + geom_point()
#p <- p + scale_size_continuous(breaks=c(1E6,3E6,5E6,7E6)) + labs(size="Sequencing depth")
p <- ggplot(anno, aes(x=gene.number, y=HCsum, colour=qc)) + geom_point()
p <- p + scale_colour_manual(values = c("darkgrey", "red"))
p <- p + theme_classic() + xlab("Gene count") + ylab("Ig Heavy chain expression") + labs(colour=NULL)
p <- p + theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12, colour = "black"))
p.qc <- p + theme(legend.position="none")


# subset by QC criteria
counts <- counts[, qc]
anno <- anno[qc, ]

# remove unexpressed genes
counts <- counts[rowSums(counts) > 0, ]

timeStamp(paste("QC done.", dim(counts)[2], "cells with", dim(counts)[1], "genes passed. Annotation summary \n"))
print(summary(anno))

if (saveResult == TRUE)
{
    counts.2 <- counts
    anno.2 <- anno
    save(counts.2, anno.2, file=paste0(RDataDir, 'secondCohort_afterQC.RData'))
}
