platform <- "home"  # or 'home'
processedData <- "firstCohort_afterQC-RUVg2.RData"
geneIDtype <- "geneSymbol"  # 'geneSymbol' or 'ESNG' or 'ESNT'

saveResult <- TRUE

# gene filter by detection rate
freq.threshold <- 0.1

# significance threshold
fdr.threshold <- 0.05
logFC.threshold <- 1

# packages
library(edgeR)
library(viridis)
library(RColorBrewer)
library(ggplot2)

# directories
workDir <- getwd()  # Or enter a path to your directory of choice.
inputDir <- file.path(workDir, "input")
imageDir <- file.path(workDir, "image")
RDataDir <- file.path(workDir, "RData")
outDir <- file.path(workDir, "DEtest_out")


# load functions
source(paste0(workDir, "firstCohort_functions.R"))

if (geneIDtype != "geneSymbol")
{
    # load gene name annotation
    anno.geneNames <- read.csv(paste0(inputDir, "geneAnnotation_human.GRCh38.csv"))
}


# import the processed data
load(paste0(RDataDir, processedData))

# round the counts
counts <- round(counts)

# remove ERCC reads
counts <- counts[grep("ERCC-", rownames(counts), invert=TRUE), ]


# subset counts and annotation
subset <- anno$classSwitch == TRUE
counts.sub <- counts[, subset]
anno.sub <- anno[subset, ]
W <- ruvg$W[subset, ]  # batch effect from the RUVg result

timeStamp("Data have been subset for edgeR")
print(summary(anno.sub))


# subset genes with detection rates
expressedGenes <- find_expressedGenes(counts.sub, freq.threshold)
counts.sub.sub <- counts.sub[expressedGenes, ]

timeStamp(paste0("Genes with detection rate below ", freq.threshold * 100, "% were removed. ", dim(counts.sub.sub)[1], " genes passed the filtering."))


# run edgeR
timeStamp("EdgeR started.")
comparison <- anno.sub$Source
design <- model.matrix(~anno.sub$Source + W)
y <- DGEList(counts=counts.sub.sub, group=anno.sub$Source)
timeStamp("Cauculating normalization factors.")
y <- calcNormFactors(y)
timeStamp("Estimating Dispersions.")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
timeStamp("Fit a generalized linear model.")
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
tag <- topTags(lrt, n=dim(counts)[1])
tag <- tag$table
timeStamp("edgeR done.")

# add a column of significance threshold
tag$threshold <- (tag$FDR <= fdr.threshold) & (abs(tag$logFC) >= logFC.threshold)
gene_symbol <- c()

# add a column of gene symbols
if (geneIDtype != "geneSymbol")
{
    for (item in rownames(tag))
    {
        newSymbol <- anno.geneNames$gene_symbol[anno.geneNames$ENSG ==item]
        gene_symbol <- c(gene_symbol, as.character(newSymbol))
    }
    tag$gene_symbol <- gene_symbol
} else {tag$gene_symbol <- rownames(tag)}

# print edgeR summary
print(paste("edgeR has been run:",
            levels(comparison)[2],
            "-",
            levels(comparison)[1]))

print(paste("FDR threshold:",
            fdr.threshold,
            "logFC threshold:",
            logFC.threshold))

print(paste(sum(tag$threshold[tag$logFC > 0]),
            "genes upregulated in",
            levels(comparison)[2]))

print(paste(sum(tag$threshold[tag$logFC < 0]),
            "genes upregulated in",
            levels(comparison)[1]))

if(saveResult == TRUE)
{
    outname <- paste0(outDir,
                      "firstCohort_edgeR_freq", freq.threshold,
                      "_fdr", fdr.threshold,
                      "_logFC", logFC.threshold,
                      levels(comparison)[2], "-", levels(comparison)[1],
                      ".csv")
    write.csv(tag, outname)
    print(paste0("tag has been saved as", outname))
}
