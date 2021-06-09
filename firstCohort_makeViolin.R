saveImage <- TRUE
jitter <- TRUE
# for t-SNE
seed <- 0
feature.num <- 500
freq.threshold=0.1

# packages
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(ggrepel)

# directory names

workDir <- getwd()  # Or enter a path to your directory of choice.
DGEDir <- file.path(workDir, "DEtest_out")
RDataDir <- file.path(workDir, "RData")
imageDir <- file.path(workDir, "image/violin")


path.to.tag <- file.path(DGEDir, "firstCohort_edgeR_freq0.1_fdr0.05_logFC1Tonsil-Kidney.csv")
path.to.clusterProfiler.out <- file.path(RDataDir, "firstCohort_clusterProfiler.RData")
path.to.counts <- file.path(RDataDir, "firstCohort_afterQC-RUVg2.RData")


# load functions
source(file.path(workDir, "firstCohort_functions.R"))


print("Load data.")



library(edgeR)

# load the normalied count and annotation data frames
load(path.to.counts)
counts <- as.data.frame(log2(cpm(ruvg$normalizedCount) + 1))
# transpose the count data frame for violin plots (without cpm)
counts4violin <- as.data.frame(t(counts))
counts <- counts[grep("ERCC-", rownames(counts), invert=TRUE), ]

timeStamp("Data loaded.")



expressedGenes <- find_expressedGenes(counts, freq.threshold)
counts <- counts[expressedGenes, ]  # subset the data
print(paste0("Genes with detection rate below ", freq.threshold * 100, "% were removed. ",
    dim(counts)[1], " genes passed the filtering."))


# modify annotation columns
cs.from <- c("TRUE", "FALSE")
cs.to <- c("Switched", "Unswitched")

for (i in 1:length(cs.from))
{
    anno$classSwitch <- gsub(pattern=cs.from[i],
                            replacement=cs.to[i],
                            x=anno$classSwitch)
}
anno$classSwitch <- factor(anno$classSwitch, levels=cs.to)



# ABC genes score
abcUp <- c("FCRL5", "ZEB2", "FCRL3", "ITGAX", "FCGR2B", "CD19", "TBX21", "IL2RA", "IL6R", "TNFRSF13B", "FAS", "SLAMF7", "IL10RA", "TLR7", "TBK1", "IRF4", "IL21R")
abcDown <- c("CR2", "CCR7", "TNFAIP3", "BACH2", "TCF7", "ETS1", "TRAF5", "CXCR5", "IRF8")

counts4violin[, "DN score"] <- calcScore(t(counts4violin), abcUp) - calcScore(t(counts4violin), abcDown)

tsne.colours <- list(Patient="Patient", classSwitch="Ig class switch", Batch="Batch", Ig.class="Ig class")
violin.genes <- c("AHNAK", "DN score", "KLF6", "RFX2", "ZBTB32")  # Add gene names to plot.

timeStamp("Make plots.")
p.violin.list <- sapply(
    violin.genes,
    USE.NAMES=TRUE,
    simplify=FALSE,
    function(gene){
        p <- plot.violin(
            counts4violin, anno, gene,
            FILL="classSwitch", FILL.LAB="Ig class switch",
            JITTER=jitter
        )
        if (gene == "DN score"){
            p <- p + theme(axis.title.y = element_text(face="plain"))
        }
        return(p)
    }
)

if (saveImage)
{
    for (name in names(p.violin.list)){
        if (jitter){
            outname <- paste0(imageDir, "/violin_", name, "_jitter.png")
        } else {
            outname <- paste0(imageDir, "/violin_", name, ".png")
        }
        pngStore(p.violin.list[[name]], outname, WIDTH=4.5, HEIGHT=3)
    }
    timeStamp("Images saved!")
}
