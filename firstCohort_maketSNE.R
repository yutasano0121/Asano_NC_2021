# packages
library(Rtsne)
library(edgeR)
library(ggplot2)
library(viridis)


# parameters
saveImage <- TRUE
freq_threshold <- 0.1
feature_num <- 500
seed <- 0


# directory names
workDir <- getwd()  # Or enter a path to your directory of choice.
RDataDir <- file.path(workDir, "RData")
imageDir <- file.path(workDir, "image/scatter")

path.to.counts <- file.path(RDataDir, "firstCohort_afterQC-RUVg2.RData")


# load functions
source(file.path(workDir, "firstCohort_functions.R"))


# load the QC-RUVg'ed counts and annotation data frames
load(path.to.counts)


# Preprocess the data for t-SNE.
# 'counts' is genes x cells
process_counts <- function(
    counts,
    FREQ_THRESHOLD=freq_threshold,
    FEATURE_NUM=feature_num
){
    # Cpm-normalize the data
    counts_norm <-  as.data.frame(log2(cpm(counts) + 1))

    # Remove ERCC genes.
    endo_genes <- grep("ERCC-", rownames(counts), invert=TRUE)
    counts_norm <- counts_norm[endo_genes, ]

    # Remove genes with detection rate lower than 'FREQ_THRESHOLD.'
    expressedGenes <- find_expressedGenes(counts_norm, FREQ_THRESHOLD)
    counts_norm <- counts_norm[expressedGenes, ]  # subset the data
    print(
        paste0(
            "Genes with detection rate below ", FREQ_THRESHOLD * 100,
            "% were removed. ", dim(counts_norm)[1], " genes passed the filtering."
        )
    )

    # Subset genes with the highest variances.
    counts_out <- subset_genes(counts_norm, FEATURE_NUM)
    print(paste("Count data is subset with", FEATURE_NUM, "genes."))

    return(counts_out)
}


# Process the data
counts_norm <-  process_counts(counts)  # data w/o RUVg
counts_ruvg_norm <- process_counts(ruvg$normalizedCount)

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

timeStamp("Data processed.")


# run t-SNE
set.seed(seed)
tsne <- Rtsne(t(counts_norm))
d.tsne <- as.data.frame(-tsne$Y)
colnames(d.tsne) <- c('V2', 'V1')

set.seed(seed)
tsne_ruvg <- Rtsne(t(counts_ruvg_norm))
d.tsne_ruvg <- as.data.frame(-tsne_ruvg$Y)

timeStamp("t-SNE done. Make plots.")

tsne_colours <- list(
    Patient="Patient", classSwitch="Ig class switch",
    Batch="Batch", Ig.class="Ig class"
)

p.tsne_list <- sapply(
    simplify=FALSE,
    names(tsne_colours),
    function(name){
        plot.tsne(d.tsne_ruvg, anno, COLOUR=name, COLOUR.LAB=tsne_colours[[name]])
    }
)
p.tsne_list[['Batch_withoutRUVg']] <- plot.tsne(d.tsne, anno, COLOUR="Batch", COLOUR.LAB="Batch")

if (saveImage){
    lapply(
        names(p.tsne_list),
        function(name){
            pngStore(
                p.tsne_list[[name]],
                paste0(imageDir, '/tsne_', name, '.png'),
                WIDTH=4.5, HEIGHT=3
            )
        }
    )
}
