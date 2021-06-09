# packages
library(edgeR)  # Use version 3.26.8
library(viridis)
library(RColorBrewer)
library(ggplot2)


processedData <- "firstCohort_afterQC-RUVg2.RData"
saveResult <- TRUE
freq_threshold <- 0.1  # gene filter by detection rate
# significance threshold
fdr_threshold <- 0.05
logFC_threshold <- 1


# directory names
workDir <- getwd()  # Or enter a path to your directory of choice.
inputDir <- file.path(workDir, "input")
outDir <- file.path(workDir, "DEtest_out")
RDataDir <- file.path(workDir, "RData")


# load functions
source(file.path(workDir, "firstCohort_functions.R"))


# import the processed data
load(file.path(RDataDir, processedData))

# round raw counts
counts <- round(counts)

# remove ERCC reads
counts <- counts[grep("ERCC-", rownames(counts), invert=TRUE), ]


# Set variables to subset the data.
kidney <- anno$Source == 'Kidney'
switched <- anno$classSwitch == TRUE


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


# subset counts and annotation
subset_data.for_edgeR <- function(
    criterion, COUNTS=counts, ANNO=anno, W=ruvg$W,
    FREQ_THRESHOLD=freq_threshold
){
    counts.sub <- COUNTS[, criterion]
    anno.sub <- ANNO[criterion, ]
    W.sub <- W[criterion, ]  # batch effect from the RUVg result

    # Subset genes with detection rates.
    expressedGenes <- find_expressedGenes(counts.sub, FREQ_THRESHOLD)
    counts.sub.sub <- counts.sub[expressedGenes, ]
    print(
        paste(
            "Genes with detection rate below", FREQ_THRESHOLD * 100,
            "% were removed.", dim(counts.sub.sub)[1], "genes passed the filtering."
        )
    )

    # Return subset data as a list.
    out_list <- list(counts = counts.sub.sub, anno = anno.sub, W = W.sub)
    return(out_list)
}


# Subset the data.
d.kidney <- subset_data.for_edgeR(kidney)
d.tonsil <- subset_data.for_edgeR(!kidney)
d.switched <- subset_data.for_edgeR(switched)
d.unswitched <- subset_data.for_edgeR(!switched)


# run edgeR
run_edgeR <- function(input_list, feature_to_compare){
    counts <- input_list$counts
    anno <- input_list$anno
    W <- input_list$W

    # Set a factor to compare and a design matrix.
    comparison <- anno[, feature_to_compare]
    design <- model.matrix(~comparison + W)

    y <- DGEList(counts=counts, group=comparison)
    timeStamp("Cauculating normalization factors.")
    y <- calcNormFactors(y)

    timeStamp("Estimating Dispersions.")
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)

    timeStamp("Fit a generalized linear model.")
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, coef=2)

    # Fetch a result data table.
    tag <- topTags(lrt, n=dim(counts)[1])
    tag <- tag$table

    # add a column of significance threshold
    tag$threshold <- (tag$FDR <= fdr_threshold) & (abs(tag$logFC) >= logFC_threshold)


    timeStamp("edgeR done.")

    # print edgeR summary
    print(paste("edgeR has been run:",
                levels(comparison)[2],
                "-",
                levels(comparison)[1]))

    print(paste("FDR threshold:",
                fdr_threshold,
                "logFC threshold:",
                logFC_threshold))

    print(paste(sum(tag$threshold[tag$logFC > 0]),
                "genes upregulated in",
                levels(comparison)[2]))

    print(paste(sum(tag$threshold[tag$logFC < 0]),
                "genes upregulated in",
                levels(comparison)[1]))

    out_list <- list(tag=tag, comparison=levels(comparison))

    return(out_list)
}


tag.kidney.SvsU <- run_edgeR(d.kidney, 'classSwitch')
tag.tonsil.SvsU <- run_edgeR(d.tonsil, 'classSwitch')
tag.switched.KvsT <- run_edgeR(d.switched, 'Source')
tag.unswitched.KvsT <- run_edgeR(d.unswitched, 'Source')


save_tag <- function(tag_list, prefix){
    # 'prefix' is either 'kidney', 'tonsil', 'switched', or 'unswitched.'
    tag <- tag_list$tag
    comparison <- tag_list$comparison

    outname <- paste0(
        outDir,
        "/firstCohort_edgeR_", prefix, "_",
        comparison[2], "-", comparison[1],
        ".csv"
    )

    write.csv(tag, outname)

    print(paste("An edgeR result has been saved as", outname))
}

if(saveResult == TRUE)
{
    save_tag(tag.kidney.SvsU, 'kidney')
    save_tag(tag.tonsil.SvsU, 'tonsil')
    save_tag(tag.switched.KvsT, 'switched')
    save_tag(tag.unswitched.KvsT, 'unswitched')
}
