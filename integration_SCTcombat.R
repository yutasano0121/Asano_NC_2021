# a code to integrate the first and second cohorts. make sure the parameters for SCTransform are set correctly.
# packages
library(dplyr)
library(Seurat)
library(edgeR)
library(sva)
library(Rtsne)
library(viridis)
library(RColorBrewer)
library(ggplot2)

saveResult <- TRUE

# for additional QC
removeHighMtRatio <- FALSE
removeLowLibrarySize <- FALSE

# for SCTransform
SCTsplitBy <- "All"  # 'Batch', 'Cohort', or 'All'
# data type to use after formalization. default 'data'
getDataType <- "counts"  # accepts either 'counts', 'data', or 'scale.data'
regressMtRatio <- FALSE

# for Seurat clustering
neighbor.dim <- 10  # number of PC to use. default 10
cluster.res <- 0.15  # default 0.5

# gene filter by detection rate
freq.threshold <- 0.1

# if save the integrated/normalized data
saveResult <- TRUE  # store QC'ed data

# directories
workDir <- getwd()  # Or enter a path to your directory of choice.
inputDir <- file.path(workDir, "input")
imageDir <- file.path(workDir, "image")
RDataDir <- file.path(workDir, "RData")

# paths to the dataset of each cohort
data.1 <- file.path(RDataDir, "firstCohort_afterQC-RUVg2.RData")
data.2 <- file.path(RDataDir, "secondCohort_afterQC.RData")




# load functions
source(file.path(workDir, "integration_functions.R"))


# import processed data
timeStamp("Data loading started.")

load(data.1)
load(data.2)
counts.1 <- as.data.frame(ruvg$normalizedCounts)
counts.2 <- round(counts.2)

# calculate mitochondrial gene ratio
anno$mt.ratio <- 100 * apply(counts.1, 2, function(x){sum(x[grep("MT-", rownames(counts))])}) / anno$library.size
anno.2$mt.ratio <- 100 * apply(counts.2, 2, function(x){sum(x[grep("MT-", rownames(counts.2))])}) / anno.2$library.size

cbind_all <- function(DF1, DF2, COMMON.ONLY=TRUE, FILTER.GENES=FALSE)
{
    if(FILTER.GENES)
    {
	DF1 <- DF1[find_expressedGenes(DF1, 0.05), ]
        DF2 <- DF2[find_expressedGenes(DF2, 0.05), ]
    }

    if (COMMON.ONLY)
    {
	common <- rownames(DF1)[rownames(DF1) %in% rownames(DF2)]
	DF1 <- DF1[common, ]
	DF2 <- DF2[common, ]
    } else {
	# combine QC/RUVg'ed first-cohort data and QC'ed second-cohort data
	# identify gene names present only in either of the datasets
	rownames.1.only <- rownames(DF1)[!rownames(DF1) %in% rownames(DF2)]
	rownames.2.only <- rownames(DF2)[!rownames(DF2) %in% rownames(DF1)]
	# fill missing rows with 0
	DF1[rownames.2.only, ] <- 0
	DF2[rownames.1.only, ] <- 0
	# make the two data have the same row order
	row.order.1 <- order(rownames(DF1))
	row.order.2 <- order(rownames(DF2))
	DF1 <- DF1[row.order.1, ]
	DF2 <- DF2[row.order.2, ]
    }
    # bind by columns
    counts <- cbind(DF1, DF2)

    return(counts)
}
counts <- cbind_all(counts.1, counts.2, FILTER.GENES=TRUE)

counts <- counts[grep("ERCC-", rownames(counts), invert=TRUE), ]  # remove ERCC reads

anno <- rbind(anno, anno.2)
rownames(anno) <- colnames(counts)

timeStamp("Data loaded.")


# plot library size and mt ratio
p <- ggplot(anno, aes(x=Cohort, y=log10(library.size))) + geom_boxplot()
p <- p + theme_classic() + xlab("") + ylab("log10 sequencing depth")
p.libSize <- p + theme(axis.text.x=element_text(colour="black", size=14),
               axis.text.y=element_text(colour="black", size=12),
               axis.title=element_text(size=14))

p <- ggplot(anno, aes(x=Cohort, y=mt.ratio)) + geom_boxplot()
p <- p + theme_classic() + xlab("") + ylab("Mitochondrial gene ratio (%)")
p.mtRatio <- p + theme(axis.text.x=element_text(colour="black", size=14),
               axis.text.y=element_text(colour="black", size=12),
               axis.title=element_text(size=14))



if (removeHighMtRatio)
{
    # remove cells with mitochondrial ratio > 8%
    counts <- counts[, anno$mt.ratio < 8]
    anno <- anno[anno$mt.ratio < 8, ]
    timeStamp("Cells with high mitochondrial gene ratio were removed.")
}

if (removeLowLibrarySize)
{
    counts <- counts[, log10(anno$library.size) > 4]
    anno <- anno[log10(anno$library.size) > 4, ]
    timeStamp("Cells with low library size were removed.")
}

# remove genes with low detection rate
expressedGenes <- find_expressedGenes(counts, freq.threshold)
counts.sub <- counts[expressedGenes, ]  # subset the data
timeStamp(paste(
    "Genes with detection rate below", freq.threshold * 100, "% were removed.",
    dim(counts.sub)[1], "genes passed the filtering."
))


# make a Seurat object
expression <- CreateSeuratObject(counts.sub, meta.data = anno)
expression <- PercentageFeatureSet(
    expression,
    pattern = "^MT-",
    col.name = "percent.mt"
)  # calculate mitochodrial gene ratio

if (SCTsplitBy != "All")
{
    expression.list <- SplitObject(expression, split.by = SCTsplitBy)  # split the data by 'Batch' or 'Cohort'
} else {expression.list <- list(expression)}

# run SCTransform
timeStamp("SCTransform run.")
for (i in 1:length(expression.list)) {
    if (regressMtRatio)
    {
    	expression.list[[i]] <- SCTransform(
            expression.list[[i]],
            vars.to.regress="mt.ratio",
    		do.center=FALSE,
    		return.only.var.genes=FALSE,
            verbose=TRUE
        )
    } else {
    	expression.list[[i]] <- SCTransform(
            expression.list[[i]],
			do.center=FALSE,
			return.only.var.genes=FALSE,
            verbose=TRUE
        )
    }
    # fetch count matrices
    expression.list[[i]] <- as.data.frame(GetAssayData(
        expression.list[[i]],
        assay="SCT",
        slot=getDataType
    ))

    # combine normalizes matrices
    if (i == 1)
    {
        counts.norm <- expression.list[[i]]
    } else {
        counts.norm <- cbind_all(counts.norm, expression.list[[i]])
    }
}
timeStamp("Data normalized by SCT.")

print(paste("Data from slot:", getDataType, "is retained."))
if (getDataType == "counts")
{
    counts.norm <- log2(counts.norm + 1)
    print("Data are log2 scaled.")
} else if (getDataType == "data"){
    print("Data are log scaled.")
} else {
    print("Data are Pearson's residuals.")
}

# remove genes with zero variance
counts.norm <- counts.norm[apply(counts.norm, 1, var) != 0, ]
print(paste("Genes:", dim(counts.norm)[1], "Cells:", dim(counts.norm)[2]))


# run ComBat
timeStamp("ComBat run.")
model.combat <-  model.matrix(~1, data=anno)
counts.norm.cb = ComBat(
    dat=as.matrix(counts.norm),
    batch=anno$Cohort,
    mod=model.combat,
    par.prior=TRUE,  # default TRUE but gave an error...
    prior.plots=FALSE
)

timeStamp("ComBat done.")


### assign clusters by Seurat ###
pbmc <- CreateSeuratObject(counts.norm.cb, meta.data = anno)
pbmc <- FindVariableFeatures(pbmc, selection.method='vst', nfeatures = 1500)  # if an error is returned, use selection.method='mean.var.plot'
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
DimHeatmap(pbmc, dims = 1:10, cells = 500, balanced = TRUE)
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.15)  # resolution 0.4 gives 5 nice clusters.
anno$seuratCluster <- as.factor(Idents(pbmc))

# assign which cluster represents which cell types
prdm1 <- aggregate(counts.norm.cb["PRDM1", ], list(anno$seuratCluster), mean)
plasmaCells <- which(prdm1$x==max(prdm1$x)) - 1
spred2 <- aggregate(counts.norm.cb["SPRED2", ], list(anno$seuratCluster), mean)
tonsilSwtiched <- which(spred2$x==max(spred2$x)) - 1
anno$clusterAssignment <- factor(rep("chunk", length(anno$seuratCluster)),
                                    levels=c("chunk", "plasmaCells", "tonsilSwitched"))
anno$clusterAssignment[anno$seuratCluster == plasmaCells] <- "plasmaCells"
anno$clusterAssignment[anno$seuratCluster == tonsilSwtiched] <- "tonsilSwitched"

print(paste("Seurat cluster for plasma cells:", plasmaCells))


if (saveResult == TRUE)
{
    save(counts.norm.cb, anno, file=file.path(RDataDir, 'integrated_SCT-ComBat.RData'))
    timeStamp("Result saved!")
}
