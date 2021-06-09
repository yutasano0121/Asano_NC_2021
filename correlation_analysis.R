# test correlation between cell population

library(ggplot2)
library(edgeR)

geneNum <- 1000  # number of high-variance genes, default 2000
remove_first_from_integrated <- FALSE

# directory names
workDir <- getwd()  # Or enter a path to your directory of choice.
inputDir <- file.path(workDir, "input")
DGEDir <- file.path(workDir, "DGE")
RDataDir <- file.path(workDir, "RData")
imageDir <- file.path(workDir, "image/heatmap")

# Load functions
source(file.path(workDir, "firstCohort_functions.R"))

timeStamp("Load the 1st cohort data.")
load(file.path(RDataDir, "firstCohort_afterQC-RUVg2.RData"))
anno1 <- anno
counts1 <- log2(cpm(ruvg$normalizedCounts) + 1)

timeStamp("Remove ERCC genes.")
endo_genes <- grep("ERCC-", rownames(counts1), invert=TRUE)
counts1 <- counts1[endo_genes, ]

timeStamp("Remove genes with detection rate lower than 0.1.")
expressedGenes <- find_expressedGenes(counts1, 0.1)
counts1 <- counts1[expressedGenes, ]  # subset the data

# subset by variance
variance <- apply(counts1, 1, var)
counts1 <- counts1[order(-variance), ]
counts1 <- counts1[1:geneNum, ]
timeStamp(paste(dim(counts1)[1], "high-variant genes are selected."))
print(counts1[1:5, 1:5])

# set group labels
groups <- list(
    Kidney.Switched = (anno1$Source == "Kidney") & (anno1$classSwitch),
    Kidney.Unswitched = (anno1$Source == "Kidney") & (!anno1$classSwitch),
    Tonsil.Switched = (anno1$Source == "Tonsil") & (anno1$classSwitch),
    Tonsil.Unswitched = (anno1$Source == "Tonsil") & (!anno1$classSwitch)
)

# test mean expression values between groups
test_cor_meanExpression <- function(COUNTS, GROUPS){
    out <- sapply(
        simplify=FALSE,
        names(GROUPS),
        function(g1){
            group1 <- GROUPS[[g1]]
            test_out_list <- sapply(
                simplify=FALSE,
                names(GROUPS),
                function(g2){
                    group2 <- GROUPS[[g2]]
                    means1 <- apply(COUNTS[, group1], 1, mean)
                    means2 <- apply(COUNTS[, group2], 1, mean)
                    test_out <- cor.test(means1, means2)
                    test_out <- data.frame(
                        group1 = g1,
                        group2 = g2,
                        p.value = test_out$p.value,
                        conf.int1 = test_out$conf.int[2],
                        conf.int2 = test_out$conf.int[1],
                        r = test_out$estimate
                    )
                    return(test_out)
                }
            )
            return(do.call(rbind, test_out_list))
        }
    )

    return(do.call(rbind, out))
}

cor_out1 <- test_cor_meanExpression(counts1, groups)


plot_cor_heatmap <- function(DF_cor){
    DF_cor.melt <- reshape2::melt(DF_cor)
    colnames(DF_cor.melt)[1:2] <- c("group1", "group2")
    # introduce newline in group names. This makes columns "characters" instead of factors
    DF_cor.melt$group1 <- gsub("\\.", "\n", DF_cor.melt$group1)
    DF_cor.melt$group2 <- gsub("\\.", "\n", DF_cor.melt$group2)
    # reverse the order of the 2nd groups
    DF_cor.melt$group2 <- factor(
        DF_cor.melt$group2,
        levels=unique(DF_cor.melt$group2)[order(unique(DF_cor.melt$group2), decreasing=TRUE)]
    )
    # round numbers
    DF_cor.melt$value <- signif(DF_cor.melt$value, 2)
    r <- DF_cor.melt[DF_cor.melt$variable=="r", ]
    p <- ggplot(r, aes(x=group1, y=group2, fill=value)) + geom_tile()
    p <- p + scale_fill_gradient(high="#006164", low="white")
    p <- p + theme.scatter + xlab("") + ylab("") + labs(fill="r") 
    p <- p + theme(
        axis.text.x = element_text(angle = -45, hjust=0, size=14),
        axis.text.y = element_text(size=14)
    )
    p <- p + geom_text(data=r, aes(label=value))
    return(p)
}

p_cor1 <- plot_cor_heatmap(cor_out1) + scale_fill_gradient(
    high="#006164", low="white",
    breaks=seq(0.5, 1, by=0.25)
)



timeStamp("Next, test correlation within the integrated data.")

load(file.path(RDataDir, "integrated_SCT-ComBat.RData"))
anno2 <- anno
counts2 <- counts.norm.cb[anno2$seuratCluster != 2, ]  # remove plasma cells
timeStamp("Remove genes with detection rate lower than 0.1.")
expressedGenes2 <- find_expressedGenes(counts2, 0.1)
counts2 <- counts2[expressedGenes2, ]

# remove cells from the 1st cohort if needed
if (remove_first_from_integrated){
    counts2 <- counts2[, anno2$Cohort == "Second"]
    anno2 <- anno2[anno2$Cohort == "Second", ]
    print("Cells from the first-cohort are removed from the integrated data.")
}

# subset by variance
variance2 <- apply(counts2, 1, var)
counts2 <- counts2[order(-variance2), ]
counts2 <- counts2[1:geneNum, ]
timeStamp(paste(dim(counts2)[1], "high-variant genes are selected."))
print(counts2[1:5, 1:5])

groups2 <- list(
    Kidney.Switched = (anno2$Source == "Kidney") & (anno2$classSwitch),
    Kidney.Unswitched = (anno2$Source == "Kidney") & (!anno2$classSwitch),
    Tonsil.Switched = (anno2$Source == "Tonsil") & (anno2$classSwitch),
    Tonsil.Unswitched = (anno2$Source == "Tonsil") & (!anno2$classSwitch)
)

cor_out2 <- test_cor_meanExpression(counts2, groups2)

p_cor2 <- plot_cor_heatmap(cor_out2) + scale_fill_gradient(
    high="#006164", low="white",
    limits=c(0.7, 1),
    breaks=seq(0.7, 1, by=0.1)
)



timeStamp("Next, test correlation between the 1st and integrated data.")

# combine annotations
anno1$dataset <- rep("first", dim(anno1)[1])
anno2$dataset <- rep("integrated", dim(anno2)[1])
anno3 <- rbind(anno1, anno2[, colnames(anno1)])

# common genes with high variance
genes <- rownames(counts1)[rownames(counts1) %in% rownames(counts2)]
timeStamp(paste("Correlation between data before/after integration is tested for", length(genes), "genes."))
counts1_2 <- counts1[genes, ]
counts2_2 <- counts2[genes, ]
counts3 <- cbind(counts1_2, counts2_2)

# set groups by dataset
groups3_first <- list(
    Kidney.Switched = ((anno3$Source == "Kidney") & (anno3$classSwitch)) & anno3$dataset == "first",
    Kidney.Unswitched = ((anno3$Source == "Kidney") & (!anno3$classSwitch)) & anno3$dataset == "first",
    Tonsil.Switched = ((anno3$Source == "Tonsil") & (anno3$classSwitch)) & anno3$dataset == "first",
    Tonsil.Unswitched = ((anno3$Source == "Tonsil") & (!anno3$classSwitch)) & anno3$dataset == "first"
)
groups3_integrated <- list(
    Kidney.Switched = ((anno3$Source == "Kidney") & (anno3$classSwitch)) & anno3$dataset == "integrated",
    Kidney.Unswitched = ((anno3$Source == "Kidney") & (!anno3$classSwitch)) & anno3$dataset == "integrated",
    Tonsil.Switched = ((anno3$Source == "Tonsil") & (anno3$classSwitch)) & anno3$dataset == "integrated",
    Tonsil.Unswitched = ((anno3$Source == "Tonsil") & (!anno3$classSwitch)) & anno3$dataset == "integrated"
)

test_cor_meanExpression2 <- function(COUNTS, GROUPS1, GROUPS2){
    out <- sapply(
        simplify=FALSE,
        names(GROUPS1),
        function(g1){
            group1 <- GROUPS1[[g1]]
            test_out_list <- sapply(
                simplify=FALSE,
                names(GROUPS2),
                function(g2){
                    group2 <- GROUPS2[[g2]]
                    means1 <- apply(COUNTS[, group1], 1, mean)
                    means2 <- apply(COUNTS[, group2], 1, mean)
                    test_out <- cor.test(means1, means2)
                    test_out <- data.frame(
                        group_first = g1,
                        group_integrated = g2,
                        p.value = test_out$p.value,
                        r = test_out$estimate
                    )
                    return(test_out)
                }
            )
            return(do.call(rbind, test_out_list))
        }
    )

    return(do.call(rbind, out))
}

cor_out3 <- test_cor_meanExpression2(counts3, groups3_first, groups3_integrated)

p_cor3 <- plot_cor_heatmap(cor_out3) + scale_fill_gradient(
    high="#006164", low="white",
    limits=c(min(cor_out3$r), 1),
    breaks=seq(0.4, 1, by=0.2)
)

pngStore(p_cor1, file.path(imageDir, "cor1.png"), WIDTH=6, HEIGHT=4.5)
pngStore(p_cor2, file.path(imageDir, "cor2.png"), WIDTH=6, HEIGHT=4.5)
pngStore(p_cor3, file.path(imageDir, "cor3.png"), WIDTH=6, HEIGHT=4.5)