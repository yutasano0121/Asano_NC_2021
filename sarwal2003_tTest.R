# Apply t-test to detect genes upregulated in rejected allografts 
# compared with normal allografts.

# directory names
workDir <- getwd()  # Or enter a path to your directory of choice.
inputDir <- file.path(workDir, "input")

# Load microarray data 
df <- read.csv(file.path(inputDir, 'sarwal2003.csv'), row.names=1)
# tissue type annotation
anno <- read.csv(file.path(inputDir, 'annotation_sarwal2003.csv'))

# aggregate rows with same gene symbols
# fetch gene symbols for each row
genes <- list(df$GENE_SYMBOL)

# remove unwanted columns
df <- df[, -1:-3]

# take a sum
df <- aggregate(df, genes, sum)

# make rownames gene symbols
rownames(df) <- df$Group.1
df <- df[, -1]

# make a function to apply t-test removing NA values
t.test_rm.na <- function(X, ANNO=anno, COMPARISON=c("Rejection", "Normal"))
{
    anno <- ANNO$Tissue_type[!is.na(X)]
    x <- X[!is.na(X)]
    result <- t.test(x[anno==COMPARISON[1]], x[anno==COMPARISON[2]])
    results <- c(result$p.value, result$statistic)
    names(results) <- c("p.value", "t")
    return(results)
}

# apply t-test for each gene
p.value <- apply(df, 1, t.test_rm.na)
t.result <- as.data.frame(t(p.value))

# calculate FDR
t.result$FDR <- p.adjust(t.result$p.value, "BH")

# fetch t-test result data frame for genes upregulated in rejected tissues
genesUpInRejection <- t.result[(t.result$t > 0) & (t.result$FDR <= 0.05), ]

# save the result
write.csv(
    genesUpInRejection, 
    file.path(inputDir, "genesUpInRejection_tTest.csv")
)
