# A code to run t-tqqests on integrated data_ordered
# make sure that the comparison is correctly defined in the code!
# packages
library(edgeR)
library(viridis)
library(RColorBrewer)
library(ggplot2)


comparisons <- c(
    "TonsilS-KidneyS",
    "TonsilScluster-KidneyS",
    "DSAneg-DSApos"
) # for defining categories to compare and the output name
# accepted 'comparison' values: "TonsilScluster-KidneyS", "TonsilS-KidneyS", "DSAneg-DSApos"

saveResult <- TRUE


# gene filter by detection rate
freq.threshold <- 0.1

# significance threshold
fdr.threshold <- 0.05
diff.threshold <- 1

# directories
workDir <- getwd()  # Or enter a path to your directory of choice.
inputDir <- file.path(workDir, "input")
imageDir <- file.path(workDir, "image")
RDataDir <- file.path(workDir, "RData")
outDir <- file.path(workDir, "DEtest_out")


# load functions
source(file.path(workDir, "integration_functions.R"))


# import the processed data
load(file.path(RDataDir, 'integrated_SCT-ComBat.RData'))
counts <- counts.norm.cb


# make an annotation column of DSA info
DSApos <- c("Kidney 4", "Kidney 6", "Kidney 7")
anno$DSA <- anno$Patient %in% DSApos


# subset genes with detection rates
expressedGenes <- find_expressedGenes(counts, freq.threshold)
counts.sub <- counts[expressedGenes, ]

timeStamp(paste(
    "Genes with detection rate below",
    freq.threshold * 100,
    "% were removed.",
    dim(counts.sub)[1],
    "genes passed the filtering."
))

# make variables to subset the data by
kidney <- anno$Source == "Kidney"
tonsil <- !kidney
switched <- anno$classSwitch
notPlasma <- anno$clusterAssignment != "plasmaCells"
tonsilDistinctCluster <- anno$clusterAssignment == "tonsilSwitched"

tTest_results <- sapply(
    simplify=FALSE,
    comparisons,
    function(comparison){
        # store the subset parameters as variables
        if (comparison == "TonsilScluster-KidneyS")
        {
            category1 <- (kidney & switched) & notPlasma
            category2 <- (tonsil & switched) & tonsilDistinctCluster
        } else if (comparison == "TonsilS-KidneyS")
        {
            category1 <- (kidney & switched) & notPlasma
            category2 <- (tonsil & switched) & notPlasma
        } else if (comparison == "DSAneg-DSApos") {
            category1 <- (kidney & anno$DSA) & notPlasma
            category2 <- (kidney & !anno$DSA) & notPlasma
        } else {stop("Value Error: 'comparison' needs to be set correctly.")}

        # make a list of categories to compare
        category <- list(category1, category2)


        timeStamp("Run t-tests.")
        print("Category 1:")
        print(summary(anno[category1, ]))
        print("Category 2:")
        print(summary(anno[category2, ]))

        t.result <- applyTtest2df(counts.sub, category)
        df.tResult <- summaryTtest(t.result)

        # print out the result
        timeStamp("T-tests run.")
        print(paste(
            "Mean difference", diff.threshold, "and FDR", fdr.threshold,
            " were used as the significance threshold."
        ))
        print(paste(
            sum(df.tResult$threshold & (df.tResult$mean.difference < 0)),
            "genes were upregulated in the category 1."
        ))
        print(paste(
            sum(df.tResult$threshold & (df.tResult$mean.difference > 0)),
            "genes were upregulated in the category 2."
        ))

        if (saveResult)
        {
            outname <- paste0("integration_tTest_", comparison, ".csv")

            write.csv(df.tResult, file.path(outDir, outname))
            print(paste("df.tResult has been saved as", outname))
        }

        return(df.tResult)
    }
)
