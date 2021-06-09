library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(edgeR)

allDEGs <- FALSE  # analyze without focusing on specific clusters
num_cluster <- 6
targetClusters <- c(2, 3, 4)  # gene clusters to analyze
print(paste('Gene clusters:', targetClusters, 'will be analyzed.'))


workDir <- getwd()  # Or enter a path to your directory of choice.
inputDir <- file.path(workDir, "input")
RDataDir <- file.path(workDir, "RData")
imageDir <- file.path(workDir, "image/violin")



# Load functions.
source(file.path(workDir, "firstCohort_functions.R"))

# Load count data.
load(file.path(RDataDir, "firstCohort_afterQC-RUVg2.RData"))

# Log2-cpm normalize the data.
counts <- log2(cpm(ruvg$normalizedCounts) + 1)

# Modify "classSwitch" column in the annotation.
classSwitch <- gsub("TRUE", "Switched", anno$classSwitch)
classSwitch <- gsub("FALSE", "Unswitched", classSwitch)
anno$classSwitch <- factor(classSwitch, levels=c("Switched", "Unswitched"))


# load a list of DEGs
load(file.path(RDataDir, paste0("hc_4comps_", num_cluster, "clusters.RData")))
gGO_list <- sapply(
    simplify=FALSE,
    names(hc_list),
    function(name){
        # list of DEGs
        tag <- hc_list[[name]]$geneClusters
        # Fetch those in clusters of interest.
        tag <- tag[tag$cluster %in% targetClusters, ]
        tag <- makeRownamesENTREZ(tag)
        entrez <- rownames(tag)

        # groupGO
        ggo5 <- groupGO(
            entrez,
            OrgDb="org.Hs.eg.db",
            ont="BP",
            level=5,
            readable=TRUE
        )
        timeStamp(paste("groupGO done for", name))

        return(ggo5)
    }
)


# fetch a result data frame
genes_list = sapply(
    simplify=FALSE,
    names(gGO_list),
    function(name){
        return(
            list(
                innate = fetchGenes(gGO_list[[name]], "innate immune response"),
                #pos = fetchGenes(gGO_list[[name]], "positive regulation of innate immune response"),
                #neg = fetchGenes(gGO_list[[name]], "negative regulation of innate immune response")
                bacteria = fetchGenes(gGO_list[[name]], "response to molecule of bacterial origin")
            )
        )
    }
)



# a list of y-axis labels
ylab.list = list(
    innate = "innate immune response",
    pos = "positive innate regulation",
    neg = "negative innate regulation",
    innate.score = "innate immune score",
    bacteria = "response to molecule of bacterial origin",
    nod = "NOD-like receptor signaling pathway",
    tlr = "Toll-like receptor signaling pathway"
)

# calculate scaled sum for each geneset
scores_list <- sapply(
    simplify=FALSE,
    names(genes_list),
    function(name){
        sapply(
            simplify=FALSE,
            names(genes_list[[name]]),
            function(name2){
                return(calcScore(counts, genes_list[[name]][[name2]]))
            }
        )
    }
)


# make a df for plot
df4plot_list <- sapply(
    simplify=FALSE,
    names(scores_list),
    function(name){
        return(as.data.frame(scores_list[[name]]))
    }
)

timeStamp("make plots")
# make a plot
plot_list <- sapply(
    simplify=FALSE,
    names(scores_list),
    function(name){
        df4plot <- df4plot_list[[name]]
        plot_list2 <- sapply(
            simplify=FALSE,
            names(df4plot),  # column names
            function(name2){
                p <- plot.violin(
                    df4plot, anno, Y=name2,
                    FILL="classSwitch", FILL.LAB="Ig class switch", JITTER=TRUE
                )
                p <- p + ylab(ylab.list[[name2]])
                p <- p + theme(axis.title.y=element_text(size=14, face="plain"))
                return(p)
            }
        )
        return(plot_list2)
    }
)


# save plots
lapply(
    names(plot_list),
    function(name){
        plot_list2 = plot_list[[name]]
        lapply(
            names(plot_list2),
            function(name2){
                p = plot_list2[[name2]]
                if (!allDEGs){
                    pngStore(
                        p,
                        file.path(
                            imageDir,
                            paste0(
                                "violin_", name,
                                "_clusters", paste(targetClusters, collapse='-'),
                                "_", name2, ".png"
                            )
                        ), HEIGHT=3., WIDTH=4.5
                    )
                } else {
                    pngStore(
                        p,
                        paste0(
                            imageDir, "/violin_", name,
                            "_allClusters", "_", name2, ".png"
                        ), HEIGHT=3., WIDTH=4.5
                    )
                }
            }
        )
    }
)

timeStamp("all done.")


# Group data by patients.
plot_list <- sapply(
    simplify=FALSE,
    names(scores_list),
    function(name){
        df4plot <- df4plot_list[[name]]
        plot_list2 <- sapply(
            simplify=FALSE,
            names(df4plot),  # column names
            function(name2){
                p <- plot.violin(
                    df4plot, anno, Y=name2, X="Patient",
                    FILL="Source", FILL.LAB="Tissue source", JITTER=TRUE
                )
                p <- p + ylab(ylab.list[[name2]])
                p <- p + theme(axis.title.y=element_text(size=14, face="plain"))
                p <- p + theme(axis.text.x=element_text(angle=-45, hjust=0))
                return(p)
            }
        )
        return(plot_list2)
    }
)

lapply(
    names(plot_list),
    function(name){
        plot_list2 = plot_list[[name]]
        lapply(
            names(plot_list2),
            function(name2){
                p = plot_list2[[name2]]
                if (!allDEGs){
                    pngStore(
                        p,
                        paste0(
                            imageDir, "/violin_", name,
                            "_clusters", paste(targetClusters, collapse='-'),
                            "_", name2, "_patient.png"
                        ), HEIGHT=4., WIDTH=6
                    )
                } else {
                    pngStore(
                        p,
                        paste0(
                            imageDir, "/violin_", name,
                            "_allClusters", "_", name2, "_patient.png"
                        ), HEIGHT=4., WIDTH=6
                    )
                }
            }
        )
    }
)

timeStamp("all done.")
