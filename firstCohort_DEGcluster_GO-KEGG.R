# packages
library(ggplot2)
library(viridis)
library(clusterProfiler)  # version 3.12.0
library(org.Hs.eg.db)  # version 3.8.2

num_cluster <- 6

workDir <- getwd()  # Or enter a path to your directory of choice.
RDataDir <- file.path(workDir, "RData/")
imageDir <- file.path(workDir, "image/GO_KEGG/")


# Load functions.
source(file.path(workDir, "firstCohort_functions.R"))

# Load gene clusters.
load(paste0(RDataDir, "/hc_4comps_", num_cluster, "clusters.RData"))


enrich_list <- sapply(
    simplify=FALSE,
    names(hc_list),
    function(name){
        # gene names with cluster ID
        tag <- hc_list[[name]]$geneClusters
        num_cluster <- length(unique(tag$cluster))

        # divide it by cluster ID
        tag.list <- lapply(
            1:num_cluster,
            function(i){
                return(tag[tag$cluster==i,])
            }
        )
        # make IDs ENTREZID
        entrez.list <- lapply(
            1:num_cluster,
            function(i){
                d <- makeRownamesENTREZ(tag.list[[i]])
                return(rownames(d))
            }
        )

        eGOlist <- lapply(
            1:num_cluster,
            function(i){
                return(test.GOenrichment(entrez.list[[i]]))
            }
        )
        timeStamp(paste("enrichGO done for", name))

        eKEGGlist <- lapply(
            1:num_cluster,
            function(i){
                return(test.KEGGenrichment(entrez.list[[i]]))
            }
        )
        timeStamp(paste("enrichKEGG done for", name))

        return(list(GO=eGOlist, KEGG=eKEGGlist))
    }
)

enrichPlot_list <- sapply(
    simplify=FALSE,
    names(enrich_list),
    function(name){
        eGOlist <- enrich_list[[name]][["GO"]]
        eKEGGlist <- enrich_list[[name]][["KEGG"]]
        plotGO_list <- lapply(
            1:length(eGOlist),
            function(i){
                eGO_simplified <- simplify(eGOlist[[i]])  # remove redundant terms
                plot.enrichedTerms(eGO_simplified, showCategory=10)
            }
        )
        plotKEGG_list <- lapply(
            1:length(eKEGGlist),
            function(i){
                plot.enrichedTerms(eKEGGlist[[i]], showCategory=10)
            }
        )
        return(list(GO = plotGO_list, KEGG = plotKEGG_list))
    }
)


for (name in names(enrichPlot_list)){
    pGO_list <- enrichPlot_list[[name]][["GO"]]
    pKEGG_list <- enrichPlot_list[[name]][["KEGG"]]
    for (i in 1:length(pGO_list)){
        pGO <- pGO_list[[i]]
        pKEGG <- pKEGG_list[[i]]
        pngStore(pGO, paste0(imageDir, "/eGO_", name, "_", i, ".png"), 3, 6)
        pngStore(pKEGG, paste0(imageDir, "/eKEGG_", name, "_", i, ".png"), 3, 6)
    }
}
