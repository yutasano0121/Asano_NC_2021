library(ggplot2)
library(ggrepel)
library(Rtsne)


saveImage <- TRUE
num_cluster <- 6

# directory names
workDir <- getwd()  # Or enter a path to your directory of choice.
inputDir <- file.path(workDir, "input")
RDataDir <- file.path(workDir, "RData")
imageDir <- file.path(workDir, "image/scatter")


# load functions
source(file.path(workDir, 'firstCohort_functions.R'))

# AHNAK-covariant genes
ahnak_genes <- read.csv(
    file.path(inputDir, 'immgen_AhnakLikeGenes.csv'),
    stringsAsFactors=FALSE
)
ahnak_genes.human <- ahnak_genes$human[!is.na(ahnak_genes$human)]
ahnak_genes.mouse <- ahnak_genes$mouse

# immgen data
immgen <- read.csv(file.path(inputDir, "immgen_Bcells.csv"), row.names=1)
geneNum.immgen <- dim(immgen)[1]

# load clusters
# Load a list of gene clusters.
load(paste0(RDataDir, "/hc_4comps_", num_cluster, "clusters.RData"))

over_representation.list <- sapply(
    simplify=FALSE,
    names(hc_list),
    function(name){
        # gene names with cluster ID
        tag <- hc_list[[name]]$geneClusters
        num_cluster <- length(unique(tag$cluster))

        # divide it by cluster ID
        enrich.list <- lapply(
            1:num_cluster,
            function(i){
                test.enrichment(
                    tag$symbol[tag$cluster==i],
                    TARGET.LIST=ahnak_genes.human,  # 293 human homologues
                    TARGET.SIZE=length(ahnak_genes.mouse),  # use 333 murine genes to calculate the background ratio of Ahnak-like genes
                    BACKGROUND.SIZE=geneNum.immgen  # total mouse gene number in Immgen
                )
            }
        )

        # summarize over-representation test results into a data frame
        enrich.df <- as.data.frame(do.call(rbind, enrich.list))

        # make columns numeric
        for (col in colnames(enrich.df))
        {
            enrich.df[, col] <- as.numeric(enrich.df[, col])
        }

        return(enrich.df)
    }
)

plot.list <- sapply(
    simplify=FALSE,
    names(over_representation.list),
    function(name){
        # plot -log10 p-value and Ahnak-like gene ratio (%)
        enrich.df <- over_representation.list[[name]]
        p <- ggplot(
            enrich.df,
            aes(x=ratio*100, y=-log10(p.value))
        ) + geom_point(aes(size=size))

        p <- p + scale_size(limits=c(0, 1000), breaks=c(200, 400, 800))
        #p <- p + geom_text(label=rownames(enrich.df), size=3)
        p <- p + geom_text_repel(
            aes(label=rownames(enrich.df), size=300),
            point.padding=0.1, box.padding=0.3, force=10, segment.alpha=0
        )
        p <- p + xlab("Ahnak-like genes (%)") + ylab("-log10 p-value") + labs(size="Total gene count")
        p <- p + theme_classic()
        p <- p + theme.scatter
        p.enrich <- p

        return(p.enrich)
    }
)


if (saveImage)
{
    for (name in names(plot.list)){
        pngStore(
            plot.list[[name]],
            paste0(imageDir, "/AHNAKgeneEnrichment_", name, ".png"),
            HEIGHT=3.5, WIDTH=5
        )
    }
    timeStamp("Image saved!")
}


timeStamp("All done!")
