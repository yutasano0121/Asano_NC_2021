library(ggplot2)
library(edgeR)
library(pheatmap)
library(viridis)

logImmgen <- TRUE  # if immgen data are log2-scaled
removeAhnakLike <- FALSE  # if calculate DEG scores without Ahnak-like genes
saveImage <- TRUE
num_cluster <- 6


# directory names
workDir <- getwd()  # Or enter a path to your directory of choice.
inputDir <- file.path(workDir, "input")
DGEDir <- file.path(workDir, "DGE")
RDataDir <- file.path(workDir, "RData")
imageDir <- file.path(workDir, "image/heatmap")


# Load functions.
source(file.path(workDir, 'firstCohort_functions.R'))

# Load the immgen data (only peripheral B cells, not average taken)
immgen <- read.csv(file.path(inputDir, "immgen_Bcells.csv"), row.names=1)
if (logImmgen){
    immgen <- log2(immgen + 1)
    timeStamp("immgen data are log2 scaled")
} else {timeStamp("immgen data are NOT log2 scaled")}

# B cell subsets of interest
cell_types <- c(
    "B.T1.Sp", "B.T2.Sp", "B.T3.Sp", "B.Fo.Sp", "B.Fo.PC",
    "B.GC.Sp","B.MZ.Sp", "B1a.Sp", "B1a.PC", "B1b.PC"
)

# Make a list of cell types.
cellTypes <- paste(colnames(immgen), collapse=" ")
cellTypes <- gsub("\\.\\d", "", cellTypes)
cellTypes <- strsplit(cellTypes, " ")

# Load the list of human/mouse gene symbols.
anno_genes <- read.csv(
    file.path(inputDir, "human_mouse_geneID.csv"),
    stringsAsFactors=FALSE
)

# Remove Ahnak-covariant genes if needed.
if (removeAhnakLike){
    ahnakLike <- read.csv(
        file.path(inputDir, "immgen_AhnakLikeGenes.csv"),
        stringsAsFactors=FALSE
    )
    anno_genes <- anno_genes[!anno_genes$Mouse.gene.name %in% ahnakLike$mouse, ]

    timeStamp("Ahnak-like genes were removed")
} else {timeStamp("DEG scores are calculated WITH Ahnak-like genes")}


# Load a list of gene clusters.
load(file.path(
    RDataDir,
    paste0("hc_4comps_", num_cluster, "clusters.RData")
))

# convert human gene symbols to mouse counterparts
mouse_genes_list <- sapply(
    simplify=FALSE,
    names(hc_list),
    function(name){
        # Fetch gene clusters.
        df <- hc_list[[name]]$geneClusters
        clusters <- unique(df$cluster)  # a list of cluster numbers.
        converted_genes_list <- sapply(
            simplify=FALSE,
            clusters,
            function(i){
                genes_in_cluster <- df$symbol[df$cluster==i]
                mouse_genes <- anno_genes$Mouse.gene.name
                human_genes <- anno_genes$Gene.name
                converted_genes <- mouse_genes[human_genes %in% genes_in_cluster]

                return(converted_genes)
            }
        )

        return(converted_genes_list)
    }
)


# Make a list of subset data frames
# and then calculate a sum of scaled expression values.
immgen_scores_list <- sapply(
    simplify=FALSE,
    names(mouse_genes_list),
    function(name){
        # Fetch gene clusters.
        mouse_genes_list_sub <- mouse_genes_list[[name]]
        scores_list <- sapply(
            simplify=FALSE,
            1:length(mouse_genes_list_sub),
            function(i){
                mouse_genes <- mouse_genes_list_sub[[i]]
                immgen_sub <- immgen[rownames(immgen) %in% mouse_genes, ]
                # Scale by row. !Note that this operation transpose the matrix!
                immgen_sub_scaled <- apply(immgen_sub, 1, function(x){scale(x, center=TRUE)})
                # Sum up scaled values by row.
                scores <- rowSums(immgen_sub_scaled)
                names(scores) <- names(immgen)

                return(scores)
            }
        )

        return(scores_list)
    }
)


make_score_df <- function(SCORE){
    score.mean <- aggregate(SCORE, by=cellTypes, FUN=mean)
    score.se <- aggregate(SCORE, by=cellTypes, FUN=function(x){sd(x) / sqrt(length(x))})
    score.df <- data.frame(
        'mean'=score.mean$x,
        'se'=score.se$x,
        'category'=factor(
            score.mean$Group.1,
            levels=cell_types
        )
    )
    # Reorder the data frame by 'category.'
    score.df <- score.df[order(score.df$category), ]
    return(score.df)
}

score_df_list <- sapply(
    simplify=FALSE,
    names(immgen_scores_list),
    function(name){
        scores_list <- immgen_scores_list[[name]]
        df_list <- sapply(
            simplify=FALSE,
            1:length(scores_list),
            function(i){
                make_score_df(scores_list[[i]])
            }
        )

        return(df_list)
    }
)


# make a data frame to plot heatmaps
df4plot_list <- sapply(
    simplify=FALSE,
    names(score_df_list),
    function(name){
        df_list <- score_df_list[[name]]
        dfs4plot <- sapply(
            simplify=FALSE,
            1:length(df_list),
            function(i){
                return(df_list[[i]]$mean)
            }
        )
        dfs4plot <- do.call(rbind, dfs4plot)
        # Scale again.
        dfs4plot <- as.data.frame(t(apply(dfs4plot, 1, scale)))
        colnames(dfs4plot) <- cell_types
        rownames(dfs4plot) <- 1:length(df_list)
        return(dfs4plot)
    }
)


if (saveImage){
    lapply(
        names(df4plot_list),
        function(name){
            df4plot <- df4plot_list[[name]]
            b <- quantile_breaks(df4plot)  # color breaks
            if (removeAhnakLike){
                png(
                    file.path(imageDir, paste0("immgen_DEGscore_heatmap_woAhnak_", name, ".png")),
                    width=5, height=3, units="in", res=600
                )
            } else {
                png(
                    file.path(imageDir, paste0("immgen_DEGscore_heatmap_", name, ".png")),
                    width=5, height=3, units="in", res=600
                )
            }

            pheatmap(
                df4plot,
                cluster_cols=FALSE,
                cluster_rows=FALSE,
                color=plasma(length(b) - 1),
                breaks=b,
                angle_col=315,
                fontsize=12,
                labels_col=colnames(df4plot)
            )

            dev.off()
        }
    )
}

timeStamp("All done!")
