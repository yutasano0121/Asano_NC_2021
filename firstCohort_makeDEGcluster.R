saveImage <- TRUE
trimSize <- 0  # percentile to trim when calculating mean
num_cluster = 6

# packages
library(edgeR)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(pheatmap)


# directory names
workDir <- getwd()  # Or enter a path to your directory of choice.
DEGDir <- file.path(workDir, "DEtest_out")
RDataDir <- file.path(workDir, "RData")
imageDir <- file.path(workDir, "image/heatmap")

# Load functions.
source(file.path(workDir, "firstCohort_functions.R"))

# Load lists of differentially expressed genes.
DEG_path_list <- list(
    'kidney_U-S' = file.path(DEGDir, "firstCohort_edgeR_kidney_Unswitched-Switched.csv"),
    'tonsil_U-S' = file.path(DEGDir, "firstCohort_edgeR_tonsil_Unswitched-Switched.csv"),
    'switched_T-K' = file.path(DEGDir, "firstCohort_edgeR_switched_Tonsil-Kidney.csv"),
    'unswitched_T-K' = file.path(DEGDir, "firstCohort_edgeR_unswitched_Tonsil-Kidney.csv")
)

DEG_list <- sapply(
    simplify=FALSE,
    names(DEG_path_list),
    function(name){read.csv(DEG_path_list[[name]], row.names=1)}
)

# Fetch a list of differential gene names.
DEGs <- lapply(
    names(DEG_list),
    function(name){extract_DEGs(DEG_list[[name]])}
)

# make it a vector
DEGs <- unique(unlist(DEGs))
print(paste("Number of DEGs:", length(DEGs)))

# Load count data.
counts_path <- file.path(RDataDir, "firstCohort_afterQC-RUVg2.RData")
load(counts_path)

# Normalize and subset the count data frame with DEGs.
counts <- log2(cpm(ruvg$normalizedCounts) + 1)[DEGs, ]

print("Data loaded.")


# Make a category vector to aggregate the counts.
category <- anno$Source
category <- paste0(category, ".", anno$classSwitch)
for (i in 1:2){
    category <- gsub(
        c("TRUE", "FALSE")[i],
        c("Switched", "Unswitched")[i],
        category
    )
}
category <- factor(
    category,
    levels=c("Kidney.Switched", "Kidney.Unswitched", "Tonsil.Switched", "Tonsil.Unswitched")
)

# Take mean values in each category
counts_mean <- aggregate(
    t(counts),
    by=list(category),
    FUN=function(x){mean(x, trim=trimSize)}
)
rownames(counts_mean) <- counts_mean$Group.1
counts_mean <- t(counts_mean[, -1])  # make it gene x category


# Scale the data across cells to have unit variance
counts_scaled <- apply(counts_mean, 1, scale)  # The result is cells x genes
rownames(counts_scaled) <- colnames(counts_mean)

# Make it genes x cells.
counts_scaled <- t(counts_scaled)

# set quantile color breaks
col_break <- quantile_breaks(counts_scaled)

print("Data processed.")


print("Make a heatmap.")

# make a hierachical cluster of rows
#method_list <- c("complete", "ward.D", "ward.D2", "mcquitty", "median", "centroid")
method_list <- c("ward.D2")


hc_list <- sapply(
    method_list,
    function(method){
        hierarchicalCluster(counts_scaled, method, num_cluster)
    },
    simplify=FALSE
)

heatmap_list <- sapply(
    simplify=FALSE,
    method_list,
    function(method){
        if (saveImage){
            png(
                paste0(imageDir, "/heatmap4comps_", method, "_", num_cluster, "_clusters.png"),
                height=8, width=4.1, units="in", res=300
            )
            plot_heatmap(
                counts_scaled, hc_list[[method]]$hc_gene, hc_list[[method]]$mycl,
                num_cluster, col_break
            )
            dev.off()
        } else {
            plot_heatmap(
                counts_scaled, hc_list[[method]]$hc_gene, hc_list[[method]]$mycl,
                num_cluster, col_break
            )
        }
    }
)

# Save a list of clustering results.
save(
    hc_list, num_cluster,
    file=paste0(RDataDir, "/hc_4comps_", num_cluster, "clusters.RData")
)

print("Done.")
