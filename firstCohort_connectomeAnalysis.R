library(igraph)
library(RCy3)

# directory names
workDir <- getwd()  # Or enter a path to your directory of choice.
inputDir <- file.path(workDir, "input")
RDataDir <- file.path(workDir, "RData")
DEGDir <- file.path(workDir, "DEtest_out")

# genes upregulated in rejected renal allografts (Sarwal et al., 2003)
tissue <- read.csv(
    file.path(inputDir, 'genesUpInRejectedTissue_tTest.csv'),
    row.names=1
)

# DEGs between kidney and tonsil B cells (gene detection rate filter: freq_threshold = 0.1)
scBcells <- read.csv(
    file.path(DEGDir, 'firstCohort_edgeR_switched_Tonsil-Kidney.csv'),
    row.names=1
)

# fetch names of genes upregulated in kidney B cells
scBcells <- rownames(scBcells)[(scBcells$logFC < 0) & (scBcells$threshold==TRUE)]

# FANTOM5 list of ligand-receptor pairs
ligRec <- read.csv(
    file.path(inputDir, 'PairsLigRec.tsv'),
    sep="\t"
)

# remove rows with a negative pair evidence
ligRec <- ligRec[grep("EXCLUDED", ligRec$Pair.Evidence, invert=TRUE), ]

# ligands and receptors in kidney-upregulated genes
scBcells.lig <- scBcells[scBcells %in% ligRec$Ligand.ApprovedSymbol]
scBcells.rec <- scBcells[scBcells %in% ligRec$Receptor.ApprovedSymbol]

# Subset the ligand-receptor list with those present in scBcells.
lig_inBcells <- ligRec$Ligand.ApprovedSymbol %in% scBcells
rec_inBcells <- ligRec$Receptor.ApprovedSymbol %in% scBcells
ligRec.sub <- ligRec[lig_inBcells | rec_inBcells, ]

# ligands and receptors present in both tissue and scBcells.
tissue.lig <- rownames(tissue)[rownames(tissue) %in% ligRec.sub$Ligand.ApprovedSymbol]
tissue.rec <- rownames(tissue)[rownames(tissue) %in% ligRec.sub$Receptor.ApprovedSymbol]

scBcells.ligRec <- lapply(
    c(scBcells.lig, scBcells.rec),
    function(x){return(c("Intra-graft Bcells", x))}
)

tissue.ligRec <- lapply(
    c(tissue.lig, tissue.rec),
    function(x){return(c("Rejected allografts\n(Sarwal et al., 2003)", x))}
)

# Detect ligand-receptor pairs present between B cells and tissues.
all.lig <- c(scBcells.lig, tissue.lig)
all.rec <- c(scBcells.rec, tissue.rec)

fetch.ligRecPairs <- function(i)
{
    if ((ligRec.sub$Ligand.ApprovedSymbol[i] %in% all.lig) &
    (ligRec.sub$Receptor.ApprovedSymbol[i] %in% all.rec))
    {
        return(c(as.character(ligRec.sub$Ligand.ApprovedSymbol[i]),
                 as.character(ligRec.sub$Receptor.ApprovedSymbol[i])
            )
        )
    }
}

# CAUTION: produces a lot of null values
ligRecPairs <- lapply(1:dim(ligRec.sub)[1], fetch.ligRecPairs)


# make a data frame to make igraph object
df <- rbind(
    do.call(rbind, scBcells.ligRec),
    do.call(rbind, tissue.ligRec),
    do.call(rbind, ligRecPairs)
)

# make an igraph object
g <- graph.data.frame(df)
V(g)$ligRecPairs <- TRUE
V(g)$ligRecPairs[1 : length(scBcells.ligRec) + length(tissue.ligRec)] <- FALSE
V(g)$category <- "Source"
V(g)$category[V(g)$name %in% all.lig] <- "Ligand"
V(g)$category[V(g)$name %in% all.rec] <- "Receptor"
V(g)$Paired <- FALSE
V(g)$Paired[V(g)$name %in% as.character(do.call(rbind, ligRecPairs))] <- TRUE

# open Cytoscape
cytoscapePing()
createNetworkFromIgraph(g, "cell-tissue connectome")
