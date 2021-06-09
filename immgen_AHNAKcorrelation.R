library(ggplot2)
library(ggrepel)
library(Rtsne)

saveImage <- FALSE

# directory names
workDir <- getwd()  # Or enter a path to your directory of choice.
inputDir <- file.path(workDir, "input")
DGEDir <- file.path(workDir, "DGE")
RDataDir <- file.path(workDir, "RData")
imageDir <- file.path(workDir, "image/heatmap")

source(file.path(workDir, 'firstCohort_functions.R'))

# load the immgen data (only peripheral B cells, not average taken)
immgen <- read.csv(file.path(inputDir, "immgen_Bcells.csv"), row.names=1)
immgen.scaled <- apply(immgen, 1, function(x){scale(x, center=FALSE)})

# make a list of cell types
cellTypes <- paste(colnames(immgen), collapse=" ")
cellTypes <- gsub("\\.\\d", "", cellTypes)
cellTypes <- strsplit(cellTypes, " ")


# aggregate by cell types
df <- aggregate(immgen.scaled, by=cellTypes, FUN=mean)

ahnak <- df$Ahnak
itgam <- df$Itgam
vim <- df$Vim
s100a6  <- df$S100a6
klf6 <- df$Klf6
rfx2  <- df$Rfx2
corrWithAhnak <- read.csv(file.path(inputDir, 'immgen_Bcells_correlationWithAhnak.csv'))
AhnakLikeMouseGenes <- as.character(corrWithAhnak$GeneSymbol[corrWithAhnak$correlationCoefficient >= 0.8])

rownames(df) <- df$Group.1
# fetch the number of genes
geneNum.immgen <- length(colnames(df))
df <- df[, AhnakLikeMouseGenes]

# calculate statistics
scaled.mean <- apply(df, 1, mean)
scaled.sd <- apply(df, 1, sd)
scaled.se<- apply(df, 1, function(x){sd(x)/sqrt(length(x))})
ci <- qnorm(0.975)*scaled.se

# make a data frame to plot
df <- data.frame(
    mean=scaled.mean,
    sd=scaled.sd,
    se=scaled.se,
    ci=ci,
    groups=names(scaled.mean),
    Ahnak=ahnak, Itgam=itgam, Vim=vim, S100a6=s100a6, Rfx2 = rfx2, Klf6=klf6
)

# reorder the 'groups' factor column
df$groups <- factor(
    df$groups,
    levels=c(
        "B.T1.Sp", "B.T2.Sp", "B.T3.Sp", "B.Fo.Sp", "B.GC.Sp",
        "B.MZ.Sp", "B1a.Sp", "B.Fo.PC", "B1a.PC", "B1b.PC"
    )
)

p <- ggplot(df, aes(groups, group=1))
p <- p + geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd, fill="grey"))
#p <- p + geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd))
p <- p + geom_line(aes(y=mean, colour="Mean")) + scale_fill_manual(values="grey")
p <- p + geom_line(aes(y=ahnak, colour="Ahnak"))
p <- p + scale_color_manual(values=c("#F8766D", "black"))
p <- p + xlab("") + ylab("Scaled expression") + labs(colour="", fill="")
p <- p + guides(fill=FALSE)  # remove the 'fill' legend
p <- p + theme_classic()
p.AhnakLikeMouseGenes <- p + theme(axis.text.x=element_text(angle=315, hjust=0, colour="black", size=12),
      axis.text.y=element_text(colour="black", size=12),
      axis.title=element_text(size=14),
      legend.text=element_text(size=14))


# barplots for Vim and Itgam
p <- ggplot(df, aes(x=groups, y=Ahnak)) + geom_bar(stat="identity")
p <- p + theme_classic() + xlab("")
p.ahnak <- p + theme(
    axis.text.x=element_text(angle=315, hjust=0, colour="black", size=12),
    axis.text.y=element_text(colour="black", size=12),
    axis.title=element_text(size=14, face="italic"),
    legend.text=element_text(size=14),
    plot.margin=margin(0.5, 0.5, 0, 0, "cm")
)
p <- ggplot(df, aes(x=groups, y=Vim)) + geom_bar(stat="identity")
p <- p + theme_classic() + xlab("")
p.vim <- p + theme(
    axis.text.x=element_text(angle=315, hjust=0, colour="black", size=12),
    axis.text.y=element_text(colour="black", size=12),
    axis.title=element_text(size=14, face="italic"),
    legend.text=element_text(size=14),
    plot.margin=margin(0.5, 0.5, 0, 0, "cm")
)
p <- ggplot(df, aes(x=groups, y=Itgam)) + geom_bar(stat="identity")
p <- p + theme_classic() + xlab("")
p.itgam <- p + theme(
    axis.text.x=element_text(angle=315, hjust=0, colour="black", size=12),
    axis.text.y=element_text(colour="black", size=12),
    axis.title=element_text(size=14, face="italic"),
    legend.text=element_text(size=14),
    plot.margin=margin(0.5, 0.5, 0, 0, "cm")
)
p <- ggplot(df, aes(x=groups, y=S100a6)) + geom_bar(stat="identity")
p <- p + theme_classic() + xlab("")
p.s100a6 <- p + theme(
    axis.text.x=element_text(angle=315, hjust=0, colour="black", size=12),
    axis.text.y=element_text(colour="black", size=12),
    axis.title=element_text(size=14, face="italic"),
    legend.text=element_text(size=14),
    plot.margin=margin(0.5, 0.5, 0, 0, "cm")
)
# width 4, height 3


# test an enrichment of Ahnak-like genes
# load the list of human/mouse gene symbols fetched from ensembl data base.
anno_genes <- read.csv(file.path(inputDir, "human_mouse_geneID.csv"))

# convert Ahnak-like mouse gene symbols to human orthologues
likeAhnak <- c()
for (gene in AhnakLikeMouseGenes)
{
    likeAhnak <- c(likeAhnak, as.character(anno_genes$Gene.name[anno_genes$Mouse.gene.name == gene]))
}
# remove duplicated names
likeAhnak <- likeAhnak[!duplicated(likeAhnak)]


# Save a list of the Ahnak/AHNAK-covariant genes.
pad <- rep(NA, length(AhnakLikeMouseGenes) - length(likeAhnak))
df_AhnakLike <- data.frame(
    mouse=AhnakLikeMouseGenes,
    human=c(likeAhnak, pad)
)
write.csv(df_AhnakLike, file.path(inputDir, 'immgen_AhnakLikeGenes.csv'))



# load 'geneClusters', annotation data frame of gene cluster
load(file.path(RDataDir, "hc_4comps_6clusters.RData"))
geneClusters <- hc_list[['ward.D2']]$geneClusters


# test enrichment of Ahnak-like genes
enrich.list <- sapply(
    unique(geneClusters$cluster),
    function(i){test.enrichment(
        geneClusters$symbol[geneClusters$cluster==i],
        TARGET.LIST=likeAhnak,  # 293 human homologues
        TARGET.SIZE=length(AhnakLikeMouseGenes),  # use 333 murine genes to calculate the background ratio of Ahnak-like genes
        BACKGROUND.SIZE=geneNum.immgen  # total mouse gene number in Immgen
    )},
    simplify=FALSE
)
enrich.df <- as.data.frame(do.call(rbind, enrich.list))

# make the columns numeric
for (col in colnames(enrich.df))
{
    enrich.df[, col] <- as.numeric(enrich.df[, col])
}

# plot -log10 p-value and Ahnak-like gene ratio (%)
p <- ggplot(enrich.df, aes(x=ratio*100, y=-log10(p.value))) + geom_point(aes(size=size))
p <- p + scale_size(limits=c(249, 1500), breaks=c(250, 500, 1000))
p <- p + ggrepel::geom_text_repel(aes(label=rownames(enrich.df), size=140), box.padding=0.3, force=10)

p <- p + xlab("Ahnak-like genes (%)") + ylab("-log10 p-value") + labs(size="Total gene count")
p <- p + theme_classic()
p <- p + theme.scatter
p.enrich <- p

if (saveImage)
{
    pngStore(p.enrich, file.path(imageDir, "immgen_firstCohort_AhnakLikeGeneEnrichment.png"), HEIGHT=3.5, WIDTH=5)
    timeStamp("Image saved!")
}




timeStamp("All done!")
