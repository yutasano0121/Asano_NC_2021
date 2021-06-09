timeStamp <- function(comment)
{
    ts <- paste(comment, Sys.time())
    print(ts)
}


pngStore <- function(plot_data, output_name, WIDTH, HEIGHT)
{
    png(output_name, height = HEIGHT, width = WIDTH, unit = "in", res = 600)
    print(plot_data)
    dev.off()
}


easyPlot <- function(p, shape.lab="", colour.lab="", fill.lab="", type="tsne", xlabel=NULL, title="")
{
    p <- p + theme_classic()
    p <- p + theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12, colour = "black"))
    p <- p + theme(legend.title = element_text(size = 14), legend.text = element_text(size = 14))
    if (type == "tsne"){
        p <- p + xlab("tSNE-1") + ylab("tSNE-2")
        if (shape.lab == "Tissue source"){p <- p + scale_shape_manual(values = c(16, 5))}
        if (colour.lab == "Patient"){
            colorPal <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(n = 7, name = "Paired"))
            p <- p + scale_color_manual(values=colorPal)
        }
        p <- p + labs(colour = colour.lab, shape = shape.lab)
    } else if (type == "volcano"){
        p <- p + xlab(xlabel) + ylab("- log10 FDR")
        p <- p + theme(legend.position="none")
    } else if (type == "violin"){
        p <- p + theme(axis.text = element_text(size = 14, colour = "black"))
        p <- p + labs(colour = colour.lab, fill = fill.lab)
    } else {print("reduction_method value error.")}


    if (!is.na(title)){
        p <- p + ggtitle(title)
        p <- p + theme(plot.title=element_text(hjust=0.5, size=18, face="bold"))
    }
    return(p)
}

# find genes expressed at a frequency higher than the threshold
find_expressedGenes <- function(DATA, FREQUENCY)  # data needs to be genes x cells
{
    gene.count <- rowSums(DATA != 0)
    detection.rate <- gene.count / dim(DATA)[2]
    expressed_genes <- detection.rate > FREQUENCY  # default 0.1
    return(expressed_genes)
}

# select high-variance genes in the expression matrix (genes x cells)
subset_genes <- function(data, geneNum=feature_num)
{
    gene_var <- apply(data, 1, var)
    data_ordered <- data[order(gene_var, decreasing=TRUE), ]
    data_ordered <- data_ordered[1:geneNum, ]
    return(data_ordered)
}


# find most highly expressed Ig isotype. pass it to apply(counts, 2, findIg)
findIg <- function(x){
    ensg <- HC[which.max(x)]
    return(anno.geneNames$gene_symbol[anno.geneNames$ENSG==ensg])
}


# Upper-quantile normalization (gene x cell)
quartile_norm <- function(DF)
{
    a <- apply(DF, 2, function(x){x / quantile(x[x != 0], 0.75)})
    return (a)
}


# run t-test on each gene (row of a count data frame). returns a list.
applyTtest2df <- function(DF, CATEGORY)  # catefory needs to be a list of 2 vectors
{
    applyTtest2row <- function(GENE)
    {
        t.result <- t.test(x=GENE[CATEGORY[[1]]], y=GENE[CATEGORY[[2]]], paired=FALSE)
        return(t.result)
    }

    t.result <- apply(DF, 1, applyTtest2row)

    return(t.result)
}


# make a summary data frame from a result of applyTtest2df
summaryTtest <- function(T.RESULT,
                        FDR.THRESHOLD=fdr.threshold, DIFF.THRESHOLD=diff.threshold)
{
    fetch.stats <- function(NAME)
    {
        t.res <- T.RESULT[[NAME]]
        p.value <- t.res$p.value
        mean.difference <- diff(as.numeric(t.res$estimate))  # remove names of the 'estimate' list

        return(c(p.value=p.value, mean.difference=mean.difference))
    }

    # make a data frame of p.values and difference of mean expression
    df.stats <- data.frame(do.call(rbind, lapply(names(T.RESULT), fetch.stats)))

    # calculate FDR
    df.stats$FDR <- p.adjust(df.stats$p.value, "BH")

    # make gene names the row names
    rownames(df.stats) <- names(T.RESULT)

    # add a column of significance threshold
    df.stats$threshold <- (df.stats$FDR <= FDR.THRESHOLD) &
                          (abs(df.stats$mean.difference) >= DIFF.THRESHOLD)

    return(df.stats)
}



# convert the rownames of a DEG data frame into ENTREZID. return a data frame.
makeRownamesENTREZ <- function(DF)
{
    ids <- bitr(rownames(DF), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    # fetch duplicated gene symbols
    duplicated.symbols <- ids$SYMBOL[duplicated(ids$SYMBOL)]
    print(paste("Gene symbols with multiple ENTREZID are removed:", paste(duplicated.symbols)))
    # remove them
    ids <- ids[!ids$SYMBOL %in% duplicated.symbols, ]
    DF.out <- DF[rownames(DF) %in% ids$SYMBOL, ]
    # change rownames to ENTREZID
    ids <- ids[order(ids$SYMBOL), ]
    DF.out <- DF.out[order(rownames(DF.out)), ]
    rownames(DF.out) <- ids$ENTREZID

    return(DF.out)
}


# run GO-enrichment test on the list of differentially expressed genes
test.GOenrichment <- function(GENES,
                              ONT='BP', P.ADJUST='BH',
                              P.CUTOFF=0.05, Q.CUTOFF=0.2,
                              MIN.SIZE=10, MAX.SIZE=500,
                              READABLE=TRUE)
{
    ego <- enrichGO(gene          = GENES,
                    OrgDb         = 'org.Hs.eg.db',
                    keyType       = 'ENTREZID',
                    ont           = ONT,
                    pAdjustMethod = P.ADJUST,
                    pvalueCutoff = P.CUTOFF,
                    qvalueCutoff  = Q.CUTOFF,
                    minGSSize     = MIN.SIZE,
                    maxGSSize     = MAX.SIZE,
                    readable      = READABLE)
    return(ego)
}

# run GO-enrichment test on the list of differentially expressed genes
test.KEGGenrichment <- function(GENES, P.ADJUST='BH',
                                MIN.SIZE=10, MAX.SIZE=500,
                                P.CUTOFF=0.05)
{
    ekegg <- enrichKEGG(gene          = GENES,
                      organism      = 'hsa',
                      pAdjustMethod = P.ADJUST,
                      pvalueCutoff  = P.CUTOFF,
                      minGSSize     = MIN.SIZE,
                      maxGSSize     = MAX.SIZE)
    return(ekegg)
}



# run go-GSEA (input needs to be a ranked gene lists)
test.GOgse <- function(rankedGENES,
                      ONT='BP', P.ADJUST='BH',
                      P.CUTOFF=0.05, PERM = 1000,
                      MIN.SIZE=100, MAX.SIZE=500)
{
    gsego <- gseGO(geneList    = rankedGENES,
                  OrgDb        = 'org.Hs.eg.db',
                  keyType      = 'ENTREZID',
                  ont          = ONT,
                  pAdjustMethod = P.ADJUST,
                  nPerm        = PERM,
                  minGSSize    = MIN.SIZE,
                  maxGSSize    = MAX.SIZE,
                  pvalueCutoff = P.CUTOFF)
    return(gsego)
}


# run KEGG-GSEA (input needs to be a ranked gene lists)
test.KEGGgse <- function(rankedGENES,
                      ONT='BP', P.ADJUST='BH',
                      P.CUTOFF=0.05, PERM = 1000,
                      MIN.SIZE=100, MAX.SIZE=500)
{
    gsekegg <- gseKEGG(gene          = rankedGENES,
                      organism      = 'hsa',
                      pAdjustMethod = P.ADJUST,
                      pvalueCutoff  = P.CUTOFF,
                      nPerm         = PERM,
                      minGSSize     = MIN.SIZE,
                      maxGSSize     = MAX.SIZE)
    return(gsekegg)
}


# calculate signal-to-noise ratio
stn <- function(VECTOR, CATEGORY)
{
    comparison <- levels(CATEGORY)
    sd.1 <- sd(VECTOR[CATEGORY==comparison[1]])
    sd.2 <- sd(VECTOR[CATEGORY==comparison[2]])
    mean.1 <- mean(VECTOR[CATEGORY==comparison[1]])
    mean.2 <- mean(VECTOR[CATEGORY==comparison[2]])

    stnRatio <- (mean.2 - mean.1) / (sd.1 + sd.2)
    return(stnRatio)
}

# calculate Baumgartner-Weiss-Schindler statistics
bws <- function(VECTOR, CATEGORY, METHOD="BWS")  # METHOD = 'BWS' or 'NEuhauser'
{
    comparison <- levels(CATEGORY)
    vec1 <- VECTOR[CATEGORY==comparison[1]]
    vec2 <- VECTOR[CATEGORY==comparison[2]]
    # one-sided BWS test
    out <- bws_test(vec1, vec2, method=METHOD)
    if ((METHOD == 'BWS') & (mean(vec2) - mean(vec1) < 0)){out$statistic <- -out$statistic}
    return(out$statistic)
}


# volcano plot
# take an edgeR result as an input
plot.volcano <- function(TAG)
{
    df <- as.data.frame(TAG)
    p <- ggplot(data = df, aes(x = mean.difference, y = -log10(FDR), colour = threshold))
    p <- p + geom_point(alpha = 0.4, size = 1.75)
    p <- p + theme(legend.position = "none")
    p <- p + xlab(paste0("log2 fold change (Kidney / Tonsil)"))
    p <- p + ylab("-log10 FDR")
    #p <- p + xlim(c(-10, 10))
    p <- p + scale_colour_manual(values = c("darkgrey", "red"))

    tmp <- df[df$threshold, ]  # significant genes
    tmp <- tmp[order(tmp$FDR), ]  # order genes by fdr
    tmp1 <- rownames(tmp[tmp$mean.difference > 0, ])[1:10]
    tmp2 <- rownames(tmp[tmp$mean.difference < 0, ])[1:10]
    tmp <- tmp[order(-abs(tmp$mean.difference)), ]  # order genes by fdr
    tmp3 <- rownames(tmp[tmp$mean.difference > 0, ])[1:10]
    tmp4 <- rownames(tmp[tmp$mean.difference < 0, ])[1:10]
    tmp <- unique(c(tmp1, tmp2, tmp3, tmp4))
    p <- p + ggrepel::geom_text_repel(
        data=TAG[tmp ,],
        aes(x=mean.difference, y=-log10(FDR), label=tmp),
        segment.color = "darkgray",
        colour='black',
        force=2
    )  # don't forget x and y aesthetics variables here!

    # further format labeling
    p <- p + theme_classic()
    p <- p + theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12, colour = "black"))
    p <- p + theme(legend.position="none")


    return(p)
}


# plot enriched GO/KEGG terms
# take an enrichGO/enrichKEGG result as an input
plot.enrichedTerms <- function(eRESULT, showCategory=NULL)
{
    df <- eRESULT@result
    df <- df[df$p.adjust < 0.05, ]  # subset significant terms
    df <- df[order(df$p.adjust), ]  # order the data frame by FDR
    if (!is.null(showCategory)){df <- df[1:showCategory, ]}
    p <- ggplot(data = df, aes(y = Count, x = reorder(Description, -p.adjust), fill = p.adjust))
    p <- p + geom_bar(stat = "identity") + coord_flip()
    p <- p + scale_fill_viridis(direction=-1, option="plasma")

    # format labeling
    p <- p + xlab("") + ylab("Gene count") + labs(fill = "FDR")
    p <- p + theme_classic()
    p <- p + theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12, colour = "black"))
    p <- p + theme(legend.title = element_text(size = 14), legend.text = element_text(size = 14))

    return(p)
}


# set color breaks for heatmap according to quantile distribution
quantile_breaks <- function(DATA, n = 100){
    breaks <- quantile(DATA, probs = seq(0, 1, length.out = n), na.rm=TRUE)
    breaks[!duplicated(breaks)]
}


# make a t-SNE plot
# change .LAB parameters as well when aesthetics parameters are changed
plot.tsne <- function(DF, ANNO,
                    X='V1', Y='V2',
                    SHAPE='Source', COLOUR='Patient',
                    SHAPE.LAB='Tissue source', COLOUR.LAB='Patient',
                    GENE=FALSE, GENE.DF=NULL)
{
    if (GENE){
        p <- ggplot(DF, aes(x=DF[, X], y=DF[, Y],
                        shape=ANNO[, SHAPE], colour=GENE.DF[, COLOUR]))
        p <- p + scale_colour_gradientn(colors=plasma(100, begin=0.1))
        COLOUR.LAB <- COLOUR
    } else {
        p <- ggplot(DF, aes(x=DF[, X], y=DF[, Y],
                        shape=ANNO[, SHAPE], colour=ANNO[, COLOUR]))
    }
    p <- p + geom_point()
    p <- p + xlab("tSNE-1") + ylab("tSNE-2")
    p <- p + labs(shape=SHAPE.LAB, colour=COLOUR.LAB)

    # format labeling
    if (SHAPE.LAB == "Tissue source")
    {
        p <- p + scale_shape_manual(values = c(16, 5))
    }

    if (COLOUR.LAB == "Patient")
    {
        colorPal <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(n = 7, name = "Paired"))
        p <- p + scale_color_manual(values=colorPal)
    }

    if (COLOUR.LAB == "Ig class")
    {
        colorPal <- c("IGHA1" = "gold1", "IGHA2" = "gold3", "IGHG1" = "brown1", "IGHG2" = "brown2", "IGHG3" = "brown3", "IGHG4" = "brown4", "IGHE" = "dodgerblue3", "IGHM" = "forestgreen", "IGHD" = "chartreuse3")
        labOrder <- c("IGHM", "IGHD", "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHE")
        p <- p + scale_color_manual(values=colorPal, breaks=labOrder)
    }

    p <- p + theme_classic()
    p <- p + theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12, colour = "black"))
    p <- p + theme(legend.title = element_text(size = 14), legend.text = element_text(size = 14))

    return(p)
}


plot.violin <- function(DF, ANNO, Y, X='Source', FILL=NULL, FILL.LAB="Add a label")
{
    DF <- as.data.frame(DF)
    if (is.null(FILL))
    {
        p <- ggplot(DF, aes(x=ANNO[, X], y=DF[, Y]))
    } else {
        p <- ggplot(DF, aes(x=ANNO[, X], y=DF[, Y], fill=ANNO[, FILL]))
        p <- p + labs(fill = FILL.LAB)
        p <- p + theme(legend.title = element_text(size = 14),
                        legend.text = element_text(size = 14))
    }

    p <- p + geom_violin(scale="width")
    p <- p + geom_boxplot(width=0.1, position=position_dodge(0.9))
    p <- p + xlab("") + ylab(Y)

    # format labeling
    p <- p + theme_classic()
    p <- p + theme(axis.title = element_text(size=14, colour='black'),
                    axis.text.x = element_text(size = 14, colour="black"),
                    axis.text.y = element_text(size = 12, colour="black"))

    return(p)
}

# make a bar plot with error bars
plot.bar <- function(DF, X, Y='mean', X.IN.DF=TRUE, ERRORBAR=TRUE, ERROR='se', Y.LAB="Score", WIDTH=NULL)
{
    if (X.IN.DF == FALSE)
    {
        # default. X-axis is determined by another data frame, list, or vector
        p <- ggplot(DF, aes(x=X, y=DF[, Y])) + geom_bar(stat="identity", width=WIDTH)
    } else {
        # if DF contains a column of a categorical variable to determine the x-axis
        p <- ggplot(DF, aes(x=DF[, X], y=DF[, Y])) + geom_bar(stat="identity", width=WIDTH)
    }

    # add error bars if needed
    if (ERRORBAR == TRUE)
    {
        p <- p + geom_errorbar(aes(ymin=DF[, Y] - DF[, ERROR], ymax=DF[, Y] + DF[, ERROR]),
                                    width=.2, position=position_dodge(.9))
    }

    # format labeling
    p <- p + xlab("") + ylab(Y.LAB)
    p <- p + theme_classic()
    p <- p + theme(axis.title = element_text(size=14, colour='black'),
                axis.text.x = element_text(size = 14, colour="black", angle=45, hjust=1),
                axis.text.y = element_text(size = 12, colour="black"))

    return(p)
}


plot.enrichNetwork <- function(enrichResult,
                                showCategory = 10,
                                kegg=FALSE,
                                saveCytoscapeObj=FALSE,
                                nameCytoscapeObj="temp")
{
    # make a data frame from an enrichResult object
    df <- as.data.frame(enrichResult)[1:showCategory, ]

    # make a list of genes for top enriched terms
    termsAndGenes <- strsplit(df$geneID, "/")
    names(termsAndGenes) <-  df$Description

    # make the list into a data frame for making a graph
    df_termsAndGenes <- lapply(names(termsAndGenes), function(name){
        data.frame(categoryID=rep(name, length(termsAndGenes[[name]])),
                   Gene=termsAndGenes[[name]])
        })
    df_termsAndGenes <- do.call('rbind', df_termsAndGenes)

    if (kegg == TRUE)
    {
        # change geneID from ENTREZID to SYMBOL
        library(org.Hs.eg.db)

        ids <- bitr(geneID=df_termsAndGenes$Gene,
                    fromType="ENTREZID", toType="SYMBOL",
                    OrgDb="org.Hs.eg.db",
                    drop=FALSE)
        df_termsAndGenes$Gene <- sapply(df_termsAndGenes$Gene, function(x){return(ids$SYMBOL[ids$ENTREZID==x])})
    }

    # convert the data frame into a graph
    g <- graph.data.frame(df_termsAndGenes, directed=FALSE)

    V(g)$size <- 1
    V(g)$size[1:showCategory] <- 4
    fdr <- df$p.adjust
    V(g)$color <- NA
    V(g)$color[1:showCategory] <- fdr
    #my.palette <- viridis(n=100, begin=0.2, end=0.9, direction=-1, option="plasma")
    my.palette <- viridis(n=100, direction=-1, option="plasma")


    if(saveCytoscapeObj==TRUE)
    {
        gg <- g
        V(gg)$category <- FALSE
        V(gg)$category[1:showCategory] <- TRUE
        fdr.scaled <- ((fdr - min(fdr)) / (max(fdr) - min(fdr))) * 99 + 1
        library(RCy3)
        cytoscapePing()
        createNetworkFromIgraph(gg, nameCytoscapeObj)
        timeStamp("igraph transfered to Cytoscape.")
    }

    V(g)$name[1:showCategory] <-  paste0("bold('", df$Description, "')")

    p <- ggraph(g, layout="kk") +
        geom_edge_link(alpha=.5, colour='darkgrey') +
    #    geom_node_point(aes_(color=~color), size=V(g)$size, shape=21) +
        geom_node_point(aes_(color=~color), size=V(g)$size) +
        scale_color_gradientn(name = "FDR", colors=my.palette, na.value = "dimgrey")
    coord <- layout_with_kk(g)
    p <- p + geom_text_repel(x=coord[, 1], y=coord[, 2], label=V(g)$name, segment.color="dimgrey", segment.linetype=3, parse=TRUE, box.padding=0.5, point.padding=0.4, max.iter=10000)
    p <- p + theme_void()

    return(p)
}


scaleColorLimit <- function(LIST.E.OBJECT)
{
    fetch.Plimits <- function(E.OBJECT)
    {
        df <- as.data.frame(E.OBJECT)
        df <- df[df$p.adjust < 0.05]  # significant ones
        maxP <- max(df$p.adjust)
        minP <- min(df$p.adjust)  # highest p-value
        return(c(max=maxP, min=minP))
    }

    Plimits <- lapply(LIST.E.OBJECT, fetch.Plimits)
    Plimits <- as.data.frame(do.call(rbind, Plimits))
    upperLimit <- max(Plimits$max)
    lowerLimit <- min(Plimits$min)
    limitList <- c(low=lowerLimit * (1 - 1E-6), high=upperLimit * (1 + 1E-6))

    if (FALSE)
    {
        # a script to make a color scale
        myPalette <- viridis(100, direction=-1, option="plasma")
        scale_fill_gradientn(colors=myPalette, limits=limitList, breaks=c())
    }

    return(limitList)
}
