# packages
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(ggrepel)

saveImage <- TRUE

# directory names
workDir <- getwd()  # Or enter a path to your directory of choice.
DEGDir <- file.path(workDir, "DEtest_out")
RDataDir <- file.path(workDir, "RData")
imageDir <- file.path(workDir, "image/scatter")

comparisons <- c(
    "TonsilS-KidneyS",
    "TonsilScluster-KidneyS",
    "DSAneg-DSApos"
)

# load functions
source(file.path(workDir, "integration_functions.R"))

path_to_tag.list <- sapply(
    simplify=FALSE,
    comparisons,
    function(comparison){
        return(paste0(DEGDir, '/integration_tTest_', comparison, '.csv'))
    }
)

# Load the data.
tag.list <- sapply(
    simplify=FALSE,
    names(path_to_tag.list),
    function(name){
        tag <- read.csv(path_to_tag.list[[name]], row.names=1)

        return(tag)
    }
)

# plot volcano
volcano.list <- sapply(
    simplify=FALSE,
    names(tag.list),
    function(name){
        tag <- tag.list[[name]]
        #tag$logFC <- -tag$logFC
        p.volcano <- plot.volcano(tag)

        return(p.volcano)
    }
)

if (saveImage){
    lapply(
        names(volcano.list),
        function(name){
            p.volcano <- volcano.list[[name]]
            pngStore(
                p.volcano,
                paste0(imageDir, '/integration_volcano_', name, '.png'),
                WIDTH=5, HEIGHT=4.5
            )
        }
    )
}
