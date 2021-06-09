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


# load functions
source(file.path(workDir, "firstCohort_functions.R"))


path_to_tag.list <- sapply(
    simplify=FALSE,
    c(
        'kidney_Unswitched-Switched',
        'tonsil_Unswitched-Switched',
        'switched_Tonsil-Kidney',
        'unswitched_Tonsil-Kidney'
    ),
    function(comparison){
        return(paste0(DEGDir, '/firstCohort_edgeR_', comparison, '.csv'))
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
        tag$logFC <- -tag$logFC
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
                paste0(imageDir, '/firstCohort_volcano_', name, '.png'),
                WIDTH=5, HEIGHT=4.5
            )
        }
    )
}
