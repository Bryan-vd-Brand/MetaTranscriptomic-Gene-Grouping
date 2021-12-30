#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
if (length(args)==0) {
  stop("No Arguments Supplied", call.=FALSE)
}

options(warn = 1)

GroupedGeneModuleTSV <- args[1]
#Genome	GeneModule	ICACount	Pval
AllGeneModuleTable <- read.table(file = GroupedGeneModuleTSV, sep = '\t', header = TRUE)

#visualize ICA stability (# of occurences in 1000 runs) vs Pval (relative clusteredness compared to random gene modules)
reliabilityPlot <- ggplot(AllGeneModuleTable, aes(x=ICACount, y=(Pval/1000), colour=Genome)) +
   geom_point() +
   xlab(label="Number of occurences in 1000 runs of ICA") +
   ylab(label="P-value")

reliabilityPlot + scale_color_brewer(palette = "Set1")    
ggsave("Stability.png", plot = reliabilityPlot, limitsize = FALSE)