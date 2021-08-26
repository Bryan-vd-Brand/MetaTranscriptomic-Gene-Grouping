#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(tidyr)
if (length(args)==0) {
  stop("No Arguments Supplied", call.=FALSE)
}

ListOfSMatrixFiles <- Sys.glob(file.path(args[1],"*_S.tsv"))

#create table out of the results, X = hcov%, Y=BinOf% (histogram), Fill = REFNAME
for (SmatrixFile in ListOfSMatrixFiles){
  fileName <- basename(SmatrixFile)
  STable <- read.table(file = SmatrixFile, sep = '\t', header = TRUE)
  #convert to long format to plot all components as a different color
  long_STable <- pivot_longer(STable,!GeneID, names_to="Component", values_to="Weight")
  #make plot and save
  STableDistro <- ggplot(data = long_STable, mapping = aes(x=Weight,colour=Component)) + geom_density() + xlab(label=sprintf("Density Plot of S-weights for Genes of %s",fileName))
  ggsave(sprintf("%s.png",fileName), plot = STableDistro, path="./results/8_ICA/", limitsize = FALSE)
}


