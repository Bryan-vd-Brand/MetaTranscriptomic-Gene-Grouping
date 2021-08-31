#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
if (length(args)==0) {
  stop("No Arguments Supplied", call.=FALSE)
}

FeatureCountTSV <- args[1]


FeatureCountTable <- read.table(file = FeatureCountTSV, sep = '\t', header = TRUE)
rpkmHeatMapTotal <- ggplot(data = FeatureCountTable, mapping = aes(x=RunAccession,y=GeneID,fill=RPKM)) + geom_tile() + xlab(label="RPKM Gene Expression using Mapped read #") +
  scale_fill_gradientn(limits= c(0,100),colours = c("navyblue","darkmagenta","darkorange1")) +
  facet_wrap(~Genome, scales = "free")
rpkmHeatMapTotal
ggsave("FeatureCountGeneExpression.png",width=180,height=100,units="cm",limitsize = FALSE)