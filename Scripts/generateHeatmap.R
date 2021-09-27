#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggdendro)

if (length(args)==0) {
  stop("No Arguments Supplied", call.=FALSE)
}

FeatureCountTSV <- args[1]

#RunAccession	Genome	GeneID	RPKM
FeatureCountTable <- read.table(file = FeatureCountTSV, sep = '\t', header = TRUE)
rpkmHeatMapTotal <- ggplot(data = FeatureCountTable, mapping = aes(x=RunAccession,y=GeneID,fill=RPKM)) + geom_tile() + xlab(label="RPKM Gene Expression using Mapped read #") +
  scale_fill_gradientn(limits= c(0,100),colours = c("navyblue","darkmagenta","darkorange1")) +
  facet_wrap(~Genome, scales = "free")
rpkmHeatMapTotal
ggsave("FeatureCountGeneExpression.png",width=180,height=100,units="cm",limitsize = FALSE)

perGenomeDF <- split(FeatureCountTable, FeatureCountTable$Genome)


for (gDF in perGenomeDF){
#Cluster
gDF_Wide <- spread(gDF, GeneID, RPKM)
gDF_Wide[is.na(gDF_Wide)] <- 0
gDF_Wide[,-c(1,2)] <- scale(gDF_Wide[,-c(1,2)])
rownames(gDF_Wide) <- gDF_Wide$RunAccession
gDF_Wide <- gDF_Wide[,-c(1,2)]
gDF_Wide <- as.data.frame(t(gDF_Wide))
d <- dist(gDF_Wide, method = "euclidean")
hc <- hclust(d, method = "complete")
dendroD <- dendro_data(hc, type = "rectangle")
#Print Plot
rpkmHeatMapGroupClustered <- ggplot(segment(dendroD)) + geom_segment(aes(x=x,y=y,xend = xend, yend = yend)) +
geom_text(data = label(dendroD),aes(x = x, y = y, label = label, hjust = 0), size = 3)  + coord_flip() + scale_y_reverse(expand = c(0.2,0))
FileName <- paste(gDF$Genome[1], "Dendrogram.png", sep="_")
ggsave(FileName,width=60,height=120,units="cm",limitsize = FALSE)
}