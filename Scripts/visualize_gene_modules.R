#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(dplyr)
if (length(args)==0) {
  stop("No Arguments Supplied", call.=FALSE)
}

FeatureCountTSV <- args[1]
GeneModulesTSV <- args[2]

FeatureCountTable <- read.table(file = FeatureCountTSV, sep = '\t', header = TRUE)
GeneModuleTable <- read.table(file = GeneModulesTSV, sep = '\t', header = TRUE)
#RunAccession	Genome	GeneID	RPKM
FeatureCountSubsetDF <- data.frame(RunAccession = character(), Genome = character(), GeneID = character(), Component = numeric(), RPKM = double())
#take from long format featurecount table all rows with a geneID contained in the geneModuleTable
#Then generate a rpkmHeatMap from those seperated by component and genome (facet_grid(c ~ g)?)
for (row in 1:nrow(GeneModuleTable)){
  genome <- substr(GeneModuleTable$Genome[row],1,nchar(GeneModuleTable$Genome[row])-6)
  component <- GeneModuleTable$Component[row]
  geneID <- GeneModuleTable$GeneID[row]
  annotation <- GeneModuleTable$Annotation[row]
  #Select all RPKM entries for this genome/gene combination
  SelectedRows <- filter(FeatureCountTable, Genome == genome & GeneID == geneID)
  #Add the component number and save it to a new table
  SelectedRows <- mutate(SelectedRows, Component = component)
  FeatureCountSubsetDF <- rbind(FeatureCountSubsetDF, SelectedRows)
}

perGenomeDF <- split(FeatureCountSubsetDF, FeatureCountSubsetDF$Genome)

for (gDF in perGenomeDF){
rpkmHeatMapGroup <- ggplot(data = gDF, mapping = aes(x=RunAccession,y=GeneID,fill=RPKM)) + geom_tile() + xlab(label="RPKM Gene Expression of Gene Modules") +
  scale_fill_gradientn(limits= c(0,100),colours = c("navyblue","darkmagenta","darkorange1")) + theme(axis.text = element_text(size = 22)) +
  facet_wrap(~Component, scales = "free")
print(gDF$Genome[1])
FileName <- paste(gDF$Genome[1], "GeneModulesExpressionHeatMap.png", sep="_")
ggsave(FileName,width=180,height=100,units="cm",limitsize = FALSE)
}