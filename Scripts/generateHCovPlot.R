#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
if (length(args)==0) {
  stop("No Arguments Supplied", call.=FALSE)
}

listOfPileUpFiles <- Sys.glob(file.path(args[1], "*", "pileup_*.txt"))
#listOfPileUpFiles contains all the pileupFiles, format = /[SAMPLENAME]/pileup_[SAMPLENAME].txt

PileUpDF <- data.frame(RunAccession = character(), Genome = character(), Hcov = numeric(), stdDev = numeric())
#create table out of the results, X = hcov%, Y=BinOf% (histogram), Fill = REFNAME
ReferenceListDF <- read.table(file = args[2], sep = "\t", header = 1)
print(head(ReferenceListDF))
ReferenceListDF <- ReferenceListDF[order(ReferenceListDF$SampleCount, decreasing = TRUE),]
#reset indexs 
rownames(ReferenceListDF) <- NULL
print(head(ReferenceListDF))

for (pileUpFile in listOfPileUpFiles){
  filename <- basename(pileUpFile)
  SampleName <- unlist(strsplit(filename, '_'))[2]
  pileupTable <- read.table(file = pileUpFile, sep = "\t", header = FALSE)
  #grab hcov of the first 5 genomes
  for (i in c(1,2,3,4,5)){
  RefName <- ReferenceListDF$Genome[i]
  OneRowPileUp <- pileupTable[which(pileupTable[,1] == RefName),]
  if(OneRowPileUp[1,5] > 50){
  rowDF <- data.frame(RunAccession = c(SampleName), Genome=c(RefName), Hcov=c(OneRowPileUp[1,5]), stdDev = c(OneRowPileUp[1,11]))
  PileUpDF <- rbind(PileUpDF, rowDF)
  }
  }
}

hcovPlot <- ggplot(PileUpDF, aes(x=Hcov, color = Genome)) +
  geom_histogram(fill="white", alpha = 0.5, position="identity", binwidth=5)
ggsave("HcovHistogram.png", plot = hcovPlot)

hcovPlot_Density <- ggplot(PileUpDF, aes(x=Hcov, color = Genome)) +
  geom_density() + xlab(label="Density Plot of Hcov %")
ggsave("HcovDensity.png", plot = hcovPlot_Density)

stdDevPlot <- ggplot(PileUpDF, aes(x=stdDev, color = Genome)) +
  geom_histogram(fill="white", alpha = 0.5, position="identity", binwidth=0.5)
ggsave("stdDevHistogram.png", plot = stdDevPlot)