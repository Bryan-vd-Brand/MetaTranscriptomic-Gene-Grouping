#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
if (length(args)==0) {
  stop("No Arguments Supplied", call.=FALSE)
}

listOfPileUpFiles <- Sys.glob(file.path(args[1], "*", "pileup_*_*.txt"))
#listOfPileUpFiles contains all the pileupFiles, format = /[SAMPLENAME]/pileup_[REFNAME]_[SAMPLENAME].txt

PileUpDF <- data.frame(RunAccession = character(), Genome = character(), Hcov = numeric(), stdDev = numeric())
#create table out of the results, X = hcov%, Y=BinOf% (histogram), Fill = REFNAME
for (pileUpFile in listOfPileUpFiles){
  split <- unlist(strsplit(pileUpFile, '_'))
  RefName <- split[4]
  SampleName <- unlist(strsplit(split[5], '[.]'))[1]
  pileupTable <- read.table(file = pileUpFile, sep = '\t', header = FALSE)
  rowDF <- data.frame(RunAccession = c(SampleName), Genome=c(RefName), Hcov=c(pileupTable[1,5]), stdDev = c(pileupTable[1,11]))
  PileUpDF <- rbind(PileUpDF, rowDF)
}

hcovPlot <- ggplot(PileUpDF, aes(x=Hcov, color = Genome)) +
  geom_histogram(fill="white", alpha = 0.5, position="identity", binwidth=0.5)
ggsave("HcovHistogram.png", plot = hcovPlot)

stdDevPlot <- ggplot(PileUpDF, aes(x=stdDev, color = Genome)) +
  geom_histogram(fill="white", alpha = 0.5, position="identity", binwidth=0.5)
ggsave("stdDevHistogram.png", plot = stdDevPlot)