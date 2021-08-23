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
for (pileUpFile in listOfPileUpFiles){
  split <- unlist(strsplit(pileUpFile, '_'))
  SampleName <- unlist(strsplit(split[4], '[.]'))[1]
  pileupTable <- read.table(file = pileUpFile, sep = '\t', header = FALSE)
  for (row in 1:nrow(pileupTable)){
    rowDF <- data.frame(RunAccession = c(SampleName), Genome=c(pileupTable[row,1]), Hcov=c(pileupTable[row,5]), stdDev = c(pileupTable[row,11]))
    PileUpDF <- rbind(PileUpDF, rowDF)
  }
}

hcovHeatMap <- ggplot(data = PileUpDF, mapping = aes(x=RunAccession,y=Genome,fill=Hcov)) + geom_tile() + xlab(label="Hcov of all sample-Genome combinations") +
  scale_fill_gradientn(limits= c(1,100),colours = c("navyblue","darkmagenta","darkorange1"))
ggsave("hcovHeatMap.png", width=180,height=100,units="cm", plot = hcovHeatMap, limitsize = FALSE)

stdDevHeatMap <- ggplot(data = PileUpDF, mapping = aes(x=RunAccession,y=Genome,fill=stdDev)) + geom_tile() + xlab(label="stdDev of all sample-Genome combinations") +
  scale_fill_gradientn(limits= c(1,100),colours = c("navyblue","darkmagenta","darkorange1"))
ggsave("stdDevHeatMap.png", width=180,height=100,units="cm", plot = stdDevHeatMap, limitsize = FALSE)