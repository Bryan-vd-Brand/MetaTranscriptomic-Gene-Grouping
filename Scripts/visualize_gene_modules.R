#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(dplyr)
library(gggenes)
if (length(args)==0) {
  stop("No Arguments Supplied", call.=FALSE)
}

FeatureCountTSV <- args[1]
GeneModulesTSV <- args[2]
GroupedGeneModuleTSV <- args[3]
GeneAnnotationTSV <- args[4]

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

#TODO: group the gene module entries by genome and redo the componentID values to have a new set for each genome (simplifying the coloring)
#Use gggenes to visualize the grouped gene modules for all gene modules with count > 1
AllGeneModuleTable <- read.table(file = GroupedGeneModuleTSV, sep = '\t', header = TRUE)
GroupedGeneModuleTable <- AllGeneModuleTable %>% group_by(AllGeneModuleTable$Genome) %>% group_split()
GeneAnnotationTable <- read.table(file= GeneAnnotationTSV, sep = '\t', header = TRUE)
LongVisualizeableGenesTable <- data.frame(Genome = character(), GeneID = character(), Start = numeric(), End = numeric(), Strand = numeric(), ComponentID = numeric())
#First take the grouped gene modules from the table, combine them with annotation information in long format
for (GeneModuleTable in GroupedGeneModuleTable){
 ComponentID <- 0
 for (row in 1:nrow(GeneModuleTable)){
  Genome <- GeneModuleTable$Genome[row]
  Count <- GeneModuleTable$Count[row]
  #If the gene module occurs less then X times of the 1000 runs ignore it
  if (Count < 500){
   next
  }
  ComponentID <- ComponentID + 1
  moduleGenes <- unlist(strsplit(substr(GeneModuleTable$GeneIDs[row], 2, nchar(GeneModuleTable$GeneIDs[row])-1), split = ", ")) #remove [] and create list out of the string
  GeneAnnotationSubset <- filter(GeneAnnotationTable, Chr == Genome)
  #For each gene look up its start, end and strand and add it to the LongVisualizeableGenesTable
  for (gene in moduleGenes){
   GeneAnnotation <- filter(GeneAnnotationSubset, GeneID == gene)
   #Incase of multiple copy number variation add them all 
   for (row in 1:nrow(GeneAnnotation)){
    Start <- GeneAnnotation$Start[row]
    End <- GeneAnnotation$End[row]
    if (GeneAnnotation$Strand[row] == '+'){
     Strand <- 1
    } else {
     Strand <- 0
    }
    newRow <- data.frame(Genome = c(Genome), GeneID = c(gene), Start = c(Start), End = c(End), Strand = c(Strand), ComponentID = c(ComponentID))
    LongVisualizeableGenesTable <- rbind(LongVisualizeableGenesTable, newRow)
   }
  }
 }
}
 
 #Then ggplot + gggenes
 #Because of the same gene being in several components make a seperate plot for each genome, facet wrap over componentID
#GroupedGeneModuleTable <- AllGeneModuleTable %>% group_by(AllGeneModuleTable$Genome) %>% group_split()
GroupedLVSGT <- LongVisualizeableGenesTable %>% group_by(LongVisualizeableGenesTable$Genome) %>% group_split()
for (GenomesTable in GroupedLVSGT){
 genome <- GenomesTable$Genome[1]
 genePlot <- ggplot(GenomesTable, aes(xmin = Start, xmax = End, y = Genome, fill = as.factor(ComponentID), forward = Strand)) +
   geom_gene_arrow() +
   facet_wrap(~ ComponentID, ncol = 1) +
   scale_fill_brewer(palette = "Paired") +
   theme_genes()
 ggsave(sprintf("500_GeneModules_%s.png",genome),width=100,units="cm", plot = genePlot, limitsize = FALSE)
}





