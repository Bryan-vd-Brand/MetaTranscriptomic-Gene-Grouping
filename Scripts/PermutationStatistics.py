import os
import pandas as pd
import argparse
import numpy as np
import os.path
import random
from statistics import mean
import multiprocessing

def parse_args():
    parser = argparse.ArgumentParser()

    requiredArgs = parser.add_argument_group("Required arguments")

    requiredArgs.add_argument(
        "-gaf",
        "--Additional_Annotation_File",
        dest = "additional_annotation_file",
        nargs="+",
        required=True,
        help="File containing additional annotations for genes, Columns:[GeneID,Genome,Annotation]"
    )
    
    requiredArgs.add_argument(
        "-ggm",
        "--Grouped_Gene_Modules",
        dest = "grouped_gene_modules",
        nargs="+",
        required=True,
        help="File containing gene modules and their occurance counts Columns:[Genome	GeneIDs	Count]"
    )
    
    return parser.parse_args()

#This script generates a distribution of distances between randomly generated gene modules selected from all genes of a genome
#Then compares the distance of a real/found gene module of the previous step to generate an indication of signifigance
#This is based on the core assumption that a genes that are co-expressed tend to be together on the genome
#Distance is defined as the mean of all distances between all genes in the module (Genes^2 distances)
AnnotationDictionary = {}

def GenerateRandomGeneModule(size, Genome):
    randoms = []
    numOfGenes = AnnotationDictionary[Genome].shape[0]
    randoms = random.sample(range(0,numOfGenes-1), size)
    #Do a N * N comparison of randoms/Genes list the resultant distances
    AllDistances = []
    for random1 in randoms:
        for random2 in randoms:          
            Gene1 = AnnotationDictionary[Genome].iloc[[random1]]
            Gene2 = AnnotationDictionary[Genome].iloc[[random2]]
            #if same skip
            if Gene1.iloc[0]["GeneID"] == Gene2.iloc[0]["GeneID"]:
               continue
            Distance = -1
            #Check if the random1 gene is before or after random2 gene
            if Gene1.iloc[0]["Start"] < Gene2.iloc[0]["Start"]:
                Distance = Gene2.iloc[0]["Start"] - Gene1.iloc[0]["End"]
            else:
                Distance = Gene1.iloc[0]["Start"] - Gene2.iloc[0]["End"]               
            if Distance == -1:
                print("Unexpected -1 Distance")
                continue
            AllDistances.append(Distance)        
    return mean(AllDistances)

def Do1000RandomGeneModules(List):
    Genome = List[0]
    size = List[1]
    MeanDistances = []
    print(F"Starting {Genome} and size {size}")
    for x in range(1,1001):
        distance = GenerateRandomGeneModule(size,Genome)
        MeanDistances.append(distance)
    return MeanDistances

def main():
    args = parse_args()
    Annotation_Table = pd.read_table(args.additional_annotation_file[0] , sep='\t', header = 0)
    Grouped_Gene_Module_Table = pd.read_table(args.grouped_gene_modules[0], sep = '\t', header = 0)
    
    #Split the annotations into a dict with key genome value annotations
    Grouped_AT = [pd.DataFrame(y) for x,y in Annotation_Table.groupby('Chr', as_index=False)]
    #split annotation file for each genome and create dictionary containing dataframe for each genome_files
    #Add a new column to the table for storing weighted read counts
    for geneTable in Grouped_AT:
        AnnotationDictionary[geneTable.iloc[0]['Chr']] = geneTable
    
    #For found gene modules, take the genomes and the sizes of their list, if found over 500 times
    Grouped_GM = [pd.DataFrame(y) for x,y in Grouped_Gene_Module_Table.groupby('Genome', as_index=False)]
    #Dictionary of Genomes containing a list of GeneModules (as lists)
    GeneModuleSizes = {}
    Found_GeneModules = {}
    for geneModuleSet in Grouped_GM:
        for row in geneModuleSet.itertuples():
            Index = row[0]
            Genome = row[1]
            #convert back to list
            GeneIDs = row[2].strip('][').replace("'","").split(", ")
            Count = row[3]
            #Skip extreme size
            if len(GeneIDs) > 75:
                continue
            #Skip if low occurance
            if Count < 500:
                continue
            #Cannot calculate distance for gene modules of size 1, skip
            if len(GeneIDs) == 1:
                continue
            #Keep track of the GeneModule sizes/Genome
            if Genome not in GeneModuleSizes.keys():
                GeneModuleSizes[Genome] = []
            if len(GeneIDs) not in GeneModuleSizes[Genome]:
                GeneModuleSizes[Genome].append(len(GeneIDs))
            #Collect [GeneModules,Count]/Genome
            if Genome not in Found_GeneModules.keys():
                Found_GeneModules[Genome] = []
                Found_GeneModules[Genome].append([GeneIDs,Count])
                continue
            Found_GeneModules[Genome].append([GeneIDs,Count])
            
    #Generate the random genemodules, calculate their internal distance and save the results as the distribution to compare against
    #Parralize
    pool = multiprocessing.Pool(processes = 50)
    RandomGeneModuleDistances = {}
    inputs = []     
    for Genome in GeneModuleSizes.keys():
        #Dict/Genome containing list of 1000 random distances for each gene module size
        RandomGeneModuleDistances[Genome] = {}
        for size in GeneModuleSizes[Genome]:
            if size == 1:
               continue
            inputs.append([Genome,size])
    print(F"{len(inputs)} jobs")
    outputs = pool.map(Do1000RandomGeneModules, inputs)
    for i in range(0,len(outputs)):
        Genome = inputs[i][0]
        size = inputs[i][1]
        RandomGeneModuleDistances[Genome][size] = outputs[i] 
    print(RandomGeneModuleDistances)
    
    #Now use the RandomGeneModuleDistances and compare to found genemodule distances and save p value
    ResultDataFrame = pd.DataFrame(columns=["Genome","GeneModule","ICACount","Pval"])
    
    for genome in Found_GeneModules.keys():
        for geneModule, Count in Found_GeneModules[genome]:
            print(F"Processing {genome} ; {geneModule} ; {Count}")
            Annotation = AnnotationDictionary[genome]
            #calc distance
            distances = []
            for geneID1 in geneModule:
                for geneID2 in geneModule:
                    #Skip same
                    if geneID1 == geneID2:
                        continue
                    Gene1 = Annotation.loc[Annotation['GeneID'] == geneID1]
                    Gene2 = Annotation.loc[Annotation['GeneID'] == geneID2]
                    Found_Distance = -1
                    #Check if gene1 is before or after gene2
                    if Gene1.iloc[0]["Start"] < Gene2.iloc[0]["Start"]:
                        Found_Distance = Gene2.iloc[0]["Start"] - Gene1.iloc[0]["End"]
                    else:
                        Found_Distance = Gene1.iloc[0]["Start"] - Gene2.iloc[0]["End"]               
                    if Found_Distance == -1:
                        print("Unexpected -1 Found_Distance")
                        continue
                    distances.append(Found_Distance)
            Mean = mean(distances)
            Pval = 0
            thisGenome = RandomGeneModuleDistances[genome]
            for random in thisGenome[len(geneModule)]:
                if random <= Mean:
                    Pval = Pval + 1
            print(F"Found Pval of {Pval}")
            rowRDF = pd.DataFrame(data={'Genome':[genome],'GeneModule':[geneModule],'ICACount':[Count],'Pval':[Pval]})
            ResultDataFrame = ResultDataFrame.append(rowRDF, ignore_index = True)
    ResultDataFrame.to_csv("results/9_PermutationStatistics/GeneModulePval.tsv", sep='\t', index = False)
                        
            
            
            
    
if __name__ == "__main__":
    main()
