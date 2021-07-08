import os
import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser()

    requiredArgs = parser.add_argument_group("Required arguments")

    requiredArgs.add_argument(
        "-i",
        "--input_files",
        dest = "input_files",
        nargs="+",
        required=True,
        help=""
    )
    return parser.parse_args()

def determineGenomes(idxstat_files):
    #this method takes a idxstat file.
    #compares the result, takes the top couple of references based on some hardcoded rules 
    #Outputs the list of new references
    allResults = {}
    
    for file in idxstat_files:
        result = []
        readCounter = 0
        sample = os.path.basename(file).split(".")[0]
        flagFile = F"results/4_ReferenceSelection/per_sample/{sample}.flagstat"
        flagStatTable = pd.read_table(flagFile, sep='\t', header = None)
        numOfUniquelyMappedReads = flagStatTable[0][4]
        
        idxStatsTable = pd.read_table(file, sep='\t', header = None)
        #Columns do not have names by default; 0 = reference name; 1=sequence length; 2=#mappedreads; 3=#unmappedreads
        idxStatsTable = idxStatsTable.sort_values(by=[2], ascending=False)
        #Take the best genomes of the file until atleast 50% of the reads are covered by those genomes
        for index, row in idxStatsTable.iterrows():
            genomeName = row[0]
            ReadCount = row[2]
            #if the reference has not enough mapped reads skip reference (likely till EOF)
            if ReadCount < 0.01 * int(numOfUniquelyMappedReads):
                continue
            
            #if we passed 50%, stop loop iteration
            if readCounter > 0.5 * int(numOfUniquelyMappedReads):
                break
            else:
                result.append(F"{genomeName}:{ReadCount}")
                
                if genomeName in allResults.keys():
                    allResults[genomeName] = allResults[genomeName] + 1
                else:
                    allResults[genomeName] = 1
                
                readCounter = ReadCount + readCounter
        with open('results/4_ReferenceSelection/ReferenceTable.tsv', 'a') as genomeTsv:
            if len(result) > 0:
                genomeTsv.write(F"{sample}\t{numOfUniquelyMappedReads}\t{result}\n")          
            else:
                #skip ; genomeTsv.write(F"{sample}\t{numOfUniquelyMappedReads}\tFAILED\n")
                print(F"Failed to select reference genomes for {sample}")
    with open('results/4_ReferenceSelection/ReferenceList.txt', 'w') as referenceList:
        for genome in allResults.keys():
            referenceList.write(F"{genome}\t{allResults[genome]} \n")

def main():

    args = parse_args()
    determineGenomes(args.input_files)

if __name__ == "__main__":
    main()
