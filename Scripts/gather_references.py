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

def determineGenomes(pileup_files):
    #this method takes a pileup file.
    #compares the result, takes the top couple of references based on some hardcoded rules 
    #Outputs the list of new references
    allResults = {}
    
    for file in pileup_files:
        result = []
        sample = os.path.basename(file).split(".")[0].split('_')[1]
        flagFile = F"results/4_ReferenceSelection/per_sample/{sample}/{sample}.flagstat"
        flagStatTable = pd.read_table(flagFile, sep='\t', header = None)
        numOfUniquelyMappedReads = flagStatTable[0][4]
        
        pileupTable = pd.read_table(file, sep='\t', header = 0)
        #Take the genomes having atleast 50% Covered_percent for this sample ; hcov Format = 0.3406
        for row in pileupTable.iterrows():
            Series = row[1]
            genomeName = Series.at['#ID']
            hCov = Series.at['Covered_percent']
            #if below 50% stop loop, next row
            if hCov < 50.0:
                continue
            else:
                result.append(F"{genomeName}:{hCov}")
                #count the number of passed samples for a genome
                if genomeName in allResults.keys():
                    allResults[genomeName] = allResults[genomeName] + 1
                else:
                    allResults[genomeName] = 1               
        with open('results/4_ReferenceSelection/ReferenceTable.tsv', 'a') as genomeTsv:
            if len(result) > 0:
                genomeTsv.write(F"{sample}\t{numOfUniquelyMappedReads}\t{result}\n")          
            else:
                print(F"No reference genomes for {sample}")
    with open('results/4_ReferenceSelection/ReferenceList.txt', 'w') as referenceList:
        for genome in allResults.keys():
            referenceList.write(F"{genome}\t{allResults[genome]} \n")

def main():

    args = parse_args()
    determineGenomes(args.input_files)

if __name__ == "__main__":
    main()
