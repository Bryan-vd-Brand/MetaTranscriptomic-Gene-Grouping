import os
import pandas as pd
import argparse
import time



def parse_args():
    parser = argparse.ArgumentParser()

    requiredArgs = parser.add_argument_group("Required arguments")
    
    requiredArgs.add_argument(
        "-rl",
        "--SelectedReferences",
        dest = "reference_list",
        nargs="+",
        required=True,
        help=".tsv file containing references selected for in previous rule"
    )
    requiredArgs.add_argument(
        "-gd",
        "--GenomesDirectory",
        dest = "genomes_directory",
        nargs="+",
        required=True,
        help="Full path to a directory containing .fasta files for each of the genomes used"
    )
    requiredArgs.add_argument(
        "-g",
        "--genome_files",
        dest = "genome_files",
        nargs="+",
        required=True,
        help="List of Fasta"
    )
    
    return parser.parse_args()

def main():

    args = parse_args()
    
    references = pd.read_table(args.reference_list[0] , sep='\t', header = None)
    references = references.sort_values(by=[1], ascending=False, ignore_index=True)
    print(references)
    #Set of best 5 references across all samples
    SelectedReferences = []
    for x in range(4):
        SelectedReferences.append(references[0][x])
        
    for refName in SelectedReferences:
        #Create Index for selected Reference genome
        for referenceFile in args.genome_files:
            if referenceFile.find(refName) > 0:
                newGenomeDirectory = F"{args.genomes_directory[0]}/{refName}"
                #if the directory does not exist the index for this genome has not been generated before. 
                if not os.path.isdir(newGenomeDirectory):
                    os.system(F"mkdir {newGenomeDirectory}")
                    os.system(F"STAR --runThreadN 5 --runMode genomeGenerate --genomeDir {newGenomeDirectory} --genomeFastaFiles {referenceFile} --genomeSAindexNbases 7")
                      
    os.system(F"touch ./results/5_SingleReferenceGeneExpression/finished_Index.touch")

if __name__ == "__main__":
    main()