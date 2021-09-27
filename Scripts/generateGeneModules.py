import os
import pandas as pd
import argparse
import glob
import os.path
from scipy import stats

def parse_args():
    parser = argparse.ArgumentParser()

    requiredArgs = parser.add_argument_group("Required arguments")

    requiredArgs.add_argument(
        "-rf",
        "--result_folder",
        dest = "result_folder",
        nargs="+",
        required=True,
        help="Directory containing tsv files respresenting S-matrixes from ICA"
    )
    requiredArgs.add_argument(
        "-gaf",
        "--Additional_Annotation_File",
        dest = "additional_annotation_file",
        nargs="+",
        required=True,
        help="File containing additional annotations for genes, Columns:[GeneID,Genome,Annotation] ; may be renamed"
    )

    return parser.parse_args()


def main():
    #Method takes the most significant outlier of the distribution of S-weights, calculates the k2 representation of the remainder 
    #Repeats until threshold is reached
    #Sign of components can be swapped, abs
    args = parse_args()
    print(F'{args.result_folder[0]}/*_S.tsv')
    Additional_Annotation_File = pd.read_table(args.additional_annotation_file[0] , sep='\t', header = 0)
    resultDF = pd.DataFrame(columns=["Genome","Component","GeneID","Annotation"])
    
    for S_MatrixFile in glob.iglob(F'{args.result_folder[0]}/*_S.tsv'):
        print(S_MatrixFile)
        sDF = pd.read_table(S_MatrixFile , sep='\t', header = 0)
        print(sDF.columns)
        #For each component in S
        for component in sDF.columns[1:]:
            #sort by largest outlier first
            descendingGeneWeights = abs(sDF[component]).sort_values(ascending=False)
            K2, p = stats.normaltest(sDF.loc[:,component])
            print(K2)
            print(p)
            i = 0
            #'take' best gene/outlier and calculate the probablity of the remainder being a normal distribution (aka noise)
            while p < 1e-3:
                i += 1
                K2, p = stats.normaltest(sDF.loc[descendingGeneWeights.index[i:],component])
            #amount of while loops (i) gives the 'taken' genes
            selectedGenes = descendingGeneWeights.iloc[:i]
            for gene in selectedGenes.index:
                GeneID = sDF.loc[gene,'GeneID']
                AnnotationRow = Additional_Annotation_File.loc[Additional_Annotation_File['GeneID'] == GeneID]
                Annotation = AnnotationRow["Annotation"].to_string(index=False)
                Start = AnnotationRow["Start"].to_string(index=False)
                End = AnnotationRow["End"].to_string(index=False)
                Strand = AnnotationRow["Strand"].to_string(index=False)
                dataDict = {'Genome':[os.path.basename(S_MatrixFile)],'Component':[component],'GeneID':[GeneID],'Start':[Start],'End':[End],'Strand':[Strand],'Annotation':[Annotation]}
                toAppendDF = pd.DataFrame(data=dataDict)
                resultDF = resultDF.append(toAppendDF, ignore_index = True)
                
    resultDF.to_csv(F"results/8_ICA/GeneModules.tsv", sep='\t', index = False)

if __name__ == "__main__":
    main()
