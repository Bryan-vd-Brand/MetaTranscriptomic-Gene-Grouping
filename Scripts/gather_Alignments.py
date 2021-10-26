import os
import pandas as pd
import argparse
import glob
import os.path

def parse_args():
    parser = argparse.ArgumentParser()

    requiredArgs = parser.add_argument_group("Required arguments")
    
    return parser.parse_args()



#select set of prefered references and do index for them ; run over all samples, if ref in n50 selection do analysis ; create refname_samplename_featurecount result
def main():
    args = parse_args()
    ResultDataFrame = pd.DataFrame(columns=["RunAccession","UniqueCount","MMCount","PercentageMM"])
    #RunAccession	Genome	GeneID	Count
    for AlignmentFile in glob.iglob('../results/4_ReferenceSelection/per_sample/*/*_Alignments.tsv'):
        print(F"Processing {AlignmentFile}")
        AlignDF = pd.read_table(AlignmentFile, sep='\t', header=None)
        RunAccession = os.path.basename(AlignmentFile).split('_')[0]
        UniqueCount = 0
        MMCount = 0
        for row in AlignDF.itertuples():
            #If gene is very large save it to a file for manual inspection
            AlignmentNumber = row[1]
            AlignmentCount = row[2]
            if AlignmentNumber == 1:
                UniqueCount = UniqueCount + AlignmentCount
                continue
            MMCount = MMCount + (AlignmentCount / AlignmentNumber)
        rowDF = pd.DataFrame(data={'RunAccession':[RunAccession],'UniqueCount':[UniqueCount],'MMCount':[MMCount],'PercentageMM':[MMCount/(UniqueCount+MMCount)]})
        ResultDataFrame = ResultDataFrame.append(rowDF, ignore_index = True)             
    ResultDataFrame.to_csv("UniqueVsMM.tsv", sep='\t', index = False)

if __name__ == "__main__":
    main()