import os
import pandas as pd
import argparse
import glob


def parse_args():
    parser = argparse.ArgumentParser()

    requiredArgs = parser.add_argument_group("Required arguments")
    
    return parser.parse_args()



#select set of prefered references and do index for them ; run over all samples, if ref in n50 selection do analysis ; create refname_samplename_featurecount result
def main():
    args = parse_args()
    ResultDataFrame = pd.DataFrame(columns=["RunAccession","Genome","GeneID","Count","Percentage"])
    #RunAccession	Genome	GeneID	Count
    for AnnotationFile in glob.iglob('../results/4_ReferenceSelection/per_sample/*/*_Annotation.tsv'):
        print(F"Processing {AnnotationFile}")
        AnnotationDF = pd.read_table(AnnotationFile, sep='\t', header=0)
        SumCount = AnnotationDF['Count'].sum()
        for row in AnnotationDF.itertuples():
            #If gene is very large save it to a file for manual inspection
            if (row.Count > (0.15  * SumCount)) & (row.Count > 100):
                print(F"Found Gene")
                rowDF = pd.DataFrame(data={'RunAccession':[row.RunAccession],'Genome':[row.Genome],'GeneID':[row.GeneID],'Count':[row.Count],'Percentage':[(row.Count/SumCount)]})
                ResultDataFrame = ResultDataFrame.append(rowDF, ignore_index = True)             
    ResultDataFrame.to_csv("FindOutliers15.tsv", sep='\t', index = False)

if __name__ == "__main__":
    main()