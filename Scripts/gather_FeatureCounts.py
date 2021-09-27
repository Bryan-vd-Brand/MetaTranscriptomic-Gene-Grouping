import os
import pandas as pd
import argparse
import glob
import os.path

def parse_args():
    parser = argparse.ArgumentParser()

    requiredArgs = parser.add_argument_group("Required arguments")

    requiredArgs.add_argument(
        "-rf",
        "--result_folder",
        dest = "result_folder",
        nargs="+",
        required=True,
        help=""
    )
    
    requiredArgs.add_argument(
        "-rl",
        "--SelectedReferences",
        dest = "reference_list",
        nargs="+",
        required=True,
        help=".tsv file containing references selected for in previous rule"
    )
    return parser.parse_args()



def main():
    #Runs over all per_sample/[Sample]/[Sample]_RPKM.tsv files and grabs the rows for which genomes were selected ; combines to one file
    args = parse_args()
    ResultDataFrame = pd.DataFrame(columns=["RunAccession","Genome","GeneID","RPKM"])
    references = pd.read_table(args.reference_list[0] , sep='\t', header = 0)
    references = references.sort_values(by=['SampleCount'], ascending=False, ignore_index=True)
    #Set of best 5 references across all samples

    for rpkmFile in glob.iglob(F'{args.result_folder[0]}*/*_RPKM.tsv'):
        basename = os.path.basename(rpkmFile)
        Sample = basename.split('_')[0]
        rpkmDF = pd.read_table(rpkmFile, sep='\t', header = 0)
        for x in range(0,5):
            ReferenceGenome = references['Genome'][x]
            #for each reference genome look in the pileupfile if the hcov is good enough
            
            PileUpFile = F"{os.path.dirname(rpkmFile)}/pileup_{Sample}.txt"
            HorizontalCoverage = 0
            if os.path.exists(PileUpFile):
                pileupTable =  pd.read_table(PileUpFile, sep='\t', header = 0)
                pileupRow = pileupTable[pileupTable['#ID'] == ReferenceGenome]
                if len(pileupRow) > 1:
                    print("Error pileupRow not 1 row")
                HorizontalCoverage = pileupRow['Covered_percent'].iloc[0]
        
            if(HorizontalCoverage < 50.0):
                #print(f"Hcov below 50 for {Sample} {ReferenceGenome} skipping ; {HorizontalCoverage}")
                continue
            #If so, take from the rpkmFile all entries with that genome and store them in the global output
            rpkmSubset = rpkmDF[rpkmDF['Genome'] == ReferenceGenome]
            #Because later pivot and no discernable benefit of CNV's condense
            condense_functions = {'RunAccession':'first','Genome':'first','GeneID':'first','RPKM':'sum'}
            rpkmSubset = rpkmSubset.groupby(rpkmSubset['GeneID']).aggregate(condense_functions)
            rpkmSubset
            ResultDataFrame = ResultDataFrame.append(rpkmSubset, ignore_index = True) 
            
    print(ResultDataFrame)
    ResultFile = "FeatureCount_table.tsv"
    ResultDataFrame.to_csv(ResultFile, sep='\t', index = False)
    

if __name__ == "__main__":
    main()
