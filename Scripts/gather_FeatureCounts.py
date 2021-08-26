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


ResultDataFrame = pd.DataFrame(columns=["RunAccession","Genome","GeneID","RPKM"])

def Read_Process_One_FeatureCount(ReferenceID, SampleRA, featureCountFile):
    SummaryFile = F"{featureCountFile}.summary"
    PileUpFile = F"{os.path.dirname(featureCountFile)}/pileup_{ReferenceID}_{SampleRA}.txt"
    NumCrAssPhageReads = -1
    DataTable = None
    HorizontalCoverage = 0
    
    print(PileUpFile)
    if os.path.exists(PileUpFile):
        pileupTable =  pd.read_table(PileUpFile, sep='\t', header = 0)
        HorizontalCoverage = pileupTable.at[0,'Covered_percent']
        
    if(HorizontalCoverage < 50.0):
        print(f"Hcov below 50 for {SampleRA} {ReferenceID} skipping ; {HorizontalCoverage}")
        return

    if os.path.exists(SummaryFile):
        with open(SummaryFile, 'r') as summFile:
            summFile.readline() # status line
            NumCrAssPhageReads = int(summFile.readline().split('\t')[1].rstrip("\n")) #Total Assigned Read Count
    else:
        print(F"ERROR: missing summary file for {SampleRA}; skipping!")
        return
    
    if os.path.exists(featureCountFile):
        DataTable = pd.read_table(featureCountFile, sep='\t', header = 1)
    else:
        print(F"ERROR: missing data file for {SampleRA}; skipping!")
        return
        
    #Example Columns for featureCountFile/Table = Geneid	Chr	Start	End	Strand	Length	ERR1356705_mapped.bam   
    GeneHitColumnName = F"{SampleRA}_mapped.bam"
    
    DataTable = DataTable.sort_values(by=GeneHitColumnName, ascending=False) #sort by gene hits
    DataTable = DataTable[DataTable[F'{GeneHitColumnName}'] > 1] #subset table, only take rows with atleast 1 hit for gene.
    DataTable = DataTable.rename(columns = {F'{GeneHitColumnName}':'GeneHitCount'}) #rename the varying sorted{SampleRA}.bam column name to static GeneHitCount

    oneRunDF = pd.DataFrame(columns=["RunAccession","Genome","GeneID","RPKM"])

    for row in DataTable.itertuples(index=False):
        GeneID = F"{row.Geneid.split('_')[0]}_{row.Geneid.split('_')[1]}"
        GeneLength = row.Length
        GeneHitCount = row.GeneHitCount
        RPKM = ((GeneHitCount/(NumCrAssPhageReads / 1000000))/GeneLength)
        dataDict = {'RunAccession':[SampleRA],'Genome':[ReferenceID],'GeneID':[GeneID],'RPKM':[RPKM]}
        toAppendDF = pd.DataFrame(data=dataDict)
        oneRunDF = oneRunDF.append(toAppendDF, ignore_index = True)
        
        
    #oneRunDF now contains a dataframe with the first column a run accession and the other columns having GeneID numbers as header and RPKM counts as row value (1row) 
    global ResultDataFrame
    ResultDataFrame = ResultDataFrame.append(oneRunDF, ignore_index = True) 

def main():
    #Runs over all per_sample/[Sample]/[Reference]_[Sample]_featurecount files and calculates the rpkm value and creates a table 
    args = parse_args()
    print(F'{args.result_folder[0]}*/*_featurecount')
    
    
    references = pd.read_table(args.reference_list[0] , sep='\t', header = None)
    references = references.sort_values(by=[1], ascending=False, ignore_index=True)
    #Set of best 5 references across all samples

    for featureCountFile in glob.iglob(F'{args.result_folder[0]}*/*_featurecount.summary'):
        featureCountTable = featureCountFile.rsplit('.',1)[0]
        basename = os.path.basename(featureCountTable)
        split = basename.split('_')
        Sample = split[len(split)-2]
        for x in range(5):
            if references[0][x] not in featureCountTable: #If the reference is not correct one for this specific globbed file, skip to the next one
                continue
            Read_Process_One_FeatureCount(references[0][x], Sample, featureCountTable)
    
    print(ResultDataFrame)
    ResultFile = "FeatureCount_table.tsv"
    ResultDataFrame.to_csv(ResultFile, sep='\t', index = False)
    

if __name__ == "__main__":
    main()
