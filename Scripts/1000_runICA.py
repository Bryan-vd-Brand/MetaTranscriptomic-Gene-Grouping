import os
import pandas as pd
import argparse
from sklearn.decomposition import FastICA,PCA
import numpy as np
import glob
import os.path
from scipy import stats

def parse_args():
    parser = argparse.ArgumentParser()

    requiredArgs = parser.add_argument_group("Required arguments")
    
    requiredArgs.add_argument(
        "-ft",
        "--FeatureCount_Table",
        dest = "FeatureCount_table",
        nargs="+",
        required=True,
        help=".tsv file containing long format expression counts for the selected references and their genes"
    )
    
    requiredArgs.add_argument(
        "-PCAvar",
        "--PCAvariance",
        dest = "PCA_Variance",
        nargs="+",
        required=True,
        help="A config variable representing the % of variance to be represented by the dimensions of the PCA"
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


def OneICA(PCA_Variance, exprDFs, Additional_Annotation_File):
    resultDF = pd.DataFrame(columns=["Genome","Component","GeneID","Annotation"])
    for exprDF in exprDFs:
        Genome = exprDF['Genome'].iloc[0]
        #run for each genome seperately, avoiding finding gene presence/notpresence patterns
        wideDF = exprDF.pivot(index="RunAccession", columns="GeneID", values="RPKM").fillna(0) #This is the data for the X matrix in ICA
        #Now the data needs to be centered and scaled, handeled by the whiten=True in PCA & ICA respectively
        #rename to X for consistency with ICA paper
        X = wideDF
        #Use PCA to reduce the dimensions of the ICA algorithm, whiten to scale
        pca = PCA(whiten=True).fit(X)
        #explained_variance_ Equal to n_components largest eigenvalues of the covariance matrix of X ; ndarray      
        cumulative_variance = np.cumsum(pca.explained_variance_ratio_)
        number_of_components = np.where(cumulative_variance > PCA_Variance)[0][0] + 1

        #Run ICA
        X_t = X.transpose()
        ica = FastICA(whiten=True,max_iter=1000,tol=1e-6,n_components = number_of_components)
        S = pd.DataFrame(ica.fit_transform(X_t),index=X_t.index)
        A = pd.DataFrame(ica.mixing_,index=X_t.columns)
        
        #generate gene module
        #For each component in S
        for component in S.columns[1:]:
            #sort by largest outlier first
            descendingGeneWeights = abs(S[component]).sort_values(ascending=False)
            K2, p = stats.normaltest(S.loc[:,component])
            i = 0
            #'take' best gene/outlier and calculate the probablity of the remainder being a normal distribution (aka noise)
            while p < 1e-3:
                i += 1
                K2, p = stats.normaltest(S.loc[descendingGeneWeights.index[i:],component])
            #amount of while loops (i) gives the 'taken' genes
            selectedGenes = descendingGeneWeights.iloc[:i]
            
            for gene in selectedGenes.index:
                GeneID = gene              
                AnnotationRow = Additional_Annotation_File.loc[Additional_Annotation_File['GeneID'] == GeneID]
                Annotation = AnnotationRow["Annotation"].to_string(index=False)
                Start = AnnotationRow["Start"].to_string(index=False)
                End = AnnotationRow["End"].to_string(index=False)
                Strand = AnnotationRow["Strand"].to_string(index=False)
                dataDict = {'Genome':[Genome],'Component':[component],'GeneID':[GeneID],'Start':[Start],'End':[End],'Strand':[Strand],'Annotation':[Annotation]}
                toAppendDF = pd.DataFrame(data=dataDict)
                resultDF = resultDF.append(toAppendDF, ignore_index = True)
                
    return resultDF


def main():
    #Generates an ICA run for each genome and all genomes inside the featureCountFile
    args = parse_args()
    PCA_variance = float(args.PCA_Variance[0])
    Additional_Annotation_File = pd.read_table(args.additional_annotation_file[0] , sep='\t', header = 0)
    
    expressionDF = pd.read_table(args.FeatureCount_table[0] , sep='\t', header = 0)
    #split the long format into several dataframes each with its own Genome
    exprDFs = [pd.DataFrame(y) for x,y in expressionDF.groupby('Genome', as_index=False)]

    for i in range (1000):
        #run ICA 1K times, save result to csv
        oneRunDF = OneICA(PCA_variance,exprDFs,Additional_Annotation_File)
        oneRunDF.to_csv(F"results/8_ICA/{i}_geneModule.tsv", sep='\t', index = False)
    #Analyze the 1000 files created ; format x_genemodules.tsv
    
    geneModuleCountsDF = pd.DataFrame(columns=["Genome","GeneIDs","Count"])
    for geneModuleFile in glob.iglob(F'./results/8_ICA/*_geneModule.tsv'):
        print(F"Processing {geneModuleFile}")
        #Store all gene module combinations found for each genome, keep ones that show up atleast 500 times. 
        geneModuleDF = pd.read_table(geneModuleFile, sep = '\t', header = 0)
        seperatedGeneModuleDF = [pd.DataFrame(y) for x, y in geneModuleDF.groupby('Genome', as_index=False)]
        for DF in seperatedGeneModuleDF:
            groupedGeneModules = [pd.DataFrame(y) for x,y in DF.groupby('Component', as_index = False)]
            for geneModules in groupedGeneModules:
                #Have gene modules for the file/genome. Need to store this combination for each genome seperately and count their occurence
                Genome = geneModules.iloc[0]['Genome']
                GeneIDs = geneModules["GeneID"].tolist()
                GeneIDc = ""
                for gene in GeneIDs:
                    GeneIDc = GeneIDc + gene + ','
                thisGenome = geneModuleCountsDF[geneModuleCountsDF["Genome"] == Genome]
                if GeneIDc in thisGenome["GeneIDs"].values:
                    #already present +1 count
                    index = thisGenome.index[thisGenome['GeneIDs'] == GeneIDc].tolist()
                    geneModuleCountsDF.at[index[0],'Count'] = geneModuleCountsDF.loc[index[0],'Count'] + 1
                    continue
                #gene module not present in DF 
                dataDict = {'Genome':[Genome],'GeneIDs':[GeneIDc],'Count':[1]}
                toAppendDF = pd.DataFrame(data=dataDict)
                geneModuleCountsDF = geneModuleCountsDF.append(toAppendDF, ignore_index = True)
    geneModuleCountsDF.to_csv(F"results/8_ICA/GroupedgeneModule.tsv", sep='\t', index = False)
        

if __name__ == "__main__":
    main()
