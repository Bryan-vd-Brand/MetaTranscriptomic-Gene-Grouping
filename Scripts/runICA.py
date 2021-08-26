import os
import pandas as pd
import argparse
from sklearn.decomposition import FastICA,PCA
import numpy as np

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
    return parser.parse_args()



def main():
    #Generates an ICA run for each genome and all genomes inside the featureCountFile
    args = parse_args()
    PCA_variance = float(args.PCA_Variance[0])
    
    expressionDF = pd.read_table(args.FeatureCount_table[0] , sep='\t', header = 0)
    #split the long format into several dataframes each with its own Genome
    exprDFs = [pd.DataFrame(y) for x,y in expressionDF.groupby('Genome', as_index=False)]
    
    for exprDF in exprDFs:
        Genome = exprDF['Genome'].iloc[0]
        #run for each genome
        wideDF = exprDF.pivot(index="RunAccession", columns="GeneID", values="RPKM").fillna(0) #This is the data for the X matrix in ICA
        #Now the data needs to be centered, i.e substract the mean for each variable (gene) from the value. 
        #because we have many genes with NA value the mean is calculated with counting NA entries as 0.
        for column in wideDF.columns:
            mean = wideDF[F"{column}"].mean()
            wideDF[F"{column}"] = wideDF[F"{column}"].apply(lambda x : x - mean)
        print(Genome)
        #rename to X for consistency with ICA paper
        X = wideDF
        #Use PCA to reduce the dimensions of the ICA algorithm
        pca = PCA().fit(X)
        #explained_variance_ Equal to n_components largest eigenvalues of the covariance matrix of X ; ndarray      
        cumulative_variance = np.cumsum(pca.explained_variance_ratio_)
        number_of_components = np.where(cumulative_variance > PCA_variance)[0][0] + 1
        print(cumulative_variance)
        print(number_of_components)
        
        #Run ICA
        X_test = X.transpose()
        ica = FastICA(whiten=True,max_iter=1000,tol=1e-6,n_components = number_of_components)
        S = pd.DataFrame(ica.fit_transform(X_test),index=X_test.index)
        A = pd.DataFrame(ica.mixing_,index=X_test.columns)
        #save S&A 
        S.to_csv(F"results/8_ICA/{Genome}_S.tsv", sep = '\t')
        A.to_csv(F"results/8_ICA/{Genome}_A.tsv", sep = '\t')

    

if __name__ == "__main__":
    main()
