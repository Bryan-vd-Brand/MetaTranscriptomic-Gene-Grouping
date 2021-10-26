import os
import pandas as pd
import argparse
from sklearn.decomposition import FastICA,PCA
from sklearn import preprocessing
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
    #Generates an ICA run for all genomes at the same time
    args = parse_args()
    PCA_variance = float(args.PCA_Variance[0])
    
    expressionDF = pd.read_table(args.FeatureCount_table[0] , sep='\t', header = 0)
    #split the long format into several dataframes each with its own Genome
    #rename to X for consistency with ICA paper
    X = expressionDF.pivot(index=["RunAccession","Genome"], columns="GeneID", values="RPKM").fillna(0)
    #scale using MinMaxScaler, not centered since many values are NA/0/Sparse
    min_max_scaler = preprocessing.MinMaxScaler()
    X_scaled = min_max_scaler.fit_transform(X)
    #Use PCA to reduce the dimensions of the ICA algorithm
    pca = PCA().fit(X_scaled)
    #explained_variance_ Equal to n_components largest eigenvalues of the covariance matrix of X ; ndarray      
    cumulative_variance = np.cumsum(pca.explained_variance_ratio_)
    number_of_components = np.where(cumulative_variance > PCA_variance)[0][0] + 1 #0-based -> 1-based
    print(cumulative_variance)
    print(number_of_components)
    
    #Run ICA
    X_t = X_scaled.transpose()
    ica = FastICA(whiten=True,max_iter=1000,tol=1e-6,n_components = number_of_components)
    S = pd.DataFrame(ica.fit_transform(X_t),index=X.columns)
    A = pd.DataFrame(ica.mixing_,index=X.index)
    #save S&A 
    S.to_csv(F"results/8_ICA/AllICA_S.tsv", sep = '\t')
    A.to_csv(F"results/8_ICA/AllICA_A.tsv", sep = '\t')

    

if __name__ == "__main__":
    main()
