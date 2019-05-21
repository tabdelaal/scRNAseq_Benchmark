import os
from sys import argv
from pathlib import Path

import rpy2.robjects as robjects
import numpy as np
import pandas as pd
from sklearn import linear_model


def rank_gene_dropouts(DataPath, CV_RDataPath, OutputDir):
    '''
    Script to rank the genes in the training set of the inputfile based on their dropout level.
    This rank is written to a file.

    Parameters
    ----------
    DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes
    CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
    OutputDir : Output directory defining the path of the exported file.
    '''

    # read the Rdata file
    robjects.r['load'](CV_RDataPath)

    nfolds = np.array(robjects.r['n_folds'], dtype = 'int')
    tokeep = np.array(robjects.r['Cells_to_Keep'], dtype = 'bool')
    train_ind = np.array(robjects.r['Train_Idx'])

    # read the data
    data = pd.read_csv(DataPath,index_col=0,sep=',')
    data = data.iloc[tokeep]
    data = np.log2(data+1)

    genes = np.zeros([np.shape(data)[1],np.squeeze(nfolds)], dtype = '>U10')

    for i in range(np.squeeze(nfolds)):
        train_ind_i = np.array(train_ind[i], dtype = 'int') - 1
        train=data.iloc[train_ind_i]
        train.columns = np.arange(len(train.columns))

        # rank genes training set
        dropout = (train == 0).sum(axis='rows')
        dropout = (dropout / train.shape[0]) * 100
        mean = train.mean(axis='rows')

        notzero = np.where((np.array(mean) > 0) & (np.array(dropout) > 0))[0]
        zero = np.where(~((np.array(mean) > 0) & (np.array(dropout) > 0)))[0]
        train_notzero = train.iloc[:,notzero]
        train_zero = train.iloc[:,zero]
        zero_genes = train_zero.columns

        dropout = dropout.iloc[notzero]
        mean = mean.iloc[notzero]

        dropout = np.log2(np.array(dropout)).reshape(-1,1)
        mean = np.array(mean).reshape(-1,1)
        reg = linear_model.LinearRegression()
        reg.fit(mean,dropout)

        residuals = dropout - reg.predict(mean)
        residuals = pd.Series(np.array(residuals).ravel(),index=train_notzero.columns)
        residuals = residuals.sort_values(ascending=False)
        sorted_genes = residuals.index
        sorted_genes = sorted_genes.append(zero_genes)

        genes[:,i] = sorted_genes.values


    genes = pd.DataFrame(genes)

    genes.to_csv(str(OutputDir / Path("rank_genes_dropouts.csv")), index = False)

rank_gene_dropouts(argv[1], argv[2], argv[3])
