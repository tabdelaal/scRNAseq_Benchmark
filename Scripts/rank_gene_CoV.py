# -*- coding: utf-8 -*-
"""
Created on Tue May  7 14:36:27 2019

@author: Lieke
"""

# -*- coding: utf-8 -*-
"""
Created on Mon May  6 14:43:16 2019

@author: Lieke
"""

import os
from pathlib import Path
import math

import rpy2.robjects as robjects
import numpy as np
import pandas as pd
import statsmodels.api as sm


def rank_gene_CoV(datafile, Rfile, output_dir, CV = True):
    '''
    Rank the genes of the training set of the inputfile based on their dropout level
    This rank is written to a file.

    Parameters
    ----------
    datafile : name of the datafile
    input_dir : directory of the files
    Rfile : file to read cross-validation indices from
    '''

    # read the Rdata file
    robjects.r['load'](Rfile)

    nfolds = np.array(robjects.r['n_folds'], dtype = 'int')
    tokeep = np.array(robjects.r['Cells_to_Keep'], dtype = 'bool')
    train_ind = np.array(robjects.r['Train_Idx'])

    # read the data
    data = pd.read_csv(datafile,index_col=0,sep=',')
    data = data.iloc[tokeep]

    if CV:

        genes = np.zeros([np.shape(data)[1],np.squeeze(nfolds)], dtype = '>U10')

        for i in range(np.squeeze(nfolds)):
            train_ind_i = np.array(train_ind[i], dtype = 'int') - 1
            train=data.iloc[train_ind_i]
            train.columns = np.arange(len(train.columns))

            # rank genes training set
            mean = train.mean(axis='rows')
            variance = train.var(axis='rows')

            notzero = np.where(np.array(mean) > 0)[0]
            zero = np.where(~(np.array(mean) > 0))[0]
            train_notzero = train.iloc[:,notzero]
            train_zero = train.iloc[:,zero]
            zero_genes = train_zero.columns

            mean = mean.iloc[notzero]
            variance = variance.iloc[notzero]

            CoV = variance.apply(math.sqrt).divide(mean)

            gamma_model = sm.GLM(mean, CoV, family=sm.families.Gamma())
            gamma_results = gamma_model.fit()

            residuals = CoV - gamma_results.predict(mean)
            residuals = pd.Series(np.array(residuals).ravel(),index=train_notzero.columns)
            residuals = residuals.sort_values(ascending=False)
            sorted_genes = residuals.index
            sorted_genes = sorted_genes.append(zero_genes)

            genes[:,i] = sorted_genes.values

    else:

        genes = np.zeros([np.shape(data)[1],1], dtype = '>U10')

        data.columns = np.arange(len(data.columns))

        # rank genes training set
        mean = data.mean(axis='rows')
        variance = data.var(axis='rows')

        notzero = np.where(np.array(mean) > 0)[0]
        zero = np.where(~(np.array(mean) > 0))[0]
        train_notzero = data.iloc[:,notzero]
        train_zero = data.iloc[:,zero]
        zero_genes = train_zero.columns

        mean = mean.iloc[notzero]
        variance = variance.iloc[notzero]

        CoV = variance.apply(math.sqrt).divide(mean)

        gamma_model = sm.GLM(mean, CoV, family=sm.families.Gamma())
        gamma_results = gamma_model.fit()

        residuals = CoV - gamma_results.predict(mean)
        residuals = pd.Series(np.array(residuals).ravel(),index=train_notzero.columns)
        residuals = residuals.sort_values(ascending=False)
        sorted_genes = residuals.index
        sorted_genes = sorted_genes.append(zero_genes)

        genes[:,0] = sorted_genes.values

    genes = pd.DataFrame(genes)

    genes.to_csv(str(output_dir / Path("rank_genes_CoV.csv")), index = False)


rank_gene_CoV(argv[1], argv[2], argv[3])
