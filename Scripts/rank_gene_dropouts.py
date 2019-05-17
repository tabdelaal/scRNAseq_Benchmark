# -*- coding: utf-8 -*-
"""
Created on Mon May  6 14:43:16 2019

@author: Lieke
"""

import os
import rpy2.robjects as robjects
import numpy as np
import pandas as pd
from sklearn import linear_model


def rank_gene_dropouts(datafile, input_dir, Rfile, CV = True):
    '''
    Rank the genes of the training set of the inputfile based on their dropout level
    This rank is written to a file.
    
    Parameters 
    ----------
    datafile : name of the datafile
    input_dir : directory of the files
    Rfile : file to read cross-validation indices from
    '''
    
    os.chdir(input_dir)
    
    # read the Rdata file
    robjects.r['load'](Rfile)

    nfolds = np.array(robjects.r['n_folds'], dtype = 'int')
    tokeep = np.array(robjects.r['Cells_to_Keep'], dtype = 'bool')
    train_ind = np.array(robjects.r['Train_Idx'])
    
    # read the data
    data = pd.read_csv(datafile,index_col=0,sep=',')
    data = data.iloc[tokeep]
    data = np.log2(data+1)
    
    if CV:
    
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
            print(np.shape(train_notzero))
            train_zero = train.iloc[:,zero]
            print(np.shape(train_zero))
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
            
    else:
        genes = np.zeros([np.shape(data)[1],1], dtype = '>U10')
        
        data.columns = np.arange(len(data.columns))
            
        # rank genes training set 
        dropout = data.shape[1] - data.astype(bool).sum(axis='rows')
        mean = data.mean(axis='rows')
            
        notzero = np.where((np.array(mean) > 0) & (np.array(dropout) > 0))[0]
        zero = np.where(~((np.array(mean) > 0) & (np.array(dropout) > 0)))[0]
        train_notzero = data.iloc[:,notzero]
        print(np.shape(train_notzero))
        train_zero = data.iloc[:,zero]
        print(np.shape(train_zero))
        zero_genes = train_zero.columns
            
        dropout = dropout.iloc[notzero]
        mean = mean.iloc[notzero]
    
        dropout = np.log2(np.array(dropout)).reshape(-1,1)
        mean = np.log2(np.array(mean)).reshape(-1,1)
        reg = linear_model.LinearRegression()
        reg.fit(mean,dropout)
    
        residuals = dropout - reg.predict(mean)
        residuals = pd.Series(np.array(residuals).ravel(),index=train_notzero.columns)
        residuals = residuals.sort_values(ascending=False)
        sorted_genes = residuals.index
        sorted_genes = sorted_genes.append(zero_genes)
            
        genes[:,0] = sorted_genes.values
  
        
    genes = pd.DataFrame(genes)
    if CV:
        genes.to_csv("rank_genes_dropouts.csv", index = False)
    else:
        genes.to_csv("rank_genes_dropouts_noCV.csv", index = False)

        

