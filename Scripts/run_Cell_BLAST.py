import os
from sys import argv
from pathlib import Path
import time as tm
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

import tensorflow as tf
tf.logging.set_verbosity(0)

import Cell_BLAST as cb
import numpy as np
from numpy import genfromtxt as gft
import rpy2.robjects as robjects


def run_Cell_BLAST(DataPath, LabelsPath, CV_RDataPath, OutputDir, GeneOrderPath = "", NumGenes = 0):
    '''
    run Cell_BLAST
    Wrapper script to run Cell_BLAST on a benchmark dataset with 5-fold cross validation,
    outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
    Parameters
    ----------
    DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
    as row names and gene names as column names.
    LabelsPath : Cell population annotations file path (.csv).
    CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
    OutputDir : Output directory defining the path of the exported file.
    GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
    defining the genes order for each cross validation fold, default is NULL.
    NumGenes : Number of genes used in case of feature selection (integer), default is 0.
    '''
        
    # read the Rdata file
    robjects.r['load'](CV_RDataPath)

    nfolds = np.array(robjects.r['n_folds'], dtype = 'int')
    tokeep = np.array(robjects.r['Cells_to_Keep'], dtype = 'bool')
    col = np.array(robjects.r['col_Index'], dtype = 'int')
    col = col - 1 
    test_ind = np.array(robjects.r['Test_Idx'])
    train_ind = np.array(robjects.r['Train_Idx'])
    
    # read the feature file
    if (NumGenes > 0):
        features = pd.read_csv(GeneOrderPath,header=0,index_col=None, sep=',')

    # read the data and labels
    data_old = cb.data.ExprDataSet.read_table(DataPath,orientation="cg", sep=",", index_col = 0, header = 0, sparsify = True).normalize()
    labels = pd.read_csv(LabelsPath, header=0,index_col=None, sep=',', usecols = col)
    
    data = cb.data.ExprDataSet(data_old.exprs[tokeep],data_old.obs.iloc[tokeep],data_old.var,data_old.uns)

    labels = gft(LabelsPath, dtype = "str", skip_header = 1, delimiter = ",", usecols = col)      
    labels = labels[tokeep]
   
    truelab = []
    pred = []
    tr_time = []
    ts_time = []
    
    for i in range(np.squeeze(nfolds)):
        test_ind_i = np.array(test_ind[i], dtype = 'int') - 1
        train_ind_i = np.array(train_ind[i], dtype = 'int') - 1

        train=data[train_ind_i,:]
        test=data[test_ind_i,:]
        y_train = labels[train_ind_i]
        y_test = labels[test_ind_i]
        
        if (NumGenes > 0):
            feat_to_use = features.iloc[0:NumGenes,i]
            train = train[:,feat_to_use]
            test = test[:,feat_to_use]

        
        train.obs['cell_type'] = y_train
                
        start = tm.time()
                
        # reduce dimensions
        num_epoch = 50
        models = []
    
        for j in range(4):
            models.append(cb.directi.fit_DIRECTi(train, epoch=num_epoch, patience=10, random_seed = j, path="%d" % j))
    
        # train model
        blast = cb.blast.BLAST(models, train).build_empirical()
        tr_time.append(tm.time()-start)
        
        # predict labels
        start = tm.time()
        test_pred = blast.query(test).annotate('cell_type')
        ts_time.append(tm.time()-start)

        truelab.extend(y_test)
        pred.extend(test_pred.values)
    
    #write results    
    truelab = pd.DataFrame(truelab)
    pred = pd.DataFrame(pred)
            
    tr_time = pd.DataFrame(tr_time)
    ts_time = pd.DataFrame(ts_time)

    truelab.to_csv(str(Path(OutputDir+"/Cell_BLAST_true.csv")),index = False)
    pred.to_csv(str(Path(OutputDir+"/Cell_BLAST_pred.csv")),index = False)
    tr_time.to_csv(str(Path(OutputDir+"/Cell_BLAST_training_time.csv")), index = False)
    ts_time.to_csv(str(Path(OutputDir+"/Cell_BLAST_test_time.csv")),index = False)


run_Cell_BLAST(argv[1], argv[2], argv[3], argv[4], argv[5], int(argv[6]))
