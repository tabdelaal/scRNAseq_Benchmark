# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 14:35:09 2019

@author: Lieke
"""

import os
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


def run_Cell_BLAST(input_dir,output_dir,datafile,labfile,Rfile,numfeat=0,featfile=''):
    '''
    Run CellBlast
	
	Parameters
	----------
	input_dir : directory of the input files
	output_dir : directory of the output files
	datafile : name of the data file
    labfile : name of the label file
    Rfile : file to read the cross validation indices from
    numfeat : number of features to select, default = 0, which means that all features are used
    featfile : file with sorted features to read
    '''
    
    os.chdir(input_dir)
    
    # read the Rdata file
    robjects.r['load'](Rfile)

    nfolds = np.array(robjects.r['n_folds'], dtype = 'int')
    tokeep = np.array(robjects.r['Cells_to_Keep'], dtype = 'bool')
    col = np.array(robjects.r['col_Index'], dtype = 'int')
    col = col - 1 
    test_ind = np.array(robjects.r['Test_Idx'])
    train_ind = np.array(robjects.r['Train_Idx'])
    
    # read the feature file
    if (numfeat > 0):
        features = pd.read_csv(featfile,header=0,index_col=None, sep=',')

    # read the data and labels
    os.chdir(input_dir)
    data_old = cb.data.ExprDataSet.read_table(input_dir + datafile,orientation="cg", sep=",", index_col = 0, header = 0)
    labels = pd.read_csv(labfile, header=0,index_col=None, sep=',', usecols = col)
    
    data = cb.data.ExprDataSet(data_old.exprs[tokeep],data_old.obs.iloc[tokeep],data_old.var,data_old.uns)

    labels = gft(input_dir + labfile, dtype = "str", skip_header = 1, delimiter = ",", usecols = col)      
    labels = labels[tokeep]

    os.chdir(output_dir)
    
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
        
        if (numfeat > 0):
            feat_to_use = features.iloc[0:numfeat,i]
            train = train[:,feat_to_use]
            test = test[:,feat_to_use]

        
        train.obs['cell_type'] = y_train
                
        start = tm.time()
        train = train.normalize()
                
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
    os.chdir(output_dir)
    
    truelab = pd.DataFrame(truelab)
    pred = pd.DataFrame(pred)
            
    tr_time = pd.DataFrame(tr_time)
    ts_time = pd.DataFrame(ts_time)
    
    if (numfeat == 0):
        truelab.to_csv("Cell_BLAST_" + str(col) +"_true.csv", index = False)
        pred.to_csv("Cell_BLAST_" + str(col) +"_pred.csv", index = False)
        tr_time.to_csv("Cell_BLAST_" + str(col) +"_training_time.csv", index = False)
        ts_time.to_csv("Cell_BLAST_" + str(col) +"_test_time.csv", index = False)
    else:
        truelab.to_csv("Cell_BLAST_" + str(col) + "_" + str(numfeat) + "_true_" + featfile, index = False)
        pred.to_csv("Cell_BLAST_" + str(col) + "_" + str(numfeat) + "_pred_" + featfile, index = False)
        tr_time.to_csv("Cell_BLAST_" + str(col) + "_" + str(numfeat) + "_training_time_" + featfile, index = False)
        ts_time.to_csv("Cell_BLAST_" + str(col) + "_" + str(numfeat) + "_test_time_" + featfile, index = False)

        