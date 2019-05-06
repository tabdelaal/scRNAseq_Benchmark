# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 16:11:55 2019

@author: Lieke
"""

from scvi.dataset import CsvDataset
import os
import numpy as np
import pandas as pd
from scvi.models import SCANVI
from scvi.inference import SemiSupervisedTrainer
import time as tm
import rpy2.robjects as robjects

def run_scVI(input_dir, output_dir, datafile, labfile, Rfile):
    '''
    Run scVI
	
	Parameters
	----------
	input_dir : directory of the input files
	output_dir : directory of the output files
	datafile : name of the data file
    labfile : name of the label file
    Rfile : file to read the cross validation indices from
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

    # read the data
    os.chdir(input_dir)
    data = pd.read_csv(datafile,index_col=0,sep=',')
    labels = pd.read_csv(labfile, header=0,index_col=None, sep=',', usecols = col)

    labels = labels.iloc[tokeep]
    data = data.iloc[tokeep]    
    
    #save labels as csv file with header and index column
    labels.to_csv('Labels_scvi.csv')
    data.to_csv('Data_scvi.csv')    
    
    train = CsvDataset('Data_scvi.csv', save_path = input_dir, sep = ",", labels_file = "Labels_scvi.csv", gene_by_cell = False)
    
    ## this semisupervised trainer automatically uses a part of the input data for training and a part for testing
    scanvi = SCANVI(train.nb_genes, train.n_batches, train.n_labels)
    trainer_scanvi = SemiSupervisedTrainer(scanvi, train, frequency=5)
    
    n_epochs = 200
    
    truelab = []
    pred = []
    tr_time = []
    ts_time = []
    
    for i in range(np.squeeze(nfolds)):
        test_ind_i = np.array(test_ind[i], dtype = 'int') - 1
        train_ind_i = np.array(train_ind[i], dtype = 'int') - 1

        trainer_scanvi.labelled_set = trainer_scanvi.create_posterior(indices=(train_ind_i).ravel(), shuffle = False)
        trainer_scanvi.labelled_set.to_monitor = ['ll','accuracy']
        trainer_scanvi.unlabelled_set = trainer_scanvi.create_posterior(indices=(test_ind_i).ravel(), shuffle = False)
        trainer_scanvi.unlabelled_set.to_monitor = ['ll','accuracy']
    
        start = tm.time()
        trainer_scanvi.train(n_epochs)
        tr_time.append(tm.time()-start)
    
        ## labels of test set are in y_pred
        ## labels are returned in numbers, should be mapped back to the real labels
        ## indices are permutated
        start = tm.time()
        y_true, y_pred = trainer_scanvi.unlabelled_set.compute_predictions()
        ts_time.append(tm.time()-start)
        
        truelab.extend(y_true)
        pred.extend(y_pred)
    
    #write results
    os.chdir(output_dir)
    
    truelab = pd.DataFrame(truelab)
    pred = pd.DataFrame(pred)
    
    tr_time = pd.DataFrame(tr_time)
    ts_time = pd.DataFrame(ts_time)

    
    truelab.to_csv("scVI_" + str(col) +"_true.csv", index = False)
    pred.to_csv("scVI_" + str(col) +"_pred.csv", index = False)
    
    tr_time.to_csv("scVI_" + str(col) +"_training_time.csv", index = False)
    ts_time.to_csv("scVI_" + str(col) +"_test_time.csv", index = False)


