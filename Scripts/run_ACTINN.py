# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 15:34:00 2019

@author: Lieke
"""

import os 
import numpy as np
import pandas as pd
import time as tm
import rpy2.robjects as robjects

def run_ACTINN(input_dir, output_dir, datafile, labfile, Rfile):
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
    data = pd.read_csv(datafile,index_col=0,sep=',')
    labels = pd.read_csv(labfile, header=0,index_col=None, sep=',', usecols = col)
    
    labels = labels.iloc[tokeep]
    data = data.iloc[tokeep]
    
    # folder with results
    os.chdir(output_dir)
    
    tot=[]
    truelab = []
    pred = []

    for i in range(np.squeeze(nfolds)):
        test_ind_i = np.array(test_ind[i], dtype = 'int') - 1
        train_ind_i = np.array(train_ind[i], dtype = 'int') - 1
    
        train=data.iloc[train_ind_i]
        test=data.iloc[test_ind_i]
        y_train=labels.iloc[train_ind_i]
        y_test=labels.iloc[test_ind_i]
        
        train = train.transpose()
        test = test.transpose()
        
        train.to_csv("train.csv")
        test.to_csv("test.csv")
        y_train.to_csv("train_lab.csv", header = False, index = True, sep = '\t')
        y_test.to_csv("test_lab.csv", header = False, index = True, sep = '\t')
            
        os.system("python /home/nfs/lcmmichielsen/classifiers/ACTINN/actinn_format.py -i train.csv -o train -f csv")
        os.system("python /home/nfs/lcmmichielsen/classifiers/ACTINN/actinn_format.py -i test.csv -o test -f csv")
        
        start = tm.time()
        os.system("python /home/nfs/lcmmichielsen/classifiers/ACTINN/actinn_predict.py -trs train.h5 -trl train_lab.csv -ts test.h5")    
        tot.append(tm.time()-start)

        truelab.extend(y_test.values)
        predlabels = pd.read_csv('predicted_label.txt',header=0,index_col=None, sep='\t', usecols = [1])            
        pred.extend(predlabels.values)
    
            
    truelab = pd.DataFrame(truelab)
    pred = pd.DataFrame(pred)
    tot_time = pd.DataFrame(tot)
    
    truelab.to_csv("ACTINN_" + str(col) +"_true.csv", index = False)
    pred.to_csv("ACTINN_" + str(col) +"_pred.csv", index = False)
    tot_time.to_csv("ACTINN_" + str(col) +"_time.csv", index = False)