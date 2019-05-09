# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 15:32:29 2019

@author: Lieke
"""

import numpy as np
import pandas as pd
import scripts.DigitalCellSorter as DigitalCellSorter
import os
import time as tm
import rpy2.robjects as robjects

def run_DigitalCellSorter(input_dir,output_dir,datafile, labfile, Rfile, numfeat = 0, featfile = ''):
    '''
    Run DigitalCellSorter
    
    NOTE: a marker genelist is needed to run DigitalCellSorter
    at the moment only the genelist for the PBMC data is publically available
    this one is therefore hardcoded here

	
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

    tokeep = np.array(robjects.r['Cells_to_Keep'], dtype = 'bool')
    
    # read the data
    data = pd.read_csv(datafile,index_col=0,sep=',')
    data = data.iloc[tokeep]

    # read the feature file
    if (numfeat > 0):
        features = pd.read_csv(featfile,header=0,index_col=None, sep=',')
        feat_to_use = features.iloc[0:numfeat,0]
        data = data.iloc[:,feat_to_use]
        
    data = data.transpose()
    
    # number of different cell types in the data?
    n_clusters = 8
    
    AvailableCPUsCount = 1
    N_samples_for_distribution = 10000
    cellTypeNameForSubclustering = None #None  'B cells'  'T cells'

    if cellTypeNameForSubclustering == 'B cells':
        clusterIndex = [0]
        geneListToUse = 'DigitalCellSorter/CIBERSORT_B_SUB.xlsx'
    elif cellTypeNameForSubclustering == 'T cells':
        clusterIndex = [2,6,7]
        geneListToUse = 'DigitalCellSorter/CIBERSORT_T_SUB.xlsx'
    else:
        clusterIndex = None
        geneListToUse = 'DigitalCellSorter/CIBERSORT.xlsx'
        
    start = tm.time()
    pred = DigitalCellSorter.DigitalCellSorter().Process(data, 'DigitalCellSorter_Zhang', 
                                                saveDir = output_dir, 
                                                geneListFileName = geneListToUse,
                                                N_samples_for_distribution = N_samples_for_distribution,
                                                AvailableCPUsCount = AvailableCPUsCount,
                                                clusterIndex=clusterIndex,
                                                clusterName=cellTypeNameForSubclustering,
                                                n_clusters=n_clusters)	
    runtime = tm.time() - start 
    
    os.chdir(output_dir)
    
    results = pd.read_excel('DigitalCellSorter_Zhang_voting.xlsx',header=0,index_col=None, usecols=[11])

    prediction = np.zeros(np.shape(pred), dtype='>U10')
    
    for i in range(len(results)):
    	prediction[np.where(pred == i)] = results.values[i]
    
    prediction = pd.DataFrame(prediction)
        
    if (numfeat == 0):
        prediction.to_csv('DigitalCellSorter_pred.csv', index = None)
        
        with open("DigitallCellSorter_time.csv", 'w') as f:
            f.write("%f\n" % runtime)
    
    else:
        prediction.to_csv('DigitalCellSorter_' + str(numfeat) + "_pred_" + featfile, index = None)
        
        with open("DigitallCellSorter_" + str(numfeat) + "_time_" + featfile, 'w') as f:
            f.write("%f\n" % runtime)

            

        