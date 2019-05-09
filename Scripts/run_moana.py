# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 11:57:26 2019

@author: Lieke
"""

import os
import pandas as pd
import numpy as np
from moana.core import ExpMatrix
from moana.classify import CellTypeClassifier
import time as tm
import rpy2.robjects as robjects

def run_moana(input_dir,output_dir,datafile,labfile,Rfile):
    '''
    Run moana
    
    NOTE: at the moment it is only possible to run moana with a pretrained classifier,
    using the PBMC classifier is therefore hardcoded here
	
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

    tokeep = np.array(robjects.r['Cells_to_Keep'], dtype = 'bool')
    
    matrix = ExpMatrix.read_tsv(datafile, sep = ',')    
    matrix = matrix.iloc[tokeep]    
    
    data = ExpMatrix(X = np.transpose(matrix.X), genes = matrix.cells, cells = matrix.genes)
    data.genes.name = 'Genes'
    data.cells.name = 'Cells'
    data.index.name = 'Genes'
    data.columns.name = 'Cells'
    
    clf = CellTypeClassifier.read_pickle("moana_pbmc_classifier.pickle")
    
    start = tm.time()
    predictions = clf.predict(data)
    runtime = tm.time() - start
    
    np.asarray(predictions)
    
    pred = pd.DataFrame(predictions)
    
    os.chdir(output_dir)
    
    pred.to_csv("Moana_pred.csv", index = False)
    
    with open("Moana_time.csv", 'w') as f:
        f.write("%f\n" % runtime)

        
    
    
    
    
    
    
