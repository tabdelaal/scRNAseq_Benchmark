# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 13:35:53 2019

@author: Lieke
"""

import os
import numpy as np
import pandas as pd
import time as tm
from pathlib import Path

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import NearestCentroid
from sklearn.neighbors import KNeighborsClassifier
import rpy2.robjects as robjects


def run_baseline(output_dir, datafile, labfile, Rfile,
        classifiers):
    '''
    Run baseline classifiers: NMC, LDA, kNN, SVM, RF

	Parameters
	----------
	output_dir : directory of the output files
	datafile : name of the data file
    labfile : name of the label file
    Rfile : file to read the cross validation indices from
    '''

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
    labels = pd.read_csv(labfile, header=0,index_col=None, sep=',',
        usecols = col)

    labels = labels.iloc[tokeep]
    data = data.iloc[tokeep]

    # normalize data
    data = np.log1p(data)

    Classifiers = {'NMC': NearestCentroid,
                   'LDA': LinearDiscriminantAnalysis,
                   'kNN': KNeighborsClassifier,
                   'SVM': LinearSVC, 'RF': RandomForestClassifier}
    Classifiers = {k: Classifiers[k] for k in classifiers}

    for c in Classifiers:
        tr_time=[]
        ts_time=[]
        truelab = []
        pred = []

        for i in range(np.squeeze(nfolds)):
            test_ind_i = np.array(test_ind[i], dtype = 'int') - 1
            train_ind_i = np.array(train_ind[i], dtype = 'int') - 1

            train=data.iloc[train_ind_i]
            test=data.iloc[test_ind_i]
            y_train=labels.iloc[train_ind_i]
            y_test=labels.iloc[test_ind_i]


            if c == 'RF': #FIXME
                Classifier=Classifiers[c](n_estimators = 50)
            elif c == 'kNN':
                Classifier=Classifiers[c](n_neighbors = 50)
            else:
                Classifier=Classifiers[c]()

            start=tm.time()
            Classifier.fit(train, y_train)
            tr_time.append(tm.time()-start)

            start=tm.time()
            predicted = Classifier.predict(test)
            ts_time.append(tm.time()-start)

            truelab.extend(y_test.values)
            pred.extend(predicted)

        mean_tr = np.mean(tr_time)
        mean_ts = np.mean(ts_time)

        truelab = pd.DataFrame(truelab)
        pred = pd.DataFrame(pred)

        output_dir = Path(output_dir)
        truelab.to_csv(str(output_dir / Path(f"{c}_true.csv")),
                       index = False, header = False)
        pred.to_csv(str(output_dir / Path(f"{c}_pred.csv")),
                    index = False, header = False)

        with (output_dir / Path(f"{c}_training_time.csv")).open(mode="w") as f:
            f.write("%f\n" % mean_tr)

        with (output_dir / Path(f"{c}_test_time.csv")).open(mode="w") as f:
            f.write("%f\n" % mean_ts)
