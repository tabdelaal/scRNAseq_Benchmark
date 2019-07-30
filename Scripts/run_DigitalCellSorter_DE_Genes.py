import numpy as np
import pandas as pd
import scripts.DigitalCellSorter as DigitalCellSorter
import os
import time as tm
import rpy2.robjects as robjects

def run_DigitalCellSorter_DE_Genes(DataPath, LabelsPath, CV_RDataPath, GeneListPath, OutputDir, NumGenes = 20):
    '''
    run DigitalCellSorter
    Wrapper script to run DigitalCellSorter on a benchmark dataset using a predefined genelist,
    outputs lists of true and predicted cell labels as csv files, as well as computation time.  
  
    Parameters
    ----------
    DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
    as row names and gene names as column names.
    LabelsPath : Cell population annotations file path (.csv).
    CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
    GeneListPath : Data file path to the genelist (the _NUMGENES_FOLD_I.csv part should be omitted as this is varied each time)
    OutputDir : Output directory defining the path of the exported file.
    '''
        
    # read the Rdata file
    robjects.r['load'](CV_RDataPath)

    nfolds = np.array(robjects.r['n_folds'], dtype = 'int')
    tokeep = np.array(robjects.r['Cells_to_Keep'], dtype = 'bool')
    col = np.array(robjects.r['col_Index'], dtype = 'int')
    col = col - 1 
    test_ind = np.array(robjects.r['Test_Idx'])
    train_ind = np.array(robjects.r['Train_Idx'])

    # read the data
    data = pd.read_csv(DataPath,index_col=0,sep=',')
    labels = pd.read_csv(LabelsPath, header=0,index_col=None, sep=',', usecols = col)
    
    labels = labels.iloc[tokeep]
    data = data.iloc[tokeep]
    
    ts_time=[]
    truelab = []
    predlab = []

    
    for i in range(np.squeeze(nfolds)):
        test_ind_i = np.array(test_ind[i], dtype = 'int') - 1
    
        test=data.iloc[test_ind_i]
        y_test=labels.iloc[test_ind_i]
            
        test = test.transpose()
        
        n_clusters = len(np.unique(labels))
        AvailableCPUsCount = 1
        N_samples_for_distribution = 10000
        
        GeneList = GeneListPath + 'DCS_Marker_Genes_' + str(NumGenes) + '_Fold_' + str(i+1) + '.xlsx'
        DataName = 'DCS_' + str(NumGenes) + '_Fold_' + str(i+1) 
        
        start = tm.time()
        pred = DigitalCellSorter.DigitalCellSorter().Process(test, DataName, 
                                                saveDir = OutputDir, 
                                                geneListFileName = GeneList,
                                                N_samples_for_distribution = N_samples_for_distribution,
                                                AvailableCPUsCount = AvailableCPUsCount,
                                                clusterIndex=None,
                                                clusterName=None,
                                                n_clusters=n_clusters,
                                                marker_expression_plot = False,
                                                marker_subplot = False,
                                                tSNE_cluster_plot = False)	
    
        os.chdir(OutputDir)
    
        results = pd.read_excel(DataName + '_voting.xlsx',header=0,index_col=None)
        results = results.iloc[:,-1]

        prediction = np.zeros(np.shape(pred), dtype='>U10')
    
        for i in range(len(results)):
        	prediction[np.where(pred == i)] = results.values[i]
        
        ts_time.append(tm.time()-start)
            
        truelab.extend(y_test.values)
        predlab.extend(prediction)
    
    truelab = pd.DataFrame(truelab)
    predlab = pd.DataFrame(predlab)
        
    ts_time = pd.DataFrame(ts_time)

    
    truelab.to_csv("DigitalCellSorter_" + str(NumGenes) + "_True_Labels.csv", index = False)
    predlab.to_csv("DigitalCellSorter_" + str(NumGenes) + "_Pred_Labels.csv", index = False)
    ts_time.to_csv("DigitalCellSorter_" + str(NumGenes) + "_Testing_Time.csv", index = False)

            

        