import numpy as np
import pandas as pd
import scripts.DigitalCellSorter as DigitalCellSorter
import os
import time as tm
import rpy2.robjects as robjects

def run_DigitalCellSorter(DataPath, LabelsPath, CV_RDataPath, GeneListPath, OutputDir, GeneOrderPath = "", NumGenes = 0):
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
    GeneListPath : Data file path to the genest.
    OutputDir : Output directory defining the path of the exported file.
    GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
    defining the genes order for each cross validation fold, default is NULL.
    NumGenes : Number of genes used in case of feature selection (integer), default is 0.
    '''
        
    # read the Rdata file
    robjects.r['load'](CV_RDataPath)

    tokeep = np.array(robjects.r['Cells_to_Keep'], dtype = 'bool')
    col = np.array(robjects.r['col_Index'], dtype = 'int')
    col = col - 1
    
    # read the data
    data = pd.read_csv(DataPath,index_col=0,sep=',')
    data = data.iloc[tokeep]
    
    truelab = pd.read_csv(LabelsPath, header=0,index_col=None, sep=',', usecols = col)
    truelab = truelab.iloc[tokeep]


    # read the feature file
    if (NumGenes > 0):
        features = pd.read_csv(GeneOrderPath,header=0,index_col=None, sep=',')
        feat_to_use = features.iloc[0:NumGenes,0]
        data = data.iloc[:,feat_to_use]
        
    data = data.transpose()
    
    # number of different cell types in the data?
    n_clusters = 8
    AvailableCPUsCount = 1
    N_samples_for_distribution = 10000
        
    start = tm.time()
    pred = DigitalCellSorter.DigitalCellSorter().Process(data, 'DigitalCellSorter_Zhang', 
                                                saveDir = OutputDir, 
                                                geneListFileName = GeneListPath,
                                                N_samples_for_distribution = N_samples_for_distribution,
                                                AvailableCPUsCount = AvailableCPUsCount,
                                                clusterIndex=None,
                                                clusterName=None,
                                                n_clusters=n_clusters)	
    runtime = tm.time() - start 
    
    os.chdir(OutputDir)
    
    results = pd.read_excel('DigitalCellSorter_Zhang_voting.xlsx',header=0,index_col=None, usecols=[11])

    prediction = np.zeros(np.shape(pred), dtype='>U10')
    
    for i in range(len(results)):
    	prediction[np.where(pred == i)] = results.values[i]
    
    prediction = pd.DataFrame(prediction)
        
    if (NumGenes == 0):  
        truelab.to_csv("DigitalCellSorter_True_Labels.csv", index = False)
        prediction.to_csv("DigitalCellSorter_Pred_Labels.csv", index = False)
        with open("DigitalCellSorter_Total_Time.csv", 'w') as f:
            f.write("%f\n" % runtime)
    else:
        truelab.to_csv("DigitalCellSorter_" + str(NumGenes) + "_True_Labels.csv", index = False)
        prediction.to_csv("DigitalCellSorter_" + str(NumGenes) + "_Pred_Labels.csv", index = False)
        with open("DigitalCellSorter_" + str(NumGenes) + "_Total_Time.csv", 'w') as f:
            f.write("%f\n" % runtime)

            

        