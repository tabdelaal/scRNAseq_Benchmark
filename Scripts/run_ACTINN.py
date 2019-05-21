import os 
import numpy as np
import pandas as pd
import time as tm
import rpy2.robjects as robjects

def run_ACTINN(DataPath, LabelsPath, CV_RDataPath, OutputDir, GeneOrderPath = "", NumGenes = 0):
    '''
    run ACTINN
    Wrapper script to run ACTINN on a benchmark dataset with 5-fold cross validation,
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

    # read the data
    data = pd.read_csv(DataPath,index_col=0,sep=',')
    labels = pd.read_csv(LabelsPath, header=0,index_col=None, sep=',', usecols = col)
    
    labels = labels.iloc[tokeep]
    data = data.iloc[tokeep]
    
    # read the feature file
    if (NumGenes > 0):
        features = pd.read_csv(GeneOrderPath,header=0,index_col=None, sep=',')
    
    # folder with results
    os.chdir(OutputDir)
    
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
        
        if (NumGenes > 0):
            feat_to_use = features.iloc[0:NumGenes,i]
            train = train.iloc[:,feat_to_use]
            test = test.iloc[:,feat_to_use]
        
        train = train.transpose()
        test = test.transpose()
        
        train.to_csv("train.csv")
        test.to_csv("test.csv")
        y_train.to_csv("train_lab.csv", header = False, index = True, sep = '\t')
        y_test.to_csv("test_lab.csv", header = False, index = True, sep = '\t')
        
        tm.sleep(60)
            
        os.system("python /home/nfs/lcmmichielsen/classifiers/ACTINN/actinn_format.py -i train.csv -o train -f csv")
        os.system("python /home/nfs/lcmmichielsen/classifiers/ACTINN/actinn_format.py -i test.csv -o test -f csv")
        
        start = tm.time()
        os.system("python /home/nfs/lcmmichielsen/classifiers/ACTINN/actinn_predict.py -trs train.h5 -trl train_lab.csv -ts test.h5")    
        tot.append(tm.time()-start)
        
        tm.sleep(60)

        truelab.extend(y_test.values)
        predlabels = pd.read_csv('predicted_label.txt',header=0,index_col=None, sep='\t', usecols = [1])            
        pred.extend(predlabels.values)
    
            
    truelab = pd.DataFrame(truelab)
    pred = pd.DataFrame(pred)
    tot_time = pd.DataFrame(tot)
    
    if (NumGenes == 0):  
        truelab.to_csv("ACTINN_True_Labels.csv", index = False)
        pred.to_csv("ACTINN_Pred_Labels.csv", index = False)
        tot_time.to_csv("ACTINN_Total_Time.csv", index = False)
    else:
        truelab.to_csv("ACTINN_" + str(NumGenes) + "_True_Labels.csv", index = False)
        pred.to_csv("ACTINN_" + str(NumGenes) + "_Pred_Labels.csv", index = False)
        tot_time.to_csv("ACTINN_" + str(NumGenes) + "_Total_Time.csv", index = False)
        
        
        
        
        
        
        
