import os
import pandas as pd
import numpy as np
from moana.core import ExpMatrix
from moana.classify import CellTypeClassifier
import time as tm
import rpy2.robjects as robjects

def run_moana(DataPath, LabelsPath, ClassifierPath, OutputDir, GeneOrderPath = "", NumGenes = 0):
    '''
    run moana
    Wrapper script to run moana on a benchmark dataset with a pretrained classifier,
    outputs lists of true and predicted cell labels as csv files, as well as computation time.  
  
    Parameters
    ----------
    DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
    as row names and gene names as column names.
    LabelsPath : Cell population annotations file path (.csv).
    CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
    ClassifierPath : Data file path to the pretrained classifier.
    OutputDir : Output directory defining the path of the exported file.
    GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
    defining the genes order for each cross validation fold, default is NULL.
    NumGenes : Number of genes used in case of feature selection (integer), default is 0.
    '''
    
#    # read the Rdata file
#    robjects.r['load'](CV_RDataPath)
#
#    tokeep = np.array(robjects.r['Cells_to_Keep'], dtype = 'bool')
#    col = np.array(robjects.r['col_Index'], dtype = 'int')
#    col = col - 1
    
    matrix = ExpMatrix.read_tsv(DataPath, sep = ',')    
#    matrix = matrix.iloc[tokeep] 
    
    truelab = pd.read_csv(LabelsPath, header=0,index_col=None, sep=',')
#    truelab = truelab.iloc[tokeep]
    
    ct_old = ['CD19+ B','CD14+ Monocyte','CD4+/CD45RA+/CD25- Naive T','CD4+/CD45RO+ Memory','CD8+/CD45RA+ Naive Cytotoxic','Dendritic', 'CD56+ NK']
    ct_new = ['B cells','CD14+ monocytes','Naive CD4+ T cells','Memory CD4+ T cells','Naive CD8+ T cells','Dendritic cells','NK cells']
    
    tokeep2 = np.isin(truelab,ct_old)
    truelab = truelab[tokeep2]
    print(len(truelab))
    matrix = matrix.iloc[np.squeeze(tokeep2)]
    
    for i in range(len(ct_old)):
        truelab.iloc[truelab == ct_old[i]] = ct_new[i]
        
    # read the feature file
    if (NumGenes > 0):
        features = pd.read_csv(GeneOrderPath,header=0,index_col=None, sep=',')
        feat_to_use = features.iloc[0:NumGenes,0]
        matrix = matrix.iloc[:,feat_to_use]

    data = ExpMatrix(X = np.transpose(matrix.X), genes = matrix.cells, cells = matrix.genes)
    data.genes.name = 'Genes'
    data.cells.name = 'Cells'
    data.index.name = 'Genes'
    data.columns.name = 'Cells'
    
    clf = CellTypeClassifier.read_pickle(ClassifierPath)
    
    start = tm.time()
    predictions = clf.predict(data)
    runtime = tm.time() - start
    
    np.asarray(predictions)
    
    pred = pd.DataFrame(predictions)
        
    os.chdir(OutputDir)
            
    if (NumGenes == 0):  
        truelab.to_csv("moana_True_Labels.csv", index = False)
        pred.to_csv("moana_Pred_Labels.csv", index = False)
        with open("moana_Total_Time.csv", 'w') as f:
            f.write("%f\n" % runtime)
    else:
        truelab.to_csv("moana_" + str(NumGenes) + "_True_Labels.csv", index = False)
        pred.to_csv("moana_" + str(NumGenes) + "_Pred_Labels.csv", index = False)
        with open("moana_" + str(NumGenes) + "_Total_Time.csv", 'w') as f:
            f.write("%f\n" % runtime)


        
    
    
    
    
    
    
