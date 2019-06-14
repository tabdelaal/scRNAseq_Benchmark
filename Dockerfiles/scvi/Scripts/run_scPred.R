run_scPred<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir,GeneOrderPath = NULL,NumGenes = NULL){
  "
  run scPred
  Wrapper script to run scPred on a benchmark dataset with 5-fold cross validation,
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
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  
  Data <- read.csv(DataPath,row.names = 1)
  Labels <- as.matrix(read.csv(LabelsPath))
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])
  Data <- Data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    GenesOrder = read.csv(GeneOrderPath)
  }
  
  #############################################################################
  #                                scPred                                     #
  #############################################################################
  library(scPred)
  library(tidyverse)
  library(SingleCellExperiment)
  True_Labels_scPred <- list()
  Pred_Labels_scPred <- list()
  Training_Time_scPred <- list()
  Testing_Time_scPred <- list()
  Data = t(as.matrix(Data))
  
  for (i in c(1:n_folds)){
    if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
      sce <- SingleCellExperiment(list(normcounts = Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]]), 
                                  colData = data.frame(cell_type1 = Labels[Train_Idx[[i]]]))
      sce_counts <- normcounts(sce)
      sce_cpm <- apply(sce_counts, 2, function(x) (x/sum(x))*1000000)
      sce_metadata <- as.data.frame(colData(sce))
      
      sce_test <- SingleCellExperiment(list(normcounts = Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]]), 
                                       colData = data.frame(cell_type1 = Labels[Test_Idx[[i]]]))
      sce_counts_test <- normcounts(sce_test)
      sce_cpm_test <- apply(sce_counts_test, 2, function(x) (x/sum(x))*1000000)
      sce_metadata_test <- as.data.frame(colData(sce_test))
    }
    else{
      sce <- SingleCellExperiment(list(normcounts = Data[,Train_Idx[[i]]]), 
                                  colData = data.frame(cell_type1 = Labels[Train_Idx[[i]]]))
      sce_counts <- normcounts(sce)
      sce_cpm <- apply(sce_counts, 2, function(x) (x/sum(x))*1000000)
      sce_metadata <- as.data.frame(colData(sce))
      
      sce_test <- SingleCellExperiment(list(normcounts = Data[,Test_Idx[[i]]]), 
                                       colData = data.frame(cell_type1 = Labels[Test_Idx[[i]]]))
      sce_counts_test <- normcounts(sce_test)
      sce_cpm_test <- apply(sce_counts_test, 2, function(x) (x/sum(x))*1000000)
      sce_metadata_test <- as.data.frame(colData(sce_test))
    }
    
    
    # scPred Training    
    start_time <- Sys.time()
    set.seed(1234)
    scp <- eigenDecompose(sce_cpm)
    scPred::metadata(scp) <- sce_metadata
    scp <- getFeatureSpace(scp, pVar = 'cell_type1')
    # plotEigen(scp, group = 'cell_type1')
    scp <- trainModel(scp)
    # plotTrainProbs(scp)
    end_time <- Sys.time()
    Training_Time_scPred[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    # scPred Prediction
    start_time <- Sys.time()
    scp <- scPredict(scp,newData = sce_cpm_test)
    end_time <- Sys.time()
    Testing_Time_scPred[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    True_Labels_scPred[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_scPred[i] <- list(getPredictions(scp)$predClass)
  }
  True_Labels_scPred <- as.vector(unlist(True_Labels_scPred))
  Pred_Labels_scPred <- as.vector(unlist(Pred_Labels_scPred))
  Training_Time_scPred <- as.vector(unlist(Training_Time_scPred))
  Testing_Time_scPred <- as.vector(unlist(Testing_Time_scPred))
  
  setwd(OutputDir)
  
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    write.csv(True_Labels_scPred,paste('scPred_',NumGenes,'_True_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Pred_Labels_scPred,paste('scPred_',NumGenes,'_Pred_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Training_Time_scPred,paste('scPred_',NumGenes,'_Training_Time.csv', sep = ''),row.names = FALSE)
    write.csv(Testing_Time_scPred,paste('scPred_',NumGenes,'_Testing_Time.csv', sep = ''),row.names = FALSE)
  }
  else{
    write.csv(True_Labels_scPred,'scPred_True_Labels.csv',row.names = FALSE)
    write.csv(Pred_Labels_scPred,'scPred_Pred_Labels.csv',row.names = FALSE)
    write.csv(Training_Time_scPred,'scPred_Training_Time.csv',row.names = FALSE)
    write.csv(Testing_Time_scPred,'scPred_Testing_Time.csv',row.names = FALSE)
  }
}
