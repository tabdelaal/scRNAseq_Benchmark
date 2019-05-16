Run_scPred<-function(DataPath,LabelsPath,CV_RDataPath, GeneOrderPath = NULL, num_of_genes = NULL){
  Data <- read.csv(DataPath,row.names = 1)
  Labels <- as.matrix(read.csv(LabelsPath))
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])
  Data <- Data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]
  if(!is.null(GeneOrderPath) & !is.null (num_of_genes)){
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
    if(!is.null(GeneOrderPath) & !is.null (num_of_genes)){
      sce <- SingleCellExperiment(list(normcounts = Data[as.vector(GenesOrder[c(1:num_of_genes),i])+1,Train_Idx[[i]]]), 
                                  colData = data.frame(cell_type1 = Labels[Train_Idx[[i]]]))
      sce_counts <- normcounts(sce)
      sce_cpm <- apply(sce_counts, 2, function(x) (x/sum(x))*1000000)
      sce_metadata <- as.data.frame(colData(sce))
      
      sce_test <- SingleCellExperiment(list(normcounts = Data[as.vector(GenesOrder[c(1:num_of_genes),i])+1,Test_Idx[[i]]]), 
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
  if(!is.null(GeneOrderPath) & !is.null (num_of_genes)){
    write.csv(True_Labels_scPred,paste('True_Labels_scPred_',num_of_genes,'.csv', sep = ''),row.names = FALSE)
    write.csv(Pred_Labels_scPred,paste('Pred_Labels_scPred_',num_of_genes,'.csv', sep = ''),row.names = FALSE)
    write.csv(Training_Time_scPred,paste('Training_Time_scPred_',num_of_genes,'.csv', sep = ''),row.names = FALSE)
    write.csv(Testing_Time_scPred,paste('Testing_Time_scPred_',num_of_genes,'.csv', sep = ''),row.names = FALSE)
  }
  else{
    write.csv(True_Labels_scPred,'True_Labels_scPred.csv',row.names = FALSE)
    write.csv(Pred_Labels_scPred,'Pred_Labels_scPred.csv',row.names = FALSE)
    write.csv(Training_Time_scPred,'Training_Time_scPred.csv',row.names = FALSE)
    write.csv(Testing_Time_scPred,'Testing_Time_scPred.csv',row.names = FALSE)
  }
}
