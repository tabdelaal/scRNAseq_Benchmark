args <- commandArgs(TRUE)

run_scID<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir,GeneOrderPath = NULL,NumGenes = NULL){
  "
  run scID
  Wrapper script to run scID on a benchmark dataset with 5-fold cross validation,
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
  #                                 scID                                      #
  #############################################################################
  library(scID)
  library(Seurat)
  True_Labels_scID <- list()
  Pred_Labels_scID <- list()
  Total_Time_scID <- list()
  Data = t(as.matrix(Data))

  for (i in c(1:n_folds)){
    if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
      Train_Labels <- list(Labels[Train_Idx[[i]]])
      names(Train_Labels[[1]]) <- colnames(Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]])
      start_time <- Sys.time()
      scID_output <- scid_multiclass(Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]],
                                     Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]],
                                     Train_Labels[[1]])
      end_time <- Sys.time()
    }
    else{
      Train_Labels <- list(Labels[Train_Idx[[i]]])
      names(Train_Labels[[1]]) <- colnames(Data[,Train_Idx[[i]]])
      start_time <- Sys.time()
      scID_output <- scid_multiclass(Data[,Test_Idx[[i]]], Data[,Train_Idx[[i]]], Train_Labels[[1]])
      end_time <- Sys.time()
    }
    Total_Time_scID[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))

    True_Labels_scID[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_scID[i] <- list(as.vector(scID_output$labels))
  }
  True_Labels_scID <- as.vector(unlist(True_Labels_scID))
  Pred_Labels_scID <- as.vector(unlist(Pred_Labels_scID))
  Total_Time_scID <- as.vector(unlist(Total_Time_scID))

  setwd(OutputDir)

  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    write.csv(True_Labels_scID,paste('scID_',NumGenes,'_true', sep = ''),row.names = FALSE)
    write.csv(Pred_Labels_scID,paste('scID_',NumGenes,'_pred', sep = ''),row.names = FALSE)
    write.csv(Total_Time_scID,paste('scID_',NumGenes,'_total_time.csv', sep = ''),row.names = FALSE)
  }
  else{
    write.csv(Pred_Labels_scID,'scID_pred.csv',row.names = FALSE)
    write.csv(True_Labels_scID,'scID_true.csv',row.names = FALSE)
    write.csv(Total_Time_scID,'scID_total_time.csv',row.names = FALSE)
  }
}
if (args[6] == "0") {
  run_scID(args[1], args[2], args[3], args[4])
} else {
  run_scID(args[1], args[2], args[3], args[4], args[5], as.numeric(args[6]))
}
