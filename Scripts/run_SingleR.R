run_SingleR<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir,GeneOrderPath = NULL,NumGenes = NULL){
  "
  run SingleR
  Wrapper script to run SingleR on a benchmark dataset with 5-fold cross validation,
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
  #                               SingleR                                     #
  #############################################################################
  library(SingleR)
  library(Seurat)
  True_Labels_SingleR <- list()
  Pred_Labels_SingleR <- list()
  Total_Time_SingleR <- list()
  Data = t(as.matrix(Data))
  
  for (i in c(1:n_folds)){
    if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
      start_time <- Sys.time()
      singler = SingleR(method = "single", Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]], 
                        Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]], 
                        Labels[Train_Idx[[i]]], numCores = 1)
      end_time <- Sys.time()
    }
    else{
      start_time <- Sys.time()
      singler = SingleR(method = "single", Data[,Test_Idx[[i]]], Data[,Train_Idx[[i]]], Labels[Train_Idx[[i]]], numCores = 1)
      end_time <- Sys.time()
    }
    Total_Time_SingleR[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    True_Labels_SingleR[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_SingleR[i] <- list(as.vector(singler$labels))
  }
  True_Labels_SingleR <- as.vector(unlist(True_Labels_SingleR))
  Pred_Labels_SingleR <- as.vector(unlist(Pred_Labels_SingleR))
  Total_Time_SingleR <- as.vector(unlist(Total_Time_SingleR))
  
  setwd(OutputDir)
  
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    write.csv(True_Labels_SingleR,paste('SingleR_',NumGenes,'_True_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Pred_Labels_SingleR,paste('SingleR_',NumGenes,'_Pred_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Total_Time_SingleR,paste('SingleR_',NumGenes,'_Total_Time.csv', sep = ''),row.names = FALSE)
  }
  else{
    write.csv(True_Labels_SingleR,'SingleR_True_Labels.csv',row.names = FALSE)
    write.csv(Pred_Labels_SingleR,'SingleR_Pred_Labels.csv',row.names = FALSE)
    write.csv(Total_Time_SingleR,'SingleR_Total_Time.csv',row.names = FALSE)
  }
}
