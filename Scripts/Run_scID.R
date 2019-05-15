Run_scID<-function(DataPath,LabelsPath,CV_RDataPath, GeneOrderPath = NULL, num_of_genes = NULL){
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
  #                                 scID                                      #
  #############################################################################
  library(scID)
  library(Seurat)
  True_Labels_scID <- list()
  Pred_Labels_scID <- list()
  Total_Time_scID <- list()
  Data = t(as.matrix(Data))
  
  for (i in c(1:n_folds)){
    if(!is.null(GeneOrderPath) & !is.null (num_of_genes)){
      Train_Labels <- list(Labels[Train_Idx[[i]]])
      names(Train_Labels[[1]]) <- colnames(Data[as.vector(GenesOrder[c(1:num_of_genes),i])+1,Train_Idx[[i]]])
      start_time <- Sys.time()
      scID_output <- scid_multiclass(Data[as.vector(GenesOrder[c(1:num_of_genes),i])+1,Test_Idx[[i]]], 
                                     Data[as.vector(GenesOrder[c(1:num_of_genes),i])+1,Train_Idx[[i]]], 
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
  if(!is.null(GeneOrderPath) & !is.null (num_of_genes)){
    write.csv(True_Labels_scID,paste('True_Labels_scID_',num_of_genes,'.csv', sep = ''),row.names = FALSE)
    write.csv(Pred_Labels_scID,paste('Pred_Labels_scID_',num_of_genes,'.csv', sep = ''),row.names = FALSE)
    write.csv(Total_Time_scID,paste('Total_Time_scID_',num_of_genes,'.csv', sep = ''),row.names = FALSE)
  }
  else{
    write.csv(True_Labels_scID,'True_Labels_scID.csv',row.names = FALSE)
    write.csv(Pred_Labels_scID,'Pred_Labels_scID.csv',row.names = FALSE)
    write.csv(Total_Time_scID,'Total_Time_scID.csv',row.names = FALSE)
  }
}
