Run_SingleR<-function(DataPath,LabelsPath,CV_RDataPath, GeneOrderPath = NULL, num_of_genes = NULL){
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
  #                               SingleR                                     #
  #############################################################################
  library(SingleR)
  library(Seurat)
  True_Labels_SingleR <- list()
  Pred_Labels_SingleR <- list()
  Total_Time_SingleR <- list()
  Data = t(as.matrix(Data))
  
  for (i in c(1:n_folds)){
    if(!is.null(GeneOrderPath) & !is.null (num_of_genes)){
      start_time <- Sys.time()
      singler = SingleR(method = "single", Data[as.vector(GenesOrder[c(1:num_of_genes),i])+1,Test_Idx[[i]]], 
                        Data[as.vector(GenesOrder[c(1:num_of_genes),i])+1,Train_Idx[[i]]], 
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
  if(!is.null(GeneOrderPath) & !is.null (num_of_genes)){
    write.csv(True_Labels_SingleR,paste('True_Labels_SingleR_',num_of_genes,'.csv', sep = ''),row.names = FALSE)
    write.csv(Pred_Labels_SingleR,paste('Pred_Labels_SingleR_',num_of_genes,'.csv', sep = ''),row.names = FALSE)
    write.csv(Total_Time_SingleR,paste('Total_Time_SingleR_',num_of_genes,'.csv', sep = ''),row.names = FALSE)
  }
  else{
    write.csv(True_Labels_SingleR,'True_Labels_SingleR.csv',row.names = FALSE)
    write.csv(Pred_Labels_SingleR,'Pred_Labels_SingleR.csv',row.names = FALSE)
    write.csv(Total_Time_SingleR,'Total_Time_SingleR.csv',row.names = FALSE)
  }
}
