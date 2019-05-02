Run_SingleR<-function(DataPath,LabelsPath,CV_RDataPath){
  Data <- read.csv(DataPath,row.names = 1)
  Labels <- as.matrix(read.csv(LabelsPath))
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])
  Data <- Data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]
  
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
    start_time <- Sys.time()
    singler = SingleR(method = "single", Data[,Test_Idx[[i]]], Data[,Train_Idx[[i]]], Labels[Train_Idx[[i]]], numCores = 1)
    end_time <- Sys.time()
    Total_Time_SingleR[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    True_Labels_SingleR[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_SingleR[i] <- list(as.vector(singler$labels))
  }
  True_Labels_SingleR <- as.vector(unlist(True_Labels_SingleR))
  Pred_Labels_SingleR <- as.vector(unlist(Pred_Labels_SingleR))
  Total_Time_SingleR <- as.vector(unlist(Total_Time_SingleR))
  write.csv(True_Labels_SingleR,'True_Labels_SingleR.csv',row.names = FALSE)
  write.csv(Pred_Labels_SingleR,'Pred_Labels_SingleR.csv',row.names = FALSE)
  write.csv(Total_Time_SingleR,'Total_Time_SingleR.csv',row.names = FALSE)
}