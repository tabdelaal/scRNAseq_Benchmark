Run_CHETAH<-function(DataPath,LabelsPath,CV_RDataPath, GeneOrderPath = NULL, num_of_genes = NULL){
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
  #                                CHETAH                                     #
  #############################################################################
  library(CHETAH)
  library(SingleCellExperiment)
  True_Labels_CHETAH <- list()
  Pred_Labels_CHETAH <- list()
  Total_Time_CHETAH <- list()
  Data = t(as.matrix(Data))
  
  for (i in c(1:n_folds)){
    if(!is.null(GeneOrderPath) & !is.null (num_of_genes)){
      sce <- SingleCellExperiment(assays = list(counts = Data[as.vector(GenesOrder[c(1:num_of_genes),i])+1,Train_Idx[[i]]]), 
                                  colData = data.frame(celltypes = Labels[Train_Idx[[i]]]))
      
      sce_test <- SingleCellExperiment(assays = list(counts = Data[as.vector(GenesOrder[c(1:num_of_genes),i])+1,Test_Idx[[i]]]), 
                                       colData = data.frame(celltypes = Labels[Test_Idx[[i]]]))
      start_time <- Sys.time()
      sce_test <- CHETAHclassifier(input = sce_test, ref_cells = sce, n_genes = num_of_genes)
      end_time <- Sys.time()
    }
    else{
      sce <- SingleCellExperiment(assays = list(counts = Data[,Train_Idx[[i]]]), 
                                  colData = data.frame(celltypes = Labels[Train_Idx[[i]]]))
      
      sce_test <- SingleCellExperiment(assays = list(counts = Data[,Test_Idx[[i]]]), 
                                       colData = data.frame(celltypes = Labels[Test_Idx[[i]]]))
      start_time <- Sys.time()
      sce_test <- CHETAHclassifier(input = sce_test, ref_cells = sce)
      end_time <- Sys.time()
    }
    
    Total_Time_CHETAH[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    True_Labels_CHETAH[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_CHETAH[i] <- list(sce_test$celltype_CHETAH)
  }
  True_Labels_CHETAH <- as.vector(unlist(True_Labels_CHETAH))
  Pred_Labels_CHETAH <- as.vector(unlist(Pred_Labels_CHETAH))
  Total_Time_CHETAH <- as.vector(unlist(Total_Time_CHETAH))
  
  if (!is.null(GeneOrderPath) & !is.null (num_of_genes)){
    write.csv(True_Labels_CHETAH,paste('True_Labels_CHETAH_',num_of_genes,'.csv', sep = ''),row.names = FALSE)
    write.csv(Pred_Labels_CHETAH,paste('Pred_Labels_CHETAH_',num_of_genes,'.csv', sep = ''),row.names = FALSE)
    write.csv(Total_Time_CHETAH,paste('Total_Time_CHETAH_',num_of_genes,'.csv', sep = ''),row.names = FALSE)
  }
  else{
    write.csv(True_Labels_CHETAH,'True_Labels_CHETAH.csv',row.names = FALSE)
    write.csv(Pred_Labels_CHETAH,'Pred_Labels_CHETAH.csv',row.names = FALSE)
    write.csv(Total_Time_CHETAH,'Total_Time_CHETAH.csv',row.names = FALSE)
  }
}
