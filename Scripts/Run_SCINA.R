Run_SCINA<-function(DataPath,LabelsPath,GeneSigPath){
  Data <- read.csv(DataPath,row.names = 1)
  Labels <- as.vector(as.matrix(read.csv(LabelsPath)))
  Data <- Data[is.element(Labels,c('CD14+ Monocyte','CD19+ B','CD56+ NK')),]
  Labels <- Labels[is.element(Labels,c('CD14+ Monocyte','CD19+ B','CD56+ NK'))]
  Labels[Labels == 'CD14+ Monocyte'] <- 'CD14_Monocyte'
  Labels[Labels == 'CD19+ B'] <- 'CD19_B'
  Labels[Labels == 'CD56+ NK'] <- 'CD56_NK'
  
  
  #############################################################################
  #                                 SCINA                                     #
  #############################################################################
  library(SCINA)
  Signature_Genes <- preprocess.signatures(GeneSigPath)
  True_Labels_SCINA <- list()
  Pred_Labels_SCINA <- list()
  Total_Time_SCINA <- list()
  
  library(preprocessCore)
  Data = t(as.matrix(Data))
  Data=log(Data+1)
  Data[]=normalize.quantiles(Data)
  
  start_time <- Sys.time()
  results = SCINA(Data, Signature_Genes)
  end_time <- Sys.time()
  
  True_Labels_SCINA <- Labels
  Pred_Labels_SCINA <- results$cell_labels
  Total_Time_SCINA <- as.numeric(difftime(end_time,start_time,units = 'secs'))
  write.csv(True_Labels_SCINA,'True_Labels_SCINA.csv',row.names = FALSE)
  write.csv(Pred_Labels_SCINA,'Pred_Labels_SCINA.csv',row.names = FALSE)
  write.csv(Total_Time_SCINA,'Total_Time_SCINA.csv',row.names = FALSE)
}