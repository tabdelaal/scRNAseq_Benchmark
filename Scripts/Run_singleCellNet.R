args <- commandArgs(TRUE)

Run_singleCellNet<-function(DataPath,LabelsPath,CV_RDataPath, output_dir){
  Data <- read.csv(DataPath,row.names = 1)
  colnames(Data) <- gsub('_','.',colnames(Data), fixed = TRUE)
  Labels <- as.matrix(read.csv(LabelsPath))
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])
  Data <- Data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]

  #############################################################################
  #                              singleCellNet                                #
  #############################################################################
  library(singleCellNet)
  library(dplyr)
  True_Labels_singleCellNet <- list()
  Pred_Labels_singleCellNet <- list()
  Training_Time_singleCellNet <- list()
  Testing_Time_singleCellNet <- list()
  Data = t(as.matrix(Data))              # deals also with sparse matrix

  for(i in c(1:n_folds)){
    DataTrain <- Data[,Train_Idx[[i]]]
    DataTest <- Data[,Test_Idx[[i]]]

    start_time <- Sys.time()
    cgenes2<-findClassyGenes(DataTrain, data.frame(Annotation = Labels[Train_Idx[[i]]]), "Annotation")
    cgenesA<-cgenes2[['cgenes']]
    grps<-cgenes2[['grps']]
    DataTrain<-as.matrix(DataTrain[cgenesA,])
    xpairs<-ptGetTop(DataTrain, grps, ncores = 1)
    pdTrain<-query_transform(DataTrain[cgenesA, ], xpairs)
    rf<-sc_makeClassifier(pdTrain[xpairs,], genes=xpairs, groups=grps)
    end_time <- Sys.time()
    Training_Time_singleCellNet[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))

    start_time <- Sys.time()
    DataTest<-query_transform(DataTest[cgenesA,], xpairs)
    classRes <-rf_classPredict(rf, DataTest)
    end_time <- Sys.time()
    Testing_Time_singleCellNet[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))

    True_Labels_singleCellNet[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_singleCellNet[i] <- list((rownames(classRes)[apply(classRes,2,which.max)])[1:length(Test_Idx[[i]])])
  }
  True_Labels_singleCellNet <- as.vector(unlist(True_Labels_singleCellNet))
  Pred_Labels_singleCellNet <- as.vector(unlist(Pred_Labels_singleCellNet))
  Training_Time_singleCellNet <- as.vector(unlist(Training_Time_singleCellNet))
  Testing_Time_singleCellNet <- as.vector(unlist(Testing_Time_singleCellNet))
  write.csv(True_Labels_singleCellNet,paste0(output_dir,'/singleCellNet_true.csv'),row.names = FALSE)
  write.csv(Pred_Labels_singleCellNet,paste0(output_dir,'/singleCellNet_pred.csv'),row.names = FALSE)
  write.csv(Training_Time_singleCellNet,paste0(output_dir,'/singleCellNet_training_time.csv'),row.names = FALSE)
  write.csv(Testing_Time_singleCellNet,paste0(output_dir,'/singleCellNet_test_time.csv'),row.names = FALSE)
}

Run_singleCellNet(args[1], args[2], args[3], args[4])
