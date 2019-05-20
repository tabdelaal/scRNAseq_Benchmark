run_singleCellNet<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir,GeneOrderPath = NULL,NumGenes = NULL){
  "
  run singleCellNet
  Wrapper script to run singleCellNet on a benchmark dataset with 5-fold cross validation,
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
  colnames(Data) <- gsub('_','.',colnames(Data), fixed = TRUE)
  Labels <- as.matrix(read.csv(LabelsPath))
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])
  Data <- Data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    GenesOrder = read.csv(GeneOrderPath)
  }
  
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
    if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
      DataTrain <- Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]]
      DataTest <- Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]]
    }
    else{
      DataTrain <- Data[,Train_Idx[[i]]]
      DataTest <- Data[,Test_Idx[[i]]]
    }
    
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
  
  setwd(OutputDir)
  
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    write.csv(True_Labels_singleCellNet,paste('singleCellNet_',NumGenes,'_True_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Pred_Labels_singleCellNet,paste('singleCellNet_',NumGenes,'_Pred_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Training_Time_singleCellNet,paste('singleCellNet_',NumGenes,'_Training_Time.csv', sep = ''),row.names = FALSE)
    write.csv(Testing_Time_singleCellNet,paste('singleCellNet_',NumGenes,'_Testing_Time.csv', sep = ''),row.names = FALSE)
  }
  else{
    write.csv(True_Labels_singleCellNet,'singleCellNet_True_Labels.csv',row.names = FALSE)
    write.csv(Pred_Labels_singleCellNet,'singleCellNet_Pred_Labels.csv',row.names = FALSE)
    write.csv(Training_Time_singleCellNet,'singleCellNet_Training_Time.csv',row.names = FALSE)
    write.csv(Testing_Time_singleCellNet,'singleCellNet_Testing_Time.csv',row.names = FALSE)
  }
}
