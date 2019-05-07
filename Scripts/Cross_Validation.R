args <- commandArgs(TRUE)

Cross_Validation <- function(LabelsFilePath, col_Index, output_dir){

  Labels <- as.matrix(read.csv(LabelsFilePath))
  Labels <- as.vector(Labels[,col_Index])

  Removed_classes <- !(table(Labels) > 10)
  Cells_to_Keep <- !(is.element(Labels,names(Removed_classes)[Removed_classes]))
  Labels <- Labels[Cells_to_Keep]

  # Getting training and testing Folds
  library(rBayesianOptimization)
  n_folds = 5
  Folds <- KFold(Labels,nfolds = n_folds, stratified = TRUE)
  Test_Folds <- c(n_folds:1)
  Train_Idx <- list()
  Test_Idx <- list()
  for (i in c(1:length(Folds))){
    Temp_Folds <- Folds
    Temp_Folds[Test_Folds[i]] <- NULL
    Train_Idx[i] <- list(unlist(Temp_Folds))
    Test_Idx[i] <- Folds[Test_Folds[i]]
  }
  remove(Temp_Folds,i,Folds)
  save(n_folds,Train_Idx,Test_Idx,col_Index,Cells_to_Keep,file = paste0(output_dir, '/CV_folds.RData'))
}

Cross_Validation(args[1], as.numeric(args[2]), args[3])
