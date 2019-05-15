Run_CaSTLe<-function(DataPath,LabelsPath,CV_RDataPath, GeneOrderPath = NULL, num_of_genes = NULL){
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
  #                                CaSTLe                                     #
  #############################################################################
  library(igraph)
  library(xgboost)
  True_Labels_Castle <- list()
  Pred_Labels_Castle <- list()
  Training_Time_Castle <- list()
  Testing_Time_Castle <- list()
  
  BREAKS=c(-1, 0, 1, 6, Inf)
  nFeatures = 100
  
  for(i in c(1:n_folds)){
    # 1. Load datasets
    if(!is.null(GeneOrderPath) & !is.null (num_of_genes)){
      ds1 = Data[Train_Idx[[i]],as.vector(GenesOrder[c(1:num_of_genes),i])+1]
      ds2 = Data[Test_Idx[[i]],as.vector(GenesOrder[c(1:num_of_genes),i])+1]
    }
    else{
      ds1 = Data[Train_Idx[[i]],]
      ds2 = Data[Test_Idx[[i]],]
    }
    
    sourceCellTypes = as.factor(Labels[Train_Idx[[i]]])
    targetCellTypes = as.factor(Labels[Test_Idx[[i]]])
    
    start_time <- Sys.time()
    # 2. Unify sets, excluding low expressed genes
    source_n_cells_counts = apply(ds1, 2, function(x) { sum(x > 0) } )
    target_n_cells_counts = apply(ds2, 2, function(x) { sum(x > 0) } )
    common_genes = intersect( colnames(ds1)[source_n_cells_counts>10], 
                              colnames(ds2)[target_n_cells_counts>10])
    remove(source_n_cells_counts, target_n_cells_counts)
    ds1 = ds1[, colnames(ds1) %in% common_genes]
    ds2 = ds2[, colnames(ds2) %in% common_genes]
    ds = rbind(ds1[,common_genes], ds2[,common_genes])
    isSource = c(rep(TRUE,nrow(ds1)), rep(FALSE,nrow(ds2)))
    remove(ds1, ds2)
    
    # 3. Highest mean in both source and target
    topFeaturesAvg = colnames(ds)[order(apply(ds, 2, mean), decreasing = T)]
    end_time <- Sys.time()
    Training_Time_Castle[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    start_time <- Sys.time()
    # for each cell - what is the most probable classification?
    L = length(levels(sourceCellTypes))
    targetClassification = as.data.frame(matrix(rep(0,L*sum(!isSource)), nrow=L), row.names = levels(sourceCellTypes))
    
    for (cellType in levels(sourceCellTypes)) {
      
      inSourceCellType = as.factor(ifelse(sourceCellTypes == cellType, cellType, paste0("NOT",cellType)))
      
      # 4. Highest mutual information in source
      topFeaturesMi = names(sort(apply(ds[isSource,],2,function(x) { compare(cut(x,breaks=BREAKS),inSourceCellType,method = "nmi") }), decreasing = T))
      
      # 5. Top n genes that appear in both mi and avg
      selectedFeatures = union(head(topFeaturesAvg, nFeatures) , head(topFeaturesMi, nFeatures) )
      
      # 6. remove correlated features
      tmp = cor(ds[,selectedFeatures], method = "pearson")
      tmp[!lower.tri(tmp)] = 0
      selectedFeatures = selectedFeatures[apply(tmp,2,function(x) any(x < 0.9))]
      remove(tmp)
      
      # 7,8. Convert data from continous to binned dummy vars
      # break datasets to bins
      dsBins = apply(ds[, selectedFeatures], 2, cut, breaks= BREAKS)
      # use only bins with more than one value
      nUniq = apply(dsBins, 2, function(x) { length(unique(x)) })
      # convert to dummy vars
      ds0 = model.matrix(~ . , as.data.frame(dsBins[,nUniq>1]))
      remove(dsBins, nUniq)
      
      cat(paste0("<h2>Classifier for ",cellType,"</h2>"))
      
      inTypeSource = sourceCellTypes == cellType
      # 9. Classify
      xg=xgboost(data=ds0[isSource,] , 
                 label=inTypeSource,
                 objective="binary:logistic", 
                 eta=0.7 , nthread=1, nround=20, verbose=0,
                 gamma=0.001, max_depth=5, min_child_weight=10)
      
      # 10. Predict
      inTypeProb = predict(xg, ds0[!isSource, ])
      
      targetClassification[cellType,] = inTypeProb
    }
    end_time <- Sys.time()
    Testing_Time_Castle[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    True_Labels_Castle[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_Castle[i] <- list(rownames(targetClassification)[apply(targetClassification,2,which.max)])
  }
  True_Labels_Castle <- as.vector(unlist(True_Labels_Castle))
  Pred_Labels_Castle <- as.vector(unlist(Pred_Labels_Castle))
  Training_Time_Castle <- as.vector(unlist(Training_Time_Castle))
  Testing_Time_Castle <- as.vector(unlist(Testing_Time_Castle))
  
  if(!is.null(GeneOrderPath) & !is.null (num_of_genes)){
    write.csv(True_Labels_Castle,paste('True_Labels_Castle_',num_of_genes,'.csv', sep = ''),row.names = FALSE)
    write.csv(Pred_Labels_Castle,paste('Pred_Labels_Castle_',num_of_genes,'.csv', sep = ''),row.names = FALSE)
    write.csv(Training_Time_Castle,paste('Training_Time_Castle_',num_of_genes,'.csv', sep = ''),row.names = FALSE)
    write.csv(Testing_Time_Castle,paste('Testing_Time_Castle_',num_of_genes,'.csv', sep = ''),row.names = FALSE)
  }
  else{
    write.csv(True_Labels_Castle,'True_Labels_CaSTLe.csv',row.names = FALSE)
    write.csv(Pred_Labels_Castle,'Pred_Labels_CaSTLe.csv',row.names = FALSE)
    write.csv(Training_Time_Castle,'Training_Time_CaSTLe.csv',row.names = FALSE)
    write.csv(Testing_Time_Castle,'Testing_Time_CaSTLe.csv',row.names = FALSE)
  }
}
