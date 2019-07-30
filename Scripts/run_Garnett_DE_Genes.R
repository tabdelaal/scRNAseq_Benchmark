run_Garnett_DE_Genes <- function(DataPath, LabelsPath, CV_RDataPath, GenesPath, OutputDir, Human = TRUE, NumGenes = 20, Normalize = FALSE, LogTransform = FALSE){
  "
  run Garnett
  Wrapper script to run Garnett on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  GenesPath : Path to the file with the genenames
  MarkerPath : Path to the file with marker genes
  OutputDir : Output directory defining the path of the exported file.
  Human : boolean indicating whether the dataset is human (TRUE) or mouse (FALSE)
  "
  
  DEgenesMAST <- function(Data, Labels, Normalize = FALSE, LogTransform = FALSE){
    # Data: genesXcells (rows = genes, columns = cells)
    
    library(Seurat)
    
    if(Normalize)
    {
      Data <- apply(Data, 2, function(x) (x/sum(x))*1000000)
    }
    
    if(LogTransform)
    {
      Data <- log(Data+1, base = 2)
    }
    SeuObj <- CreateSeuratObject(raw.data = Data, project = "DEgenes")
    SeuObj <- SetIdent(SeuObj, ident.use = Labels)
    DEgenes <- FindAllMarkers(SeuObj, test.use = "MAST")
    Markers <- matrix(nrow = 20,ncol = length(unique(Labels)))
    colnames(Markers) <- unique(Labels)
    for (i in unique(Labels)){
      i
      TempList <- DEgenes$gene[((DEgenes$cluster == i) & (DEgenes$avg_logFC > 0))]
      MarkerGenes <- DEgenes$p_val_adj[DEgenes$cluster == i]
      print(MarkerGenes[1:20])
      if (length(TempList) >= 20){
        Markers[,i] <- TempList[1:20]
      }
      else{
        if(length(TempList) > 0){
          Markers[c(1:length(TempList)),i] <- TempList
        }
      }
    }
    return(Markers)
  }
  
  # load needed libraries
  library(garnett)
  if (Human) {
    library(org.Hs.eg.db)
  } else {
    library(org.Mm.eg.db)
  }
  
  # load the CVFile
  load(CV_RDataPath)
  
  # read the labels
  labels <- as.matrix(read.csv(LabelsPath))
  labels <- as.vector(labels[,col_Index])
  labels <- labels[Cells_to_Keep]
  
  # read the data
  mat <- read.table(DataPath, sep = ",")
  data <- mat[-1,-1]
  data <- data[Cells_to_Keep,]
  data <- t(data) #ensure that the genes are rows, and the cells are columns

  data_markergenes <- read.csv(DataPath,row.names = 1)
  data_markergenes <- data_markergenes[Cells_to_Keep,]
  data_markergenes <- t(data_markergenes)
  
  cells <- mat[-1,1]
  cells <- cells[Cells_to_Keep]
  
  # read the genefile 
  fdata <- read.table(GenesPath, header = FALSE)
  names(fdata) <- 'gene_short_name'
  row.names(fdata) <- fdata$gene_short_name
  fd <- new("AnnotatedDataFrame", data = fdata)
  
  true_labels <- list()
  pred_labels <- list()
  train_time <- list()
  test_time <- list()
  
  setwd(OutputDir)
  
  for (i in c(1:n_folds)){
    lab_train = labels[Train_Idx[[i]]]
    lab_test = labels[Test_Idx[[i]]]
    
    train = data[,Train_Idx[[i]]]
    test = data[,Test_Idx[[i]]]
    
    mat_train = data_markergenes[,Train_Idx[[i]]]
    
    ## Select the marker genes based on the training set and write the marker file
    Markers = DEgenesMAST(mat_train, lab_train, Normalize = Normalize, LogTransform = LogTransform)
     
    if(NumGenes < 20){
      Markers <- Markers[1:NumGenes,]
    }
     
    MarkerPath = paste('Garnett_Marker_Genes_' , NumGenes , '_Fold_' , i , '.txt', sep = "")
     
    sink(MarkerPath)
     
    for (j in c(1:length(unique(labels)))){
      cat('>')
      cat(colnames(Markers)[j])
      cat('\n')
      
      cat('expressed: ')
      
      for (k in c(1:nrow(Markers))){
        cat(Markers[k,j])
        
        if (k != nrow(Markers)){
          cat(', ')
        }
      }
      
      cat('\n')
      cat('\n')
      
    }
    
    sink()
    
    cells_train = cells[Train_Idx[[i]]]
    cells_test = cells[Test_Idx[[i]]]
    
    pdata_train = data.frame(cells_train)
    pdata_test = data.frame(cells_test)
    
    row.names(train) <- row.names(fdata)
    row.names(test) <- row.names(fdata)
    colnames(train) <- row.names(pdata_train)
    colnames(test) <- row.names(pdata_test)
    
    pd_train <- new("AnnotatedDataFrame", data = pdata_train)
    pd_test <- new("AnnotatedDataFrame", data = pdata_test)
    
    pbmc_cds_train <- newCellDataSet(as(train, "dgCMatrix"), phenoData = pd_train, featureData = fd)
    pbmc_cds_test <- newCellDataSet(as(test, "dgCMatrix"), phenoData = pd_test, featureData = fd)
    
    pbmc_cds_train <- estimateSizeFactors(pbmc_cds_train)
    pbmc_cds_test <- estimateSizeFactors(pbmc_cds_test)
    
    # training
    start_train <- Sys.time()
    
    if (Human){
      pbmc_classifier <- train_cell_classifier(cds = pbmc_cds_train, 
                                               marker_file = MarkerPath,
                                               db=org.Hs.eg.db,
                                               cds_gene_id_type = "SYMBOL",
                                               num_unknown = 50,
                                               marker_file_gene_id_type = "SYMBOL")
    } else {
      pbmc_classifier <- train_cell_classifier(cds = pbmc_cds_train, 
                                               marker_file = MarkerPath,
                                               db=org.Mm.eg.db,
                                               cds_gene_id_type = "SYMBOL",
                                               num_unknown = 50,
                                               marker_file_gene_id_type = "SYMBOL")
      
    }
    end_train <- Sys.time()
    train_time[i] <- as.numeric(end_train - start_train)
    
    # testing
    start_test <- Sys.time()
    
    if (Human) {
      pbmc_cds_test <- classify_cells(pbmc_cds_test, 
                                      pbmc_classifier, 
                                      db = org.Hs.eg.db, 
                                      cluster_extend = TRUE,
                                      cds_gene_id_type = "SYMBOL")
    } else {
      pbmc_cds_test <- classify_cells(pbmc_cds_test, 
                                      pbmc_classifier, 
                                      db = org.Mm.eg.db, 
                                      cluster_extend = TRUE,
                                      cds_gene_id_type = "SYMBOL")
    }
    end_test <- Sys.time()
    test_time[i] <- as.numeric(end_test - start_test)
    
    true_labels[i] <- list(lab_test)
    pred_labels[i] <- list(pData(pbmc_cds_test)$cluster_ext_type)
    
    
  }
  
  true_labels <- as.vector(unlist(true_labels))
  pred_labels <- as.vector(unlist(pred_labels))
  train_time <- as.vector(unlist(train_time))
  test_time <- as.vector(unlist(test_time))
  
  write.csv(train_time,paste('Garnett_',NumGenes,'_Testing_Time.csv',sep = ""),row.names = FALSE)
  write.csv(test_time,paste('Garnett_',NumGenes,'_Training_Time.csv',sep = ""),row.names = FALSE)
  write.csv(true_labels,paste('Garnett_',NumGenes,'_True_Labels.csv',sep = ""), row.names = FALSE)
  write.csv(pred_labels,paste('Garnett_',NumGenes,'_Pred_Labels.csv',sep = ""), row.names = FALSE)
  
  
}