run_Garnett_CV <- function(input_dir, output_dir, datafile, labelfile, markerfile, genefile, CVFile){
  "
  Run Garnett
  
  Wrapper script to run Garnett with using crossvalidation
  
  Parameters
  ----------
  input_dir : directory of the input files
  output_dir : directory of the output files
  datafile : name of the data file
  labelfile : name of the file containing the labels
  markerfile : name of the marker file
  genefile : name of the file containing the genes
  CVfile : name of the file containing cross validation indices
  " 
  setwd(input_dir)
  
  # load needed libraries
  library(garnett)
  library(org.Hs.eg.db)
  
  # load the CVFile
  load(CVFile)
  
  # read the labels
  labels <- as.matrix(read.csv(labelfile))
  labels <- as.vector(labels[,col_Index])
  labels <- labels[Cells_to_Keep]
  
  # read the data
  mat <- read.table(datafile, sep = ",")
  data <- mat[-1,-1]
  data <- data[Cells_to_Keep,]
  data <- t(data) #ensure that the genes are rows, and the cells are columns
  
  cells <- mat[-1,1]
  cells <- cells[Cells_to_Keep]
  
  # read the genefile 
  fdata <- read.table(genefile)
  names(fdata) <- 'gene_short_name'
  row.names(fdata) <- fdata$gene_short_name
  fd <- new("AnnotatedDataFrame", data = fdata)
  
  true_labels <- list()
  pred_labels <- list()
  train_time <- list()
  test_time <- list()
  
  for (i in c(1:n_folds)){
    lab_train = labels[Train_Idx[[i]]]
    lab_test = labels[Test_Idx[[i]]]
    
    train = data[,Train_Idx[[i]]]
    test = data[,Test_Idx[[i]]]
    
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
    pbmc_classifier <- train_cell_classifier(cds = pbmc_cds_train, 
                                             marker_file = markerfile,
                                             db=org.Hs.eg.db,
                                             cds_gene_id_type = "SYMBOL",
                                             num_unknown = 50,
                                             marker_file_gene_id_type = "SYMBOL")
    end_train <- Sys.time()
    train_time[i] <- as.numeric(end_train - start_train)
    
    # testing
    start_test <- Sys.time()
    pbmc_cds_test <- classify_cells(pbmc_cds_test, 
                                    pbmc_classifier, 
                                    db = org.Hs.eg.db, 
                                    cluster_extend = TRUE,
                                    cds_gene_id_type = "SYMBOL")
    end_test <- Sys.time()
    test_time[i] <- as.numeric(end_test - start_test)
    
    true_labels[i] <- list(lab_test)
    pred_labels[i] <- list(pData(pbmc_cds_test)$cluster_ext_type)
    
    
  }
  
  true_labels <- as.vector(unlist(true_labels))
  pred_labels <- as.vector(unlist(pred_labels))
  train_time <- as.vector(unlist(train_time))
  test_time <- as.vector(unlist(test_time))

  setwd(output_dir)
  
  write.csv(train_time,'Garnett_test_time.csv',row.names = FALSE)
  write.csv(test_time,'Garnett_training_time.csv',row.names = FALSE)
  write.csv(true_labels, 'Garnett_true.csv', row.names = FALSE)
  write.csv(pred_labels, 'Garnett_pred.csv', row.names = FALSE)
  
  
}