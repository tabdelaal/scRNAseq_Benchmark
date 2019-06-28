args <- commandArgs(TRUE)

run_Garnett_CV <- function(DataPath, LabelsPath, CV_RDataPath, GenesPath, MarkerPath, OutputDir, Human){
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

  cells <- mat[-1,1]
  cells <- cells[Cells_to_Keep]

  # read the genefile
  fdata <- read.table(GenesPath)
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

  write.csv(true_labels,paste0(OutputDir,'/Garnett_CV_true.csv'),row.names = FALSE)
  write.csv(pred_labels,paste0(OutputDir,'/Garnett_CV_pred.csv'),row.names = FALSE)
  write.csv(train_time,paste0(OutputDir,'/Garnett_CV_training_time.csv'),row.names = FALSE)
  write.csv(test_time,paste0(OutputDir,'/Garnett_CV_test_time.csv'),row.names = FALSE)

}

run_Garnett_CV(args[1], args[2], args[3], args[4], args[5], args[6], args[7])
