args <- commandArgs(TRUE)

run_Garnett_Pretrained <- function(DataPath, LabelsPath, GenesPath, CV_RDataPath, ClassifierPath, OutputDir, Human){
  "
  run Garnett
  Wrapper script to run Garnett on a benchmark dataset with a pretrained classifier,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.

  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  GenesPath : Path to the file with the genenames
  ClassifierPath : Path to the pretrained classifier
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

  # load data, genes, and marker file
  load(CV_RDataPath)

  load(ClassifierPath)

  labels <- as.matrix(read.csv(LabelsPath))
  labels <- labels[Cells_to_Keep]

  mat <- read.table(DataPath, sep = ",")
  data <- mat[-1,-1]
  data <- data[Cells_to_Keep,]
  data <- t(data) #ensure that the genes are rows, and the cells are columns

  barcodes <- mat[-1,1]

  pdata = data.frame(barcodes)
  fdata <- read.table(GenesPath)
  names(fdata) <- 'gene_short_name'
  row.names(fdata) <- fdata$gene_short_name

  row.names(data) <- row.names(fdata)
  colnames(data) <- row.names(pdata)

  pd <- new("AnnotatedDataFrame", data = pdata)
  fd <- new("AnnotatedDataFrame", data = fdata)
  pbmc_cds <- newCellDataSet(as(data, "dgCMatrix"),
                             phenoData = pd,
                             featureData = fd)

  start_time <- Sys.time()

  pbmc_cds <- estimateSizeFactors(pbmc_cds)

  if (Human){
    pbmc_cds <- classify_cells(pbmc_cds, hsPBMC, db = org.Hs.eg.db, cluster_extend = TRUE, cds_gene_id_type = "SYMBOL")
  } else {
    pbmc_cds <- classify_cells(pbmc_cds, mmLung, db = org.Mm.eg.db, cluster_extend = TRUE, cds_gene_id_type = "SYMBOL")
  }

  end_time <- Sys.time()

  test_time <- as.numeric(end_time - start_time)

  write.table(pData(pbmc_cds)$cluster_ext_type,
              file = paste0(OutputDir, "/Garnett_Pretrained_pred.csv"), append = FALSE, quote = TRUE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              qmethod = c("escape", "double"),
              fileEncoding = "")

  write.csv(labels,paste0(OutputDir,"/Garnett_Pretrained_true.csv"), row.names = FALSE)
  write.csv(test_time,paste0(OutputDir,'/Garnett_Pretrained_test_time.csv'),row.names = FALSE)
}

run_Garnett_Pretrained(args[1], args[2], args[3], args[4], args[5], args[6], args[7])
