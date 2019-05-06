run_Garnett <- function(input_dir, output_dir, datafile, classifier, genefile){
  "
  Run Garnett

  Wrapper script to run Garnett with a pretrained classifier
  
  Parameters
  ----------
  input_dir : directory of the input files
  output_dir : directory of the output files
  datafile : name of the data file
  classifier : name of the file containing the pretrained classifier
  genefile : name of the file containing the genes
  " 
  # load needed libraries
  library(garnett)
  library(org.Hs.eg.db)
  
  # load data, genes, and marker file
  setwd(input_dir)
  load(classifier)
  
  mat <- read.table(datafile, sep = ",")
  data <- mat[-1,-1]
  data <- t(data) #ensure that the genes are rows, and the cells are columns
  
  barcodes <- mat[-1,1]
  
  pdata = data.frame(barcodes)
  fdata <- read.table(genefile)
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
  pbmc_cds <- classify_cells(pbmc_cds, hsPBMC, db = org.Hs.eg.db, cluster_extend = TRUE, cds_gene_id_type = "SYMBOL")
  
  end_time <- Sys.time()
  
  test_time <- as.numeric(end_time - start_time)
  
  setwd(output_dir)
  
  write.table(pData(pbmc_cds)$cluster_ext_type, file = "Garnett_pred.csv", append = FALSE, quote = TRUE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              qmethod = c("escape", "double"),
              fileEncoding = "")
  
  write.csv(test_time,'Garnett_test_time.csv',row.names = FALSE)
  
  
  
}