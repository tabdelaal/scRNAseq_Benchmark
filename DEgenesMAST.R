DEgenesMAST <- function(Data, Labels, Normalize = FALSE, LogTransform = FALSE){
  # This functions applies a differential expression test to the data using one vs all
  # The training data should be used a an input
  # The output is a matrix with marker genes where the columns are the cell populations and the rows are the top20 marker genes
  # This output can be rewritten to the format of the prior-knowledge-supervised classifiers and afterwards be used to classify the test set.
  
  # Data: genes X cells (rows = genes, columns = cells)
  # Labels: labels of the data
  # Normalize: the input for MAST should be cpm normalized data, 
  #            if the data is not normalized yet, this should be set to TRUE
  # LogTransform: the input for MAST should be logtransformed,
  #            if the data is not logtransformed yet, this should be set to TRUE
  
  
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
