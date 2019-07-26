selectMarkerGenes_DCS <- function(DataPath, LabelsPath, CV_RDataPath, OutputDir, Normalize = FALSE, LogTransform = FALSE, NumGenes = 20){
  "
  Wrapper script to select the marker genes for DigitalCellSorter on a benchmark dataset with 5-fold cross validation,
  outputs excel file in the format that DigitalCellSorter can read as input.
  
  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported files
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
  
  library(xlsx)
  
  Data <- read.csv(DataPath,row.names = 1)
  Labels <- as.matrix(read.csv(LabelsPath))
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])
  Data <- Data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]
  
  CellTypes <- unique(Labels)
  CellTypes <- data.frame(CellTypes)
  names(CellTypes) <- 'CellType'
  CellTypes$CellTypeGrouped <- unique(Labels)
  
  setwd(OutputDir)
  
  Data <- t(Data)
    
  for (i in c(1:n_folds)){
    Train <- Data[,Train_Idx[[i]]]
    Lab_Train <- Labels[Train_Idx[[i]]]
    
    Markers <- DEgenesMAST(Train, Lab_Train, Normalize = Normalize, LogTransform = LogTransform)
    
    File = paste('Markers_',i,'.RData',sep='')
    save(Markers, file = File)
    
      if(NumGenes[k] < 20){
        Markers <- Markers[1:NumGenes[k],]
      }

      Unique_Markers <- unique(na.omit(as.vector(Markers)))
      Res <- matrix(data = 0, nrow = length(Unique_Markers), ncol = ncol(Markers))
      Res <- data.frame(Res, row.names = Unique_Markers)
      names(Res) <- colnames(Markers)

      for (j in c(1:ncol(Markers))){
        Res[is.element(row.names(Res),Markers[,j]),j] = 1
      }
      
      Res$Marker <- rownames(Res)
      Res <- Res[,c(ncol(Res),1:(ncol(Res)-1))]

      Excel_File <- paste('DCS_Marker_Genes_',NumGenes[k],'_Fold_',i,'.xlsx',sep = '')

      # write MarkerCellType sheet
      write.xlsx(Res, Excel_File, sheetName = 'MarkerCellType', col.names = TRUE, row.names = FALSE)

      # # append CellTypesGrouped sheet
      write.xlsx(CellTypes, Excel_File, sheetName = 'CellTypesGrouped', col.names = TRUE, row.names = FALSE, append = TRUE)

  }
  
}