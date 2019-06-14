withCallingHandlers({
  install.packages("BiocManager", repos="https://cloud.r-project.org/")
  BiocManager::install(version = "devel", ask = FALSE);
  BiocManager::install(c("scater", "MAST"))
  install.packages("devtools", repos="https://cloud.r-project.org/")
  devtools::install_github("satijalab/seurat")
  devtools::install_github("BatadaLab/scID")
},
warning = function(w) stop(w))
