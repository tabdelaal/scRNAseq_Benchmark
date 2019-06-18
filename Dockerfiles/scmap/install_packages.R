withCallingHandlers({
  install.packages("BiocManager", repos="https://cloud.r-project.org/")
  BiocManager::install(ask = FALSE)
  BiocManager::install("SingleCellExperiment")
  install.packages("devtools", repos="https://cloud.r-project.org/")
  devtools::install_github("hemberg-lab/scmap")
},
warning = function(w) stop(w))
