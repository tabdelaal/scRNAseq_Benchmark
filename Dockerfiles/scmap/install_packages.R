withCallingHandlers({
  install.packages("BiocManager", repos="https://cloud.r-project.org/")
  BiocManager::install(ask = FALSE)
  BiocManager::install("SingleCellExperiment")
  install.packages("devtools", repos="https://cloud.r-project.org/")
  devtools::install_github("hemberg-lab/scmap", ref="v1.1.5")
},
warning = function(w) stop(w))
