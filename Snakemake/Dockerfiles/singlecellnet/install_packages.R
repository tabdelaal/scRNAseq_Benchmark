withCallingHandlers({
  install.packages("devtools", repos="https://cloud.r-project.org/")
  install.packages("BiocManager", repos="https://cloud.r-project.org/")
  BiocManager::install("fgsea")
  devtools::install_github("thomasp85/patchwork", ref="fd7958bae3e7a1e30237c751952e412a0a1d1242")
  devtools::install_github("pcahan1/singleCellNet", ref="4279a68112743b783cc82628421dd703261ec117")
},
warning = function(w) stop(w))
