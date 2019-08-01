withCallingHandlers({
  install.packages("BiocManager", repos="https://cloud.r-project.org/")
  BiocManager::install(c("monocle", "DelayedArray", "DelayedMatrixStats",
                       "org.Hs.eg.db", "org.Mm.eg.db"))
  install.packages("devtools", repos="https://cloud.r-project.org/")
  devtools::install_github("cole-trapnell-lab/garnett", ref="0.1.4")
},
warning = function(w) stop(w))
