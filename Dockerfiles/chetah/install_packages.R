withCallingHandlers({
  install.packages("devtools", repos="https://cloud.r-project.org/")
  install.packages("BiocManager", repos="https://cloud.r-project.org/")
  BiocManager::install(c("bioDist", "ggplot2", "gplots", "cowplot",
                         "dendextend", "corrplot", "reshape2", "plotly"))
  devtools::install_github("jdekanter/CHETAH", ref="b777e6f671bff3c434842adb655869a52bc9e368")
},
warning = function(w) stop(w))
