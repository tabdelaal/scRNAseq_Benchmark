withCallingHandlers({
  install.packages("devtools", repos="https://cloud.r-project.org/")
  install.packages("Seurat", repos="https://cloud.r-project.org/")
  devtools::install_github("dviraran/SingleR", ref="db4823b380ba2c3142c857c8c0695200dd1736f6")
},
warning = function(w) stop(w))
