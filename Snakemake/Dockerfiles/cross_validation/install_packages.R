withCallingHandlers({
  install.packages("lhs", repos="https://cloud.r-project.org/")
  install.packages("rBayesianOptimization", repos="https://cloud.r-project.org/")
},
warning = function(w) stop(w))
