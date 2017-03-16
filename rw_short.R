library(TMB)
obj <- MakeADFun(data, parameters, random="X")
opt <- do.call("optim", obj)
