summaryData <- function(object){
  nFeatures = unlist(lapply(object$data, ncol))
  mapped = unlist(lapply(object$metadata, nrow))
  nRecords = unlist(lapply(object$data, nrow))
  nMods = unlist(lapply(object$modifications, nrow))
  res = as.data.frame(matrix(0,nrow=5, ncol=5))
  colnames(res) <- c("nFeatures","mapped","nRecords","nDataPoints","modifications")
  rownames(res) <- c("pheno","geno","weather","pedigree","qtl")
  ## add 
  res[names(nFeatures),"nFeatures"] = nFeatures
  res[names(mapped),"mapped"] = mapped
  res[names(nRecords),"nRecords"] = nRecords
  res$nDataPoints = res$mapped * res$nRecords
  res[names(nMods),"modifications"] = nMods
  return(res)
}