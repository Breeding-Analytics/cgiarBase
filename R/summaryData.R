summaryData <- function(object){
  nFeatures = unlist(lapply(object$data, ncol))
  mapped = unlist(lapply(object$metadata, nrow))
  nRecords = unlist(lapply(object$data, nrow))
  nMods = unlist(lapply(object$modifications, nrow))
  res = data.frame(nFeatures=nFeatures, mapped=mapped, nRecords=nRecords )
  res$nDataPoints = res$mapped * res$nRecords
  res$modifications <- NA
  res[names(nMods),"modifications"] = nMods
  return(res)
}