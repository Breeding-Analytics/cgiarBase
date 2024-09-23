summaryData <- function(object){
  nFeatures = unlist(lapply(object$data, ncol))
  nRecords = unlist(lapply(object$data, nrow))
  res = data.frame(nFeatures=nFeatures, nRecords=nRecords )
  return(res)
}