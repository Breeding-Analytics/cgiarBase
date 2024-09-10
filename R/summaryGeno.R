summaryGeno <- function(object,analysisId ){
  
  xx=object$modifications$geno
  xx <- xx[which(xx$analysisId == analysisId),]
  suppressedMarkers <- length(which(is.na(xx$row) & !is.na(xx$col) & is.na(xx$value) ))
  suppressedInds <- length(which(!is.na(xx$row) & is.na(xx$col) & is.na(xx$value) ))
  suppressedCells <- length(which(!is.na(xx$row) & !is.na(xx$col) & is.na(xx$value) )) + (suppressedInds*ncol(object$data$geno)) + (suppressedMarkers*nrow(object$data$geno))
  imputedCells <- length(which(!is.na(xx$row) & !is.na(xx$col) & !is.na(xx$value) ))
  final <- data.frame(analysisId=analysisId, suppressedMarkers=suppressedMarkers,
                      suppressedInds=suppressedInds, totalSuppressedCells=suppressedCells,
                      imputedCells=imputedCells)
  return(final)
  
}