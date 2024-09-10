summaryGeno <- function(modifications ){
  
  # xx=object$modifications$geno
  # xx <- xx[which(xx$analysisId == analysisId),]
  suppressedMarkers <- length(which(is.na(modifications$row) & !is.na(modifications$col) & is.na(modifications$value) ))
  suppressedInds <- length(which(!is.na(modifications$row) & is.na(modifications$col) & is.na(modifications$value) ))
  suppressedCells <- length(which(!is.na(modifications$row) & !is.na(modifications$col) & is.na(modifications$value) )) + (suppressedInds*ncol(object$data$geno)) + (suppressedMarkers*nrow(object$data$geno))
  imputedCells <- length(which(!is.na(modifications$row) & !is.na(modifications$col) & !is.na(modifications$value) ))
  final <- data.frame(suppressedMarkers=suppressedMarkers,
                      suppressedInds=suppressedInds, totalSuppressedCells=suppressedCells,
                      imputedCells=imputedCells)
  return(final)
  
}