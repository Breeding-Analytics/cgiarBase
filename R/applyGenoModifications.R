applyGenoModifications <- function(M, modifications){
  
  # apply the modifications that only apply to row (individuals with too much missing data)
  onlyRow <- which(!is.na(modifications$row) & is.na(modifications$col))
  if(length(onlyRow) > 0){
    for(iRow in 1:length(onlyRow)){
      M[ modifications$row[onlyRow[iRow]],] <-  modifications$value[onlyRow[iRow]]
    }
  }
  # apply the modifications that only apply to col (markers with too much missing data)
  onlyCol <- which(is.na(modifications$row) & !is.na(modifications$col))
  if(length(onlyCol) > 0){
    for(iCol in 1:length(onlyCol)){
      M[, modifications$col[onlyCol[iCol]]] <-  modifications$value[onlyCol[iCol]]
    }
  }
  # apply the modifications that only apply to col (imputations)
  both <- which(!is.na(modifications$row) & !is.na(modifications$col))
  if(length(both) > 0){
    M[ as.matrix(modifications[both,c("row","col")]) ] <- modifications[both,"value"]
  }
  ## now remove all rows and columns that do not have 100 % complete data
  propNaIndividual <- apply(M,1,function(x){(length(which(is.na(x)))/length(x))})
  badInds <- which(propNaIndividual == 1) # individuals that should be totally removed
  if(length(badInds)>0){M <- M[-badInds,,drop=FALSE]}
  propNaMarker <- apply(M,2,function(x){(length(which(is.na(x)))/length(x))})
  badMarker <- which(propNaMarker == 1)
  if(length(badMarker)>0){M <- M[,-badMarker]}
  return(M)
}