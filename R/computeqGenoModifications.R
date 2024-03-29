computeGenoModifications <- function(
    M = NULL, # markers are assumed to come centered
    propNaUpperThreshForMarker=.3,
    propNaUpperThreshForInds=.3,
    propHetUpperThreshForMarker = 1, 
    propFisUpperThreshForMarker = 1,
    maf=0.05, ploidy=2,
    imputationMethod="median"
){
  ## THIS FUNCTION IMPUTES A MARKER MATRIX
  ## IS USED IN THE BCLEAN APP UNDER THE QA/QC MODULES
  qamAnalysisId <- as.numeric(Sys.time())
  modifications <- as.data.frame(matrix(nrow=0, ncol=6))
  colnames(modifications) <- c("module" ,     "analysisId"  ,"reason"  ,     "row" ,"col"  , "value" )
  ############################
  # loading the dataset
  ## identify which individuals should be removed according to the threshold
  propNaIndividual <- apply(M,1,function(x){(length(which(is.na(x)))/length(x))})
  badIndividual <- which(propNaIndividual > propNaUpperThreshForInds)
  if(length(badIndividual) > 0){
    modifications <- rbind(modifications, data.frame(module="qaGeno", analysisId=qamAnalysisId, reason="%missing", row=badIndividual, col=NA, value=NA) )
  }
  ## identify which markers should be removed according to the threshold
  propNaMarker <- apply(M,2,function(x){(length(which(is.na(x)))/length(x))})
  badMarker <- which(propNaMarker > propNaUpperThreshForMarker)
  if(length(badMarker) > 0){
    modifications <- rbind(modifications, data.frame(module="qaGeno", analysisId=qamAnalysisId, reason="%missing", row=NA, col=badMarker, value=NA) )
  }
  ## identify bad MAF
  MAF <- apply(M+1, 2, function(x) {
    AF <- mean(x, na.rm = T)/ploidy
    MAF <- ifelse(AF > 0.5, 1 - AF, AF)
  })
  badMarker2 <- which(MAF < maf)
  if(length(badMarker2) > 0){
    modifications <- rbind(modifications, data.frame(module="qaGeno", analysisId=qamAnalysisId, reason="MAF", row=NA, col=badMarker2, value=NA) )
  }
  ## heterozigosity per marker
  q <- MAF; p <- 1-q
  he <- 2 * p * q
  ho <- colMeans((M+1) == 1, na.rm = TRUE) # Get the obseved heterozygosity.
  badMarker3 <- which(ho > propHetUpperThreshForMarker)
  if(length(badMarker3) > 0){
    modifications <- rbind(modifications, data.frame(module="qaGeno", analysisId=qamAnalysisId, reason="heterozygosity", row=NA, col=badMarker3, value=NA) )
  }
  ## inbreeding per marker Fis
  Fis <- ifelse(he == 0, yes = 0, no = 1 - (ho / he))
  badMarker4 <- which(Fis > propFisUpperThreshForMarker)
  if(length(badMarker4) > 0){
    modifications <- rbind(modifications, data.frame(module="qaGeno", analysisId=qamAnalysisId, reason="inbreeding", row=NA, col=badMarker4, value=NA) )
  }
  ###################
  ###################
  # imputation track
  toImpute <- which(is.na(M), arr.ind = TRUE)
  remove1 <- which(toImpute[,1] %in% badIndividual)
  if(length(remove1) > 0){toImpute <- toImpute[-remove1,]}
  remove2 <- which(toImpute[,2] %in% badMarker)
  if(length(remove2) > 0){toImpute <- toImpute[-remove2,]}
  remove3 <- which(toImpute[,2] %in% badMarker2)
  if(length(remove3) > 0){toImpute <- toImpute[-remove3,]}
  remove4 <- which(toImpute[,2] %in% badMarker3)
  if(length(remove4) > 0){toImpute <- toImpute[-remove4,]}
  remove5 <- which(toImpute[,2] %in% badMarker4)
  if(length(remove5) > 0){toImpute <- toImpute[-remove5,]}
  # impute only markers and indivuals that will mot be removed given the current cleaning parameters
  if(imputationMethod == "median"){
    M <- apply(M,2,sommer::imputev)
    if(nrow(toImpute) > 0){
      modifications <- rbind(modifications, data.frame(module="qaGeno", analysisId=qamAnalysisId, reason="impute", row=toImpute[,1], col=toImpute[,2], value=M[toImpute] ) )
    }
  }else{
    stop("Method not implemented yet.",call. = FALSE)
  }
  # M$modifications$geno <- rbind(M$modifications$geno, modifications)
  #########################################
  ## update the rest of the data tables
  ## write the parameters to the parameter database
  # M$status <- rbind( M$status, data.frame(module="qaGeno", analysisId=qamAnalysisId))

  return(modifications)
}
