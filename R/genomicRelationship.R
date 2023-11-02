grm <- function(
    M= NULL,#wd=NULL,
    verbose=FALSE
){
  ## THIS FUNCTION CALCULATES A GENOMIC RELATIONSHIP MATRIX USING GENETIC MARKERS
  ## IS USED IN THE BANAL APP UNDER THE STRUCTURE MODULES
  ############################
  # loading the dataset
  if (is.null(M)) stop("No input marker data file specified.")
  
  ############################
  # calculate the relationship matrix
  missingData <- apply(M,2,function(x){length(which(is.na(x)))/length(x)}) # check percentage of missing data
  M <- M[,which(missingData < .50)] # keep markers with at least 50 % of information
  M <- apply(M,2,function(x){x-(mean(c(min(x),max(x))))}) # center the matrix at 0
  A <- sommer::A.mat(M)
  

  return(A)
}
