tabToMat <- function(Adf){
  A <- matrix(NA, nrow=max(Adf$Var1n,na.rm = TRUE), ncol = max(Adf$Var2n, na.rm = TRUE))
  A[as.matrix(Adf[,c("Var1n","Var2n")])] = Adf[,"Freq"]
  A[lower.tri(A)] <- t(A)[lower.tri(A)] # fill the lower triangular
  rownames(A) <- colnames(A) <- levels(Adf$Var1)
  return(A)
}


