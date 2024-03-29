tabToMat <- function(Adf, symmetric=TRUE, V1="Var1n", V2="Var2n", Freq="Freq"){
  A <- matrix(NA, nrow=max(Adf[,V1],na.rm = TRUE), ncol = max(Adf[,V2], na.rm = TRUE))
  A[as.matrix(Adf[,c(V1,V2)])] = Adf[,Freq]
  if(symmetric){
    A[lower.tri(A)] <- t(A)[lower.tri(A)] # fill the lower triangular
    rownames(A) <- colnames(A) <- levels(Adf[,V1])
  }
  return(A)
}


