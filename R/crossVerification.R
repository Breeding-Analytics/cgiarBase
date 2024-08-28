

crossVerification <- function(Mf,Mm,Mp, ploidy=2){
  
  consensus <- function(M, FUN){
    # take a set of samples and get the consensus sequence
    Mnew <- apply(M,2,FUN)
    return(Mnew)
  }
  # take the markers of the progeny (Mp) and 
  # calculate the probability of belonging to a given male and female
  if(any(!is.matrix(Mf) | !is.matrix(Mm) | !is.matrix(Mp) )){
    stop("One or more of the provided arguments is not a matrix", call. = FALSE)
  }
  Mf1 <- (Mf + Mm)/2 #apply(rbind(Mf,Mm),2,mean) # expected genotype
  informativeMarkers <- which( apply(rbind(Mf,Mm),2,var) > 0 ) # markers with variance in parents at the pop level

  # matrix to save results
  resultMatch <- Matrix::Matrix(0, nrow = nrow(Mf1), ncol=ncol(Mf1))
  # matching matrix
  match1 <- abs(Mf1 - Mp)
  
  # when both expected and realized progeny have same call is a full match
  resultMatch[which(match1 == 0, arr.ind = TRUE)]=1 
  # when expected and realized progeny was of the form P1 x P2 = 1 x 0
  # expected progeny could be either 0 or 1 so Mf1=0.5 so if Mp=0 or Mp=1 
  # Mf1 - Mp will be 0.5 or -0.5
  # we get a 0.5, it is a match but we have some uncertainty
  resultMatch[which(match1 == (0.5*(ploidy/2)), arr.ind = TRUE)]=0.5
  # when expected and realized progeny was of the form P1 x P2 = 2 x 1
  # expected progeny could be either 2 or 1 so Mf1=1.5 so if Mp=2 or Mp=1 
  # Mf1 - Mp will be 0.5 or -0.5
  # we get a 0.5, it is a match but we have some uncertainty
  
  # when expected and realized progeny was of the form P1 x P2 = 1 x 1
  # expected progeny could be either 2 or 1 or 0 so Mf1=1 so if Mp=2 or Mp=1 or Mp=0
  # Mf1 - Mp will be 1 or 0 or 0, no way to track it, intead we go back to those
  # specific cases and assign directly a 0.25 probability
  # we get a 0.25, it is a match but we have some uncertainty
  resultMatch[which(Mf1 == (1*(ploidy/2)), arr.ind = TRUE)]=0.25
  
  # compute metrics
  probMatch <- apply(resultMatch,1,function(x){sum(x)/length(x)})  # proportion of markers matching
  heteroMp <- apply(Mp,1,function(x){length(which(x == 1))/length(x)}) # heterozigosity found in progeny
  heteroMf1 <- apply(Mf1,1,function(x){length(which(x == 1))/length(x)}) #heterozigosity expected in progeny
  heteroDeviation <- abs(heteroMf1 - heteroMp) # deviation from expectedheterosigosity
  propComplete <- apply(Mp,1,function(x){length(which(!is.na(x)))/length(x)}) # proportion of complete cases
  
  res <- data.frame(designation=rownames(Mp),probMatch,heteroMp,heteroMf1,heteroDeviation,propComplete)
  rownames(res) <- NULL
  return(list(metrics=res,matchMat=resultMatch))
  
}

