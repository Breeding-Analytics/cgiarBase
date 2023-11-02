nrm2 <- function(
    pedData= NULL,
    verbose=FALSE
){
  ## THIS FUNCTION CALCULATES A NUMERATOR RELATIONSHIP MATRIX
  ## IS USED IN THE BANAL APP UNDER THE STRUCTURE MODULES
  ############################
  # loading the dataset
  if (is.null(pedData)) stop("No input marker data file specified.")
  pedData <- unique(pedData)
  colnames(pedData) <- c("indiv","dam","sire")
  for(j in c("sire","dam","indiv")){
    pedData[,j] <- cgiarBase::replaceValues(pedData[,j], Search = "", Replace = NA)
    pedData[,j] <- cgiarBase::replaceValues(pedData[,j], Search = " ", Replace = NA)
    pedData[,j] <- cgiarBase::replaceValues(pedData[,j], Search = "  ", Replace = NA)
  }
  pedData <- unique(pedData)
  if(nrow(pedData) == length(which(is.na(pedData$dam)))){ # no pedigree information available
    A <- diag(nrow(pedData)); rownames(A) <- colnames(A) <- as.character(pedData$indiv)
  }else{
    # "IR78761-B-SATB2-17-1"
    # the ones that should have NA in dam and sire
    badSire <-  which(is.na(pedData[,2]) & (pedData[,1] == pedData[,3]) )
    if(length(badSire) > 0){ pedData[badSire,3]=NA}
    badDam <- which(is.na(pedData[,3]) & (pedData[,1] == pedData[,2]) )
    if(length(badDam) > 0){ pedData[badDam,2]=NA}
    # the ones who cannot be true
    unrealisticDam <- which(!is.na(pedData[,3]) & (pedData[,1] == pedData[,2]) )
    if(length(unrealisticDam) > 0){ pedData[unrealisticDam,1]=paste0(pedData[unrealisticDam,1],"-2")}
    unrealisticSire <- which(!is.na(pedData[,2]) & (pedData[,1] == pedData[,3]) )
    if(length(unrealisticSire) > 0){ pedData[unrealisticSire,1]=paste0(pedData[unrealisticSire,1],"-2")}
    # the ones that are identical id to mother and father
    unrealisticInd <- which((pedData[,1] == pedData[,2]) & (pedData[,1] == pedData[,3]) )
    if(length(unrealisticInd) > 0){ pedData[unrealisticInd,2]=NA;pedData[unrealisticInd,3]=NA}
    # unique
    pedData <- unique(pedData)
    first <- which(is.na(pedData[,2]) & is.na(pedData[,3])) # no dams and sire
    # second <- which(is.na(pedData[,3]) & !is.na(pedData[,2])) # no sire
    # third <- which(is.na(pedData[,2]) & !is.na(pedData[,3])) # no dam 9313
    foundersName <- na.omit(unique(pedData[first,"indiv"]))
    damsToAddName <- na.omit(setdiff(unique(pedData[,"dam"]),c(unique(pedData[,"indiv"]),foundersName ) )) # dams not present in inviduals or founders
    siresToAddName <- na.omit(setdiff(unique(pedData[,"sire"]),c(damsToAddName, unique(pedData[,"indiv"]), foundersName) ))
    progenyToAddName <- na.omit(setdiff(unique(pedData[,"indiv"]),c(damsToAddName, siresToAddName, foundersName) ))
    ids <- na.omit(unique(c(foundersName,damsToAddName, siresToAddName,progenyToAddName)))
    #
    idsDf <- data.frame(idsn = 1:length(ids)); rownames(idsDf) <- as.character(ids)
    #
    idsDfInverse <- data.frame(ids = ids); rownames(idsDfInverse) <- 1:length(ids)
    # matrix of numerical ids for individuals with information for father, mother or both
    mydataWithPedigreeOnly <- pedData[which(pedData[,"indiv"] %in% progenyToAddName),]
    mydataWithPedigreeOnly <- mydataWithPedigreeOnly[which(!duplicated(mydataWithPedigreeOnly[,"indiv"])),]
    orPedN <- apply(mydataWithPedigreeOnly,2,function(x){ # not sure why is working for the  "" genotypes but is working :)
      idsDf[as.character(x),]
    })
    orPedN <- orPedN[ order(orPedN[,3],orPedN[,2]), ]
    # matrix of numerical ids for individuals without father and mother
    orPedNfounders <- matrix(NA, ncol=3, nrow=length(c(foundersName,damsToAddName, siresToAddName)))
    orPedNfounders[,1] <- 1:length(c(foundersName,damsToAddName, siresToAddName))
    colnames(orPedNfounders) <- colnames(orPedN)
    # bind both
    orPedN2 <- rbind(orPedNfounders,orPedN)
    pede <- data.frame(sire = as.character(orPedN2[,3]),
                       dam  = as.character(orPedN2[,2]),
                       label= as.character(orPedN2[,1])
    )
    pede2<- pedigreemm::editPed(sire=pede$sire, dam= pede$dam, label=pede$label, verbose = FALSE)
    ped<- with(pede2, pedigreemm::pedigree(label=label, sire=sire, dam=dam))
    
    A <- as.matrix(pedigreemm::getA(ped))
    rownames(A) <- colnames(A) <- as.character(idsDfInverse[rownames(A),])
  }
  return(A)
}
