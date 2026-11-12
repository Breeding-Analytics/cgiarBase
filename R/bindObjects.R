bindObjects <- function(
    object1= NULL,
    object2= NULL,
    verbose=FALSE
){
  ## THIS FUNCTION ROW BINDS TWO MARKER MATRICES AND ASSOCIATED TABLES
  if(is.null(object1)){stop("Please provide the name of the file to be used for analysis", call. = FALSE)}
  if(is.null(object2)){stop("Please provide the name of the file to be used for analysis", call. = FALSE)}

  # base0 <- list(pheno=data.frame(), pedigree=data.frame(), geno=data.frame(), weather=data.frame(), qtl=data.frame()  )
  # mainElements <- list(base0,base0,base0); names(mainElements) <- c("data", "metadata", "modifications")
  # mainElements <- lapply(vector(mode="list",3), function(x,nn){resx <- vector(mode="list",nn); names(resx) <- c("pheno","pedigree","weather","qtl","geno");return(resx)}, nn=5)
  # names(mainElements) <- c("data", "metadata", "modifications")
  # mainElements$status <- mainElements$modeling <- mainElements$metrics <- mainElements$predictions <- data.frame()
  mainElements <- cgiarBase::create_getData_object()
  '%!in%' <- function(x,y)!('%in%'(x,y))
  ###################################
  # loading the dataset
  for(iMain in 1:length(mainElements)){ #  for each element of data, metadata, modifications, modeling, status, etc ... bind   iMain=1
    if( names(mainElements)[iMain] %in% c("data","metadata","modifications") ){ # user may not have same format
      subItems <- names(mainElements[[names(mainElements)[iMain]]])
      for(iSub in 1:length(mainElements[[names(mainElements)[iMain]]]) ){ # for each nested element   iSub=1
        if(names(mainElements)[iMain] == "data"){ # this are complex, data can be in completely different formats, column names, etc
          # we have to modify both data and metadata at the same time
          if(subItems[iSub] %in% c("pheno","pedigree","qtl") ){
            # bind data tables
            object1$metadata[[subItems[iSub]]] <- object1$metadata[[subItems[iSub]]][which(object1$metadata[[subItems[iSub]]]$parameter %!in% c("", "source", "sourceId")  & object1$metadata[[subItems[iSub]]]$value %!in% c("", "source", "sourceId")),]
            residTrait1 <- paste0(object1$metadata[[subItems[iSub]]][which(object1$metadata[[subItems[iSub]]]$parameter == "trait"),"value"],"-residual")
            residTrait1 <- residTrait1[which(residTrait1 %in% colnames(object1$data[[subItems[iSub]]]))]
            provPheno1 <- object1$data[[subItems[iSub]]][ ,c(object1$metadata[[subItems[iSub]]]$value,residTrait1), drop=FALSE]
            if(!is.null(provPheno1)){
              if(subItems[iSub] %in% c("pheno","pedigree") ){
                colnames(provPheno1) <- cgiarBase::replaceValues(Source = colnames(provPheno1),
                                                                 Search = object1$metadata[[subItems[iSub]]]$value[object1$metadata[[subItems[iSub]]]$parameter %!in% c("trait","source","sourceId")],
                                                                 Replace = object1$metadata[[subItems[iSub]]]$parameter[object1$metadata[[subItems[iSub]]]$parameter %!in% c("trait","source","sourceId")])
              }
            }
            object2$metadata[[subItems[iSub]]] <- object2$metadata[[subItems[iSub]]][which(object2$metadata[[subItems[iSub]]]$parameter %!in% c("", "source", "sourceId")  & object2$metadata[[subItems[iSub]]]$value %!in% c("", "source", "sourceId") ),]
            residTrait2 <- paste0(object2$metadata[[subItems[iSub]]][which(object2$metadata[[subItems[iSub]]]$parameter == "trait"),"value"],"-residual")
            residTrait2 <- residTrait2[which(residTrait2 %in% colnames(object2$data[[subItems[iSub]]]))]
            provPheno2 <- object2$data[[subItems[iSub]]][ ,c(object2$metadata[[subItems[iSub]]]$value,residTrait2), drop=FALSE]
            if(!is.null(provPheno2)){
              if(subItems[iSub] %in% c("pheno","pedigree") ){
              colnames(provPheno2) <- cgiarBase::replaceValues(Source = colnames(provPheno2),
                                                               Search = object2$metadata[[subItems[iSub]]]$value[object2$metadata[[subItems[iSub]]]$parameter %!in% c("trait","source","sourceId")],
                                                               Replace = object2$metadata[[subItems[iSub]]]$parameter[object2$metadata[[subItems[iSub]]]$parameter %!in% c("trait","source","sourceId")])
              }
            }
            allNames <- unique(c(colnames(provPheno1), colnames(provPheno2)))
            commonNames <- intersect(colnames(provPheno1), colnames(provPheno2)) #use intersect in a list of vectors
            differentNames <- c( setdiff(colnames(provPheno1), colnames(provPheno2)) , setdiff(colnames(provPheno2), colnames(provPheno1)) )
            if(!is.null(provPheno1)){
              dfg1 <- data.frame(matrix(ncol = length(allNames), nrow = nrow(provPheno1))); colnames(dfg1) <- allNames
              dfg1[,colnames(provPheno1)] <- provPheno1
              dfg1 <- dfg1[,c(commonNames,differentNames), drop=FALSE]
            }else{dfg1 <- data.frame()}
            if(!is.null(provPheno2)){
              dfg2 <- data.frame(matrix(ncol = length(allNames), nrow = nrow(provPheno2))); colnames(dfg2) <- allNames
              dfg2[,colnames(provPheno2)] <- provPheno2
              dfg2 <- dfg2[,c(commonNames,differentNames), drop=FALSE]
            }else{dfg2 <- data.frame()}
            newData <- unique( rbind(dfg1,dfg2) )
            rownames(newData) <- NULL
            mainElements[[names(mainElements)[iMain]]][[subItems[iSub]]] <- newData
            # update metadata
            newMetadata <- rbind(object1$metadata[[subItems[iSub]]],object2$metadata[[subItems[iSub]]])
            traits <- newMetadata[which(newMetadata$parameter %in% c("trait")),]; traits <- traits[which(!duplicated(traits$value)),]
            if(subItems[iSub] %in% c("pheno","pedigree") ){
            newMetadata <- newMetadata[which(!duplicated(newMetadata$parameter)),]
            newMetadata <- newMetadata[which(!newMetadata$parameter == "trait"),]
            newMetadata$value <- newMetadata$parameter
            }
            newMetadata <- unique(rbind(newMetadata, traits))
            mainElements$metadata[[subItems[iSub]]] <- newMetadata
          }else if(subItems[iSub] == "geno"){ # if we need to fill the genotype slot
            if( inherits(object1$data$geno, "genlight") || inherits(object2$data$geno, "genlight") ){ # genlight format (current Bioflow)
              genoRes <- bindGenlightSlots(object1, object2)
              mainElements$data$geno <- genoRes$geno
              mainElements$metadata$geno <- genoRes$metadata
            }else if( !is.null(object1$data$geno) &  !is.null(object2$data$geno) ){ # ifboth objects have genotype data
              uniqueMarkers <- unique( c( colnames(object1$data$geno), colnames(object2$data$geno) ))
              # only keep new individuals, ignore the others
              newInds <- setdiff( rownames(object2$data$geno), rownames(object1$data$geno) )
              Mx <- matrix(NA, nrow=nrow(object1$data$geno)+nrow(object2$data$geno), ncol=length(uniqueMarkers) ) # new marker matrix
              colnames(Mx) <- uniqueMarkers
              metaX <- data.frame(marker=uniqueMarkers, chr=NA, pos=NA, refAllele=NA, altAllele=NA ) # new metadata for geno
              if(length(newInds) > 0){ploidyInObject2 <- round(diff(range( object2$data$geno, na.rm=TRUE)))}else{ploidyInObject2 <- round(diff(range( object1$data$geno, na.rm=TRUE)))}
              for(iMarker in uniqueMarkers){ # for each unique marker  iMarker = uniqueMarkers[2]
                presentIn1 <- intersect(colnames(object1$data$geno), iMarker)
                presentIn2 <- intersect(colnames(object2$data$geno), iMarker)
                if( (length(presentIn1)+length(presentIn2)) == 2 ){ # if marker is present in both objects
                  # print(iMarker)
                  meta1 <- object1$metadata$geno[object1$metadata$geno$marker == iMarker,]
                  if(is.na(meta1$refAllele)){meta1$refAllele <- meta1$altAllele}
                  if(is.na(meta1$altAllele)){meta1$altAllele <- meta1$refAllele}
                  meta2 <- object2$metadata$geno[object2$metadata$geno$marker == iMarker,]
                  if(is.na(meta2$refAllele)){meta2$refAllele <- meta2$altAllele}
                  if(is.na(meta2$altAllele)){meta2$altAllele <- meta2$refAllele}
                  differentialAlleles <- setdiff(c(meta1$refAllele, meta1$altAllele), c(meta2$refAllele, meta2$altAllele) )
                  if(length(differentialAlleles) > 0){ # multi-allelic marker, we set to NA all the calls
                    Mx[,iMarker] <- NA # we ignore this ones
                    metaX[metaX$marker == iMarker,] <- meta1
                    metaX[metaX$marker == iMarker,"refAllele"] <- paste(meta1$refAllele, meta2$refAllele, collapse = ",")
                    metaX[metaX$marker == iMarker,"altAllele"] <- paste(meta1$altAllele, meta2$altAllele, collapse = ",")
                  }else{ # confirmed bi-allelic marker, we proceed
                    if(meta1$refAllele == meta2$refAllele){ # we have same reference allele, we do a straigth bind
                      Mx[,iMarker] <- c(object1$data$geno[,iMarker], object2$data$geno[,iMarker])
                    }else{ # we have different reference allele, transform marker matrix
                      Mx[,iMarker] <- c(object1$data$geno[,iMarker], (object2$data$geno[,iMarker] - ploidyInObject2) )
                    }
                    metaX[metaX$marker == iMarker,] <- meta1
                  }
                }else{ # if marker is only present in one object
                  if(length(presentIn1) > 0){ # marker is in object 1
                    Mx[,iMarker] <- c(object1$data$geno[,iMarker], rep(NA, nrow(object2$data$geno) ) )
                    metaX[metaX$marker == iMarker,] <-  object1$metadata$geno[object1$metadata$geno$marker == iMarker,]
                  }else if(length(presentIn2) > 0){ # marker is in object 2
                    Mx[,iMarker] <- c( rep(NA, nrow(object1$data$geno) ), object2$data$geno[,iMarker] )
                    metaX[metaX$marker == iMarker,] <-  object2$metadata$geno[object2$metadata$geno$marker == iMarker,]
                  }
                } # end of:  if( (length(presentIn1)+length(presentIn2)) == 2)
              } # end of:  for(iMarker in uniqueMarkers)
              ## there may be duplicated individuals, keep the best record from each individual
              concatenatedRownames <-  c( rownames(object1$data$geno), rownames(object2$data$geno) )
              checkIndsReps <- table( concatenatedRownames )
              indsToCorrect <- unique(names(which(checkIndsReps > 1)))
              if(length(indsToCorrect) > 0){
                toRemove <- vector(mode = "list", length = length(indsToCorrect) ); counter <- 1
                for(iInd in indsToCorrect){ # iInd = indsToCorrect[1]
                  provInd <- which(concatenatedRownames %in% iInd)
                  propNa <- apply(Mx[provInd,],1,function(x){length(which(is.na(x)))/length(x)})
                  keepRow <- provInd[which(propNa == min(propNa))[1]] # lowest mssing data
                  toRemove[[counter]] <- setdiff(provInd, keepRow); counter <- counter + 1 # higher missing data
                }
                Mx <- Mx[-unlist(toRemove),]
                concatenatedRownames <- concatenatedRownames[-unlist(toRemove)]
              }
              rownames(Mx) <- concatenatedRownames
              mainElements[[names(mainElements)[iMain]]][[subItems[iSub]]] <- Mx # store new markers
              mainElements$metadata[[subItems[iSub]]] <- metaX # store new metadata
            }else{ # only one object has genotype data
              if(!is.null(object1$data$geno)){ # markers are in object 1
                mainElements[[names(mainElements)[iMain]]][[subItems[iSub]]] <- object1$data$geno
                mainElements$metadata[[subItems[iSub]]] <- object1$metadata$geno
              }else if(!is.null(object2$data$geno)){ # markers are in object 2
                mainElements[[names(mainElements)[iMain]]][[subItems[iSub]]] <- object2$data$geno
                mainElements$metadata[[subItems[iSub]]] <- object2$metadata$geno
              }
            } # end of: if( !is.null(object1$data$geno) &  !is.null(object2$data$geno) )
          }else if(subItems[iSub] == "weather"){
            if(!is.null(object1$data[[subItems[iSub]]]) | !is.null(object2$data[[subItems[iSub]]])){
              mainElements$data[[subItems[iSub]]] <- unique( rbind(object1$data[[subItems[iSub]]], object2$data[[subItems[iSub]]] ) )
              mainElements$metadata[[subItems[iSub]]] <- unique( rbind(object1$metadata[[subItems[iSub]]], object2$metadata[[subItems[iSub]]] ) )
            }
          } # end of: if(subItems[iSub] %in% c("pheno","pedigree","weather","qtl") )
        }else if(names(mainElements)[iMain] == "modifications"){ #  modifications  can be easily binded, is just long format tables
          if(!is.null(object1[[names(mainElements)[iMain]]][[subItems[iSub]]]) | !is.null(object2[[names(mainElements)[iMain]]][[subItems[iSub]]])){
            mainElements[[names(mainElements)[iMain]]][[subItems[iSub]]] <-  unique( rbind( object1[[names(mainElements)[iMain]]][[subItems[iSub]]], object2[[names(mainElements)[iMain]]][[subItems[iSub]]] ) )
          }
        }
      } # end of: for(iSub in 1:length(mainElements[[names(mainElements)[iMain]]]) )
    }else{ # user wants to do a pure row bind (predictions, metrics, modeling, status tables)
      if(!is.null(object1[[names(mainElements)[iMain]]])){
        if(names(mainElements)[iMain] == "status"){
          if ("analysisIdName" %!in% colnames(object1[[names(mainElements)[iMain]]])){
            object1[[names(mainElements)[iMain]]]$analysisIdName <- ""
          }
        }else if(names(mainElements)[iMain] == "predictions"){
          if ("effectType" %!in% colnames(object1[[names(mainElements)[iMain]]])){
            object1[[names(mainElements)[iMain]]]$effectType <- NA
          }
        }
      }

      if(!is.null(object2[[names(mainElements)[iMain]]])){
        if(names(mainElements)[iMain] == "status"){
          if ("analysisIdName" %!in% colnames(object2[[names(mainElements)[iMain]]])){
            object2[[names(mainElements)[iMain]]]$analysisIdName <- ""
          }
        }else if(names(mainElements)[iMain] == "predictions"){
          if ("effectType" %!in% colnames(object2[[names(mainElements)[iMain]]])){
            object2[[names(mainElements)[iMain]]]$effectType <- NA
          }
        }
      }

      if(!is.null(object1[[names(mainElements)[iMain]]]) | !is.null(object2[[names(mainElements)[iMain]]])){
        mainElements[[names(mainElements)[iMain]]] <-  unique( rbind( object1[[names(mainElements)[iMain]]], object2[[names(mainElements)[iMain]]] ) )
        # print(mainElements[[names(mainElements)[iMain]]])
        if(nrow(mainElements[[names(mainElements)[iMain]]]) == 0){mainElements[[names(mainElements)[iMain]]] <- NULL}
      }
    }
  }

  # carry over slots not defined in create_getData_object() (e.g. modifications$selection,
  # data$geno_imp, data$TPP, metadata$TPP, metadata$tpp_analysis_config) so they are not
  # silently dropped. Named lists (keyed by analysisId / TPP id) are unioned by name and
  # object1 wins on duplicated names.
  for(iMain in c("data", "metadata", "modifications")){
    extraSlots <- setdiff(unique(c(names(object1[[iMain]]), names(object2[[iMain]]))),
                          names(mainElements[[iMain]]))
    for(iSlot in extraSlots){
      x1 <- object1[[iMain]][[iSlot]]; x2 <- object2[[iMain]][[iSlot]]
      if(is.data.frame(x1) && is.data.frame(x2) && nrow(x1) > 0 && nrow(x2) > 0){
        allCols <- unique(c(colnames(x1), colnames(x2)))
        for(cc in setdiff(allCols, colnames(x1))) x1[[cc]] <- NA
        for(cc in setdiff(allCols, colnames(x2))) x2[[cc]] <- NA
        newX <- unique(rbind(x1[, allCols, drop = FALSE], x2[, allCols, drop = FALSE]))
        rownames(newX) <- NULL
      }else if(is.list(x1) && is.list(x2) && !is.data.frame(x1) && !is.data.frame(x2)){
        newX <- c(x1, x2[setdiff(names(x2), names(x1))]) # named lists keyed by analysisId
      }else if(is.data.frame(x1) && nrow(x1) > 0){
        newX <- x1
      }else if(!is.null(x2) && !(is.data.frame(x2) && nrow(x2) == 0)){
        newX <- x2
      }else{
        newX <- x1
      }
      if(!is.null(newX)) mainElements[[iMain]][[iSlot]] <- newX
    }
  }

  return(mainElements)
}

## Internal: bind the genotype slots of two Bioflow objects when data$geno is an
## adegenet genlight and metadata$geno is a parameter/value table.
## - Same loci in both objects: rows are stacked (rbind.genlight).
## - Different loci: union of loci; calls missing in one object are NA. For shared
##   loci with reversed alleles (A/G vs G/A) the second object's dosage is flipped;
##   for incompatible alleles the locus is set to NA in the second object.
## - Duplicated individuals keep the record from object1.
bindGenlightSlots <- function(object1, object2){
  if(!requireNamespace("adegenet", quietly = TRUE)){
    stop("Package 'adegenet' is required to bind genlight genotype data.", call. = FALSE)
  }
  g1 <- object1$data$geno; g2 <- object2$data$geno
  m1 <- object1$metadata$geno; m2 <- object2$metadata$geno
  isGl1 <- inherits(g1, "genlight") && adegenet::nInd(g1) > 0
  isGl2 <- inherits(g2, "genlight") && adegenet::nInd(g2) > 0
  # only one side has genotype data
  if(!isGl2){ return(list(geno = g1, metadata = m1)) }
  if(!isGl1){ return(list(geno = g2, metadata = m2)) }

  # metadata: parameter/value table, object1 wins on duplicated parameters
  metaX <- rbind(m1, m2)
  if(!is.null(metaX) && all(c("parameter", "value") %in% colnames(metaX))){
    p1 <- m1$value[m1$parameter == "ploidity"]; p2 <- m2$value[m2$parameter == "ploidity"]
    if(length(p1) && length(p2) && as.character(p1[1]) != as.character(p2[1])){
      stop("Cannot bind genotype data with different ploidy levels (", p1[1], " vs ", p2[1], ").", call. = FALSE)
    }
    metaX <- metaX[!duplicated(metaX$parameter), , drop = FALSE]
    rownames(metaX) <- NULL
  }

  ploidy1 <- unique(adegenet::ploidy(g1)); ploidy2 <- unique(adegenet::ploidy(g2))
  if(length(ploidy1) == 1 && length(ploidy2) == 1 && ploidy1 != ploidy2){
    stop("Cannot bind genotype data with different ploidy levels (", ploidy1, " vs ", ploidy2, ").", call. = FALSE)
  }

  loc1 <- adegenet::locNames(g1); loc2 <- adegenet::locNames(g2)
  newInds <- setdiff(adegenet::indNames(g2), adegenet::indNames(g1))

  if(identical(loc1, loc2)){ # same marker set and order: straight row bind
    if(length(newInds) == 0){ return(list(geno = g1, metadata = metaX)) }
    gx <- adegenet::rbind.genlight(g1, g2[adegenet::indNames(g2) %in% newInds, ])
    return(list(geno = gx, metadata = metaX))
  }

  # different marker sets: union of loci through dosage matrices
  M1 <- adegenet::as.matrix(g1)
  M2 <- adegenet::as.matrix(g2)[adegenet::indNames(g2) %in% newInds, , drop = FALSE]
  uniqueMarkers <- unique(c(loc1, loc2))
  Mx <- matrix(NA_real_, nrow = nrow(M1) + nrow(M2), ncol = length(uniqueMarkers),
               dimnames = list(c(rownames(M1), rownames(M2)), uniqueMarkers))
  Mx[seq_len(nrow(M1)), loc1] <- M1
  if(nrow(M2) > 0){
    rows2 <- nrow(M1) + seq_len(nrow(M2))
    Mx[rows2, loc2] <- M2
    al1 <- adegenet::alleles(g1); al2 <- adegenet::alleles(g2)
    if(!is.null(al1) && !is.null(al2)){
      names(al1) <- loc1; names(al2) <- loc2
      shared <- intersect(loc1, loc2)
      a1 <- al1[shared]; a2 <- al2[shared]
      differ <- shared[a1 != a2]
      if(length(differ) > 0){
        flipped <- vapply(strsplit(al2[differ], "/", fixed = TRUE),
                          function(x) paste(rev(x), collapse = "/"), character(1))
        toFlip <- differ[flipped == al1[differ]]
        toNA <- setdiff(differ, toFlip)
        pl <- if(length(ploidy2) == 1) ploidy2 else 2
        if(length(toFlip)) Mx[rows2, toFlip] <- pl - Mx[rows2, toFlip]
        if(length(toNA)) Mx[rows2, toNA] <- NA
      }
    }
  }
  gx <- new("genlight", Mx, ploidy = if(length(ploidy1) == 1) ploidy1 else 2)
  # carry over locus annotation (object1 first, then loci only present in object2)
  getAnn <- function(f, g, loc){ v <- f(g); if(is.null(v)) NULL else stats::setNames(as.character(v), loc) }
  chr1 <- getAnn(adegenet::chromosome, g1, loc1); chr2 <- getAnn(adegenet::chromosome, g2, loc2)
  if(!is.null(chr1) || !is.null(chr2)){
    adegenet::chromosome(gx) <- unname(c(chr1, chr2)[uniqueMarkers])
  }
  pos1 <- getAnn(adegenet::position, g1, loc1); pos2 <- getAnn(adegenet::position, g2, loc2)
  if(!is.null(pos1) || !is.null(pos2)){
    adegenet::position(gx) <- unname(c(pos1, pos2)[uniqueMarkers])
  }
  al1 <- getAnn(adegenet::alleles, g1, loc1); al2 <- getAnn(adegenet::alleles, g2, loc2)
  if(!is.null(al1) || !is.null(al2)){
    adegenet::alleles(gx) <- unname(c(al1, al2)[uniqueMarkers])
  }
  list(geno = gx, metadata = metaX)
}
