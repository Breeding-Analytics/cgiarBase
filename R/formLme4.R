goodLevels <- function(object, analysisId, includeCovars=TRUE){
  '%!in%' <- function(x,y)!('%in%'(x,y))
  metaPheno <- object$metadata$pheno
  predictions <- object$predictions[which(object$predictions$analysisId == analysisId ),]
  # add missing columns
  keep <- which( metaPheno$parameter %!in% c("trait","designation","environment","rep","row","col","iBlock") )
  if(length(keep) > 0){
    toExtractFromData <- metaPheno[keep, "value"]
    tpe <- unique(object$data$pheno[,c("environment",toExtractFromData)])
    colnames(tpe) <- cgiarBase::replaceValues(colnames(tpe), Search = metaPheno$value, Replace = metaPheno$parameter )
    predictions <- merge(predictions, tpe, by="environment", all.x = TRUE)
  }
  # add missing weather
  weather <- cgiarPipeline::summaryWeather(object)# object$metadata$weather
  if(nrow(weather) > 0){
    weather$traitParameter <- paste(weather$trait, weather$parameter, sep="_")
    wide <- reshape(weather[,-which(colnames(weather)%in%c("trait","parameter"))], direction = "wide", idvar = "environment",
                    timevar = c("traitParameter"), v.names = "value", sep= "_")
    colnames(wide) <- gsub("value_","",colnames(wide))
    predictions <- merge(predictions, wide, by="environment", all.x = TRUE)
  }
  
  # extract
  availableTraits <- intersect( metaPheno[which(metaPheno$parameter %in% c("trait")),"value"], unique(predictions$trait) )
  factorPerTrait <- vector(mode="list", length = length(availableTraits)); names(factorPerTrait) <- availableTraits
  availableFactors <- metaPheno[which(metaPheno$parameter %in% c("environment","year","season","location","trial","study","management")),"parameter"]
  
  for(iTrait in availableTraits){ # iTrait = availableTraits[1]
    provPred <- predictions[which(predictions$trait == iTrait),]
    provPred <- provPred[which(!is.na(provPred$predictedValue)),]
    for(iFactor in availableFactors){ # iFactor = availableFactors[1]
      if(length(na.omit(unique(provPred[,iFactor]))) > 1){factorPerTrait[[iTrait]] <- c(factorPerTrait[[iTrait]] , iFactor)}
    }
    if(includeCovars){
      for(iWeather in unique(weather$traitParameter)){ # iWeather = unique(weather$traitParameter)[1]
        if( var(provPred[,iWeather], na.rm=TRUE ) > 0 ){ factorPerTrait[[iTrait]] <- c(factorPerTrait[[iTrait]] , iWeather) }
      }
    }
  }
  return(factorPerTrait)
}

formLme4 <- function(input0,object, analysisId, trait){
  # if(is.null(metaPheno)){stop("Metadata for phenotypes needs to be provided.", call. = FALSE)}
  metaPheno <- object$metadata$pheno
  weather <- cgiarPipeline::summaryWeather(object)# object$metadata$weather # object$metadata$weather
  predictions <- object$predictions[which(object$predictions$analysisId == analysisId ),]
  '%!in%' <- function(x,y)!('%in%'(x,y))
  if(nrow(weather) > 0){
    weather$traitParameter <- paste(weather$trait, weather$parameter, sep="_")
    wide <- reshape(weather[,-which(colnames(weather)%in%c("trait","parameter"))], direction = "wide", idvar = "environment",
                    timevar = c("traitParameter"), v.names = "value", sep= "_")
    colnames(wide) <- gsub("value_","",colnames(wide))
    predictions <- merge(predictions, wide, by="environment", all.x = TRUE)
  }
  keep <- which(metaPheno$parameter %!in% c("trait","designation","environment","rep","row","col","iBlock","entryType") )
  if(length(keep) > 0){
    toExtractFromData <- metaPheno[keep, "value"]
    tpe <- unique(object$data$pheno[,c("environment",toExtractFromData)])
    colnames(tpe) <- cgiarBase::replaceValues(colnames(tpe), Search = metaPheno$value, Replace = metaPheno$parameter )
    predictions <- merge(predictions, tpe, by="environment", all.x = TRUE)
  }
  ## make the formula
  availableTraits <- intersect( metaPheno[which(metaPheno$parameter %in% c("trait")),"value"], unique(predictions$trait) )
  form <- preds <- list()
  for(iTrait in availableTraits){ # iTrait = availableTraits[1]
    predictionsTrait <- predictions[which(predictions$trait == iTrait),]
    columnClass <- unlist(lapply(predictionsTrait, class))
    formulas <- list()
    for(i in 1:length(input0)){ # for each effect to fit   i=1
      input <- input0[[i]]
      
      if( is.null(input$left) ){ # fixed effect
        
        if(input$center == ""){ # fixed effect
          if(is.null(input$right)){
            warning(paste("The term",i, "will be ignored since it was misspecified. Empty intercept and empty slope"))
          }else{
            formulas[[i]] <- paste(input$right, collapse = ":")
          }
        }else{ # | or  ||
          warning(paste("The term",i, "will be ignored since it was misspecified. Empty intercept but structure specified"))
        }
        
      }else if( length(input$left) == 1 ){ # only one term in the left side of the equation
        
        if(input$left == "0"){ # invalid intercept
          warning(paste("The term",i, "will be ignored since it was misspecified. Empty in both sides of the equation"))
        }else if(input$left == "1"){ # simple intercept
          if(input$center == ""){
            warning(paste("The term",i, "will be ignored since it was misspecified. Intercept specified but no structure specified"))
          }else if(input$center == "||"){ # if user assigns 1 on intercept the structure can only be |
            input$center <- "|"
            warning(paste("The term",i, "will shift from || to | since it was misspecified."))
          }else{ # |
            
            newLeftCol <- paste(input[["left"]], collapse = "_")
            if(newLeftCol == "1"){ # simple structure
              formulas[[i]] <-  paste( "( ", newLeftCol, input[["center"]], paste(input[["right"]], collapse=":") ,")" )
            }else{ # complex structure
              if(input$nPC == 0){ # no FA
                predictionsTrait[, newLeftCol] <- apply( predictionsTrait[ , input[["left"]] ] , 1 , function(x){paste(na.omit(x), collapse = "-")}  )
                predictionsTrait[ which(predictionsTrait[,newLeftCol] == ""), newLeftCol] <- NA
                Z <- Matrix::sparse.model.matrix(as.formula(paste("~", newLeftCol,"-1")), na.action = na.pass, data=predictionsTrait)
                colnames(Z) <- gsub(newLeftCol, "", colnames(Z))
                colnames(Z) <- gsub(" ","", colnames(Z))
                for(j in 1:ncol(Z)){predictionsTrait[,colnames(Z)[j]] <- Z[,j]}
                formulas[[i]] <-  paste( "( 0 +", paste(colnames(Z), collapse = " + "), input[["center"]], paste(input[["right"]], collapse=":") ,")" )
              }else{ # FA model
                formulas[[i]] <-  paste( "( 0 +", paste(paste0("PC",1:input$nPC), collapse = " + "), input[["center"]], paste(input[["right"]], collapse=":") ,")" )
              }
            }
            
          }
        }else{ # complex intercept (other than "1")
          
          if(input$center == ""){
            warning(paste("The term",i, "will be ignored since it was misspecified. Intercept specified but no structure specified"))
          }else{ # |
            
            columnClassLeft <- columnClass[setdiff( input$left, c("0","1"))]
            leftFactor <- names(columnClassLeft[which(columnClassLeft %in% c("character","factor") )])
            leftNumeric <- names(columnClassLeft[which(columnClassLeft %in% c("integer","numeric") )])
            
            newLeftCol <- paste(input[["left"]], collapse = "_")
            
            if(newLeftCol == "1"){ # simple structure
              formulas[[i]] <-  paste( "( ", newLeftCol, input[["center"]], paste(input[["right"]], collapse=":") ,")" )
            }else{ # complex structure
              if(input$nPC == 0){ # no FA
                
                if( length(leftFactor) > 0 ){  
                  predictionsTrait[, newLeftCol] <- apply( predictionsTrait[ , input[["left"]], drop=FALSE ] , 1 , function(x){paste(na.omit(x), collapse = "-")}  )
                  predictionsTrait[ which(predictionsTrait[,newLeftCol] == ""), newLeftCol] <- NA
                  if(length(unique(na.omit(predictionsTrait[, newLeftCol]))) > 1){
                    Z <- Matrix::sparse.model.matrix(as.formula(paste("~", newLeftCol,"-1")), na.action = na.pass, data=predictionsTrait)
                    colnames(Z) <- gsub(newLeftCol, "L.", colnames(Z))
                    colnames(Z) <- gsub(" ","", colnames(Z))
                  }else{
                    Z <- Matrix::Matrix(1,ncol=1,nrow=nrow(predictionsTrait))
                    colnames(Z) <- paste0("L.",unique(na.omit(predictionsTrait[, newLeftCol])) )
                    colnames(Z) <- gsub(" ","", colnames(Z))
                  }
                  if( length(leftNumeric) > 0 ){  
                    Z <- cbind(Z, predictionsTrait[, leftNumeric, drop=FALSE])
                    colnames(Z) <- gsub(" ","", colnames(Z))
                  }
                }else{
                  if( length(leftNumeric) > 0 ){  
                    Z <- predictionsTrait[, leftNumeric, drop=FALSE]
                    colnames(Z) <- gsub(" ","", colnames(Z))
                  }
                }
                for(j in 1:ncol(Z)){predictionsTrait[,colnames(Z)[j]] <- Z[,j]}
                formulas[[i]] <-  paste( "( 0 +", paste(colnames(Z), collapse = " + "), input[["center"]], paste(input[["right"]], collapse=":") ,")" )
                
              }else{ # FA model
                
                predictionsTrait[, newLeftCol] <- apply( predictionsTrait[ , setdiff( input$left, c("0","1")), drop=FALSE ] , 1 , function(x){paste(na.omit(x), collapse = "-")}  )
                predictionsTrait[ which(predictionsTrait[,newLeftCol] == ""), newLeftCol] <- NA
                
                H0 <- reshape(predictionsTrait[,c(input$right, newLeftCol, "predictedValue" )], direction = "wide", idvar = input$right,
                              timevar = newLeftCol, v.names = "predictedValue", sep= "_")
                H0 <- apply(H0[,2:ncol(H0),drop=FALSE],2,sommer::imputev)
                colnames(H0) <- gsub("predictedValue_","",colnames(H0))
                Z <- with(predictionsTrait, lme4breeding::smm(lme4breeding::rrm(predictionsTrait[,newLeftCol], H = H0, nPC = input$nPC)) )
                colnames(Z) <- paste(colnames(Z), newLeftCol, sep="_")
                colnames(Z) <- gsub(" ","", colnames(Z))
                
                for(j in 1:ncol(Z)){predictionsTrait[,colnames(Z)[j]] <- Z[,j]}
                formulas[[i]] <- paste( "( 0 +", paste(colnames(Z), collapse = " + "), input[["center"]], paste(input[["right"]], collapse=":") ,")" )
                
              } # ond of nPC
            } # end of complex structure
          }# end of | structure
          
        }
        
      }else{ # multiple terms in the left side (intercepts)
        
        if( length( setdiff( input$left, c("0","1")) ) == 0 ){ # misspecified
          warning(paste("The term",i, "will be ignored since it was misspecified. Intercept 0 + 1 is not allowed"))
        }else{
          if(input$center == ""){
            warning(paste("The term",i, "will be ignored since it was misspecified. Intercept specified but no structure specified"))
          }else{ # || or  |
            
            if(is.null(input$right)){ # no right side
              warning(paste("The term",i, "will be ignored since it was misspecified. Intercept and structure specified but no right side"))
            }else{
              
              columnClassLeft <- columnClass[setdiff( input$left, c("0","1"))]
              leftFactor <- names(columnClassLeft[which(columnClassLeft %in% c("character","factor") )])
              leftNumeric <- names(columnClassLeft[which(columnClassLeft %in% c("integer","numeric") )])
              
              newLeftCol <- paste( setdiff( input$left, c("0","1")) , collapse = "_")
              if(input$nPC == 0){ # no FA
                if( length(leftFactor) > 0 ){  
                  predictionsTrait[, newLeftCol] <- apply( predictionsTrait[ , setdiff( input$left, c("0","1")), drop=FALSE ] , 1 , function(x){paste(na.omit(x), collapse = "-")}  )
                  predictionsTrait[ which(predictionsTrait[,newLeftCol] == ""), newLeftCol] <- NA
                  if(length(unique(na.omit(predictionsTrait[, newLeftCol]))) > 1){
                    Z <- Matrix::sparse.model.matrix(as.formula(paste("~", newLeftCol,"-1")), na.action = na.pass, data=predictionsTrait)
                    colnames(Z) <- gsub(newLeftCol, "L.", colnames(Z))
                    colnames(Z) <- gsub(" ","", colnames(Z))
                  }else{
                    Z <- Matrix::Matrix(1,ncol=1,nrow=nrow(predictionsTrait))
                    colnames(Z) <- paste0( "L.", unique(na.omit(predictionsTrait[, newLeftCol])) )
                    colnames(Z) <- gsub(" ","", colnames(Z))
                  }
                  if( length(leftNumeric) > 0 ){  
                    Z <- cbind(Z, predictionsTrait[, leftNumeric, drop=FALSE])
                    colnames(Z) <- gsub(" ","", colnames(Z))
                  }
                }else{
                  if( length(leftNumeric) > 0 ){  
                    Z <- predictionsTrait[, leftNumeric, drop=FALSE]
                    colnames(Z) <- gsub(" ","", colnames(Z))
                  }
                }
                for(j in 1:ncol(Z)){predictionsTrait[,colnames(Z)[j]] <- Z[,j]}
                formulas[[i]] <-  paste( "( 0 +", paste(colnames(Z), collapse = " + "), input[["center"]], paste(input[["right"]], collapse=":"),")" )
              }else{ # FA model
                
                predictionsTrait[, newLeftCol] <- apply( predictionsTrait[ , setdiff( input$left, c("0","1")), drop=FALSE ] , 1 , function(x){paste(na.omit(x), collapse = "-")}  )
                predictionsTrait[ which(predictionsTrait[,newLeftCol] == ""), newLeftCol] <- NA
                
                H0 <- reshape(predictionsTrait[,c(input$right, newLeftCol, "predictedValue" )], direction = "wide", idvar = input$right,
                              timevar = newLeftCol, v.names = "predictedValue", sep= "_")
                H0 <- apply(H0[,2:ncol(H0),drop=FALSE],2,sommer::imputev)
                colnames(H0) <- gsub("predictedValue_","",colnames(H0))
                Z <- with(predictionsTrait, lme4breeding::smm(lme4breeding::rrm(predictionsTrait[,newLeftCol], H = H0, nPC = input$nPC)) )
                colnames(Z) <- paste(colnames(Z), newLeftCol, sep="_")
                colnames(Z) <- gsub(" ","", colnames(Z))
                for(j in 1:ncol(Z)){predictionsTrait[,colnames(Z)[j]] <- Z[,j]}
                formulas[[i]] <- paste( "( 0 +", paste(colnames(Z), collapse = " + "), input[["center"]], paste(input[["right"]], collapse=":") ,")" )
                
              } # ond of nPC
              # } # end of complex structure
              
            }# end of right-side exist
            
          } # end of | or ||
        }# end of multiple terms in left side properly specified
        
      }# end of multiple terms in the left side
      
    } # end of for each effect
    form[[iTrait]] <- paste( unique(unlist(formulas)), collapse = " + ")
    preds[[iTrait]] <- predictionsTrait
  }
  newPredictions <- do.call(rbind, preds)
  return(list(form=form, predictions=newPredictions, input=input0))
}