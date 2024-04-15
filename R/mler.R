mler <- function(fixed, 
                    random=NULL, ginverse=NULL,
                    group, data,
                    yourLayers=NULL, epochs_M=50, 
                    batch_size=50,
                    returnMatrices=FALSE){
  # print("usingML")
  `%>%` <- dplyr::`%>%`
  ## get matrices for random terms
  counter <- 1
  if(!missing(group)){ # we are using rrBLUP
    Zlist <- list()
    for(iTerm in 1:length(group)){ # iTerm=1
      Zlist[[counter]] <- as.matrix(data[,group[[iTerm]]])
      counter <- counter + 1
    }
    names(Zlist) <- names(group)
    # Z <- do.call(cbind, Zlist)
  }else{Zlist <- list()}
  
  
  if(!is.null(random)){
    randomTerms <- gsub(" ", "", strsplit(as.character(random[2]), split = "[+-]")[[1]])
    for(iTerm in randomTerms){ # iTerm <- randomTerms[1]
      Za <- Matrix::sparse.model.matrix(as.formula(paste("~",iTerm,"-1")), data)
      removeNames <- strsplit(iTerm,":")[[1]]
      for(iName in 1:length(removeNames)){
        if(iName == 1){
          colnames(Za) <- gsub(removeNames[iName], "", colnames(Za) )
        }else{
          colnames(Za) <- gsub(paste0(":",removeNames[iName]), "", colnames(Za) )
        }
      }
      if("genoF" %in% removeNames){
        common <- intersect(colnames(ginverse$genoF), colnames(Za))
        different <- setdiff(colnames(ginverse$genoF), colnames(Za))
        if(length(different) > 0){
          Zb <- Matrix::Matrix(0, nrow=nrow(Za), ncol=length(different))
          colnames(Zb) <- different
          Za <- cbind(Za,Zb)
        }
        Za <- Za[,colnames(ginverse$genoF)] %*% ginverse$genoF
      }
      Zlist[[counter]] <- Za
      names(Zlist)[counter] <- iTerm
      counter <- counter + 1
    }
  }
  if(length(Zlist) == 0){
    stop("Please fit a random formula or provide the grouping argument")
  }
  Z <- do.call(cbind, Zlist)
  
  
  ## get response
  response <- strsplit(as.character(fixed[2]), split = "[+]")[[1]]
  responsef <- as.formula(paste(response,"~1"))
  mfna <- try(model.frame(responsef, data = data, na.action = na.pass), silent = TRUE)
  if (is(mfna, "try-error") ) { # class(mfna) == "try-error"
    stop("Please provide the 'data' argument for your specified variables.\nYou may be specifying some variables in your model not present in your dataset.", call. = FALSE)
  }
  mfna <- eval(mfna, data, parent.frame())
  yvar <- Matrix::sparse.model.matrix(as.formula(paste("~",response,"-1")),data)
  ## scale response
  yvarScaled = scale(yvar) # scale TP
  means=apply(yvar,2,mean) # mean for each trait
  sds=apply(yvar,2,sd) # sd for each trait
  ## put them in a list
  yvarList <- vector(mode = "list", length = ncol(yvar))
  for(k in 1:ncol(yvar)){yvarList[[k]] <- yvarScaled[,k]}
  ## get fixed effect model matrices
  data$`1` <-1
  newfixed=fixed
  # fixedTerms <- gsub(" ", "", strsplit(as.character(fixed[3]), split = "[+-]")[[1]])
  mf <- try(model.frame(newfixed, data = data, na.action = na.pass), silent = TRUE)
  mf <- eval(mf, parent.frame())
  X <-  Matrix::sparse.model.matrix(newfixed, mf)
  ## bind matrices
  W <- cbind(X,Z)
  
  if(returnMatrices){
    return(list(yvar=yvar, yvarScaled=yvarScaled, means=means, sds=sds, X=X, Z=Z, Zlist=Zlist, W=W ))
  }else{
    # define dimensions of covariates (independent variables)
    input<- keras::layer_input(shape=dim(W)[2],name="covars")
    # add hidden layers
    base_model <- input 
    for(iRow in 1:nrow(yourLayers)){ # iRow=1
      base_model <- base_model %>% keras::layer_dense(units = yourLayers[iRow,"units"], 
                                                      activation= yourLayers[iRow,"activation"] )%>%
        keras::layer_dropout(rate = yourLayers[iRow,"rate"] )
    }
    ## add output(s) units
    yhatList <- list()
    for(j in 1:ncol(yvar)){
      yhatList[[j]] <- base_model %>% keras::layer_dense(units = 1, name=paste0("yhat",j))
    }
    # build multi-output model
    model <- keras::keras_model(input, outputs=yhatList ) %>%
      keras::compile(optimizer ="rmsprop",
                     loss="mse",
                     metrics="mae",
                     loss_weights=rep(1/nrow(yourLayers), nrow(yourLayers))
      )
    #fitmodel  # ?fit.keras.engine.training.Model
    model_fit <- model %>%
      keras::fit(x=W,
                 y=yvarList,
                 epochs=epochs_M,
                 batch_size = batch_size,
                 verbose=0)
    
    # Yhat <-  predict(model, W) %>%
    #   data.frame() %>%
    #   setNames(colnames(yvar))
    # plot(Yhat$predictedValue)
    
    return(list(model_fit=model_fit,model=model))
  }
}