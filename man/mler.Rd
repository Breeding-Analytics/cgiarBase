\name{mler}
\alias{mler}

\title{
Machine learning model
}
\description{
This function forms the linear equations to fit a machine learning model based on back propagation using the keras and tensorflow libraries.
}
\usage{
mler(fixed, 
        random=NULL, ginverse=NULL,
        group, data,
        yourLayers=NULL, epochs_M=50, 
        batch_size=50,
        returnMatrices=FALSE)
}

\arguments{
  \item{fixed}{
Fixed formula for the desired model fitting.
}
\item{random}{
Random formula for the desired model fitting.
}
\item{ginverse}{
A list where each element contains a relationship matrix for the levels of a random effect specified in the random formula. It should include the individuals in the dataset plus the individuals to be predicted.
}
\item{group}{
A list where each element contains a numeric vector indicating which columns in the dataset contain the incidence matrix to be fitted.
}
\item{data}{
Input dataset.
}
\item{yourLayers}{
A data frame with as many rows as layers of neurons to be fitted. The data frame should contain 3 columns named 'units', 'activation' and 'rate' to specify the number of neurons in each layer, the activation function for the layer and the rate at which the neurons should be dropped.
}
\item{epochs_M}{
The number of epochs or iterations to train the model on. This only makes sense when the batch size is smaller than the number of records (which normally is).
}
\item{batch_size}{
A numeric value indicating the size of each batch to be used in the training process or the epochs.
}
\item{returnMatrices}{
A TRUE/FALSE statement to indicate if we should fit the model or just return the equations.
}

}
\details{

This function forms the mixed model equations developed by Henderson but uses the tensorflow and keras libraries to fit the predictive model using the number of layers specified by the user with their respective parameters such as number of units or neurons, the activation function, the rate of neuron dropping along with training model parameters such as number of iterations or epochs and batch size.

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{model}{a model fitted for prediction.}
}
}

\references{

Allaire J, Chollet F (2023). _keras: R Interface to 'Keras'_. R package version 2.13.0.9000

Allaire J, Tang Y (2023). _tensorflow: R Interface to 'TensorFlow'_. R package version 2.14.0.9000

}

\examples{
####=========================================####

# data("DT_example", package = "cgiarBase")
# idSta <- result$status$analysisId[2]
# iTrait <- "Plant_Height_cm"
# DT <- result$predictions[result$predictions$analysisId == idSta &
#                            result$predictions$trait == iTrait,]
# GT <- result$data$geno
# DT <- DT[which(DT$designation %in% rownames(GT)),]
# 
# #################################
# ## model with relationship matrix
# A <- sommer::A.mat(GT) 
# GTi=A # t(chol(A))
# ## model with markers
# GTi <- apply(GT, 2, sommer::imputev)
# `%>%` <- dplyr::`%>%`
# 
# DT2 <- cbind(DT, GTi[DT$designation,])
# ncol(DT); ncol(DT2)
# myLayers <- data.frame(units=c(100,100,100), activation="relu", rate=0.3);
# myLayers
# 
# ss <- sample(1:nrow(DT2), round(nrow(DT2)/2)) # divide in 2 samples
# TP <- DT2[ss,] # training population
# VP <- DT2[setdiff(1:nrow(DT2), ss),] # validation population
# dim(TP); dim(VP)
# # machine learning model using the grouping effect 
# model1 <- mler(fixed= predictedValue~1,
#                   random=~1,
#                   group=list(id=14:112),
#                   ginverse=NULL,
#                   data=TP, # model is trained with training population
#                   yourLayers=myLayers, epochs_M=50,
#                   batch_size=50, returnMatrices=FALSE)
# # matrices for the validation population
# mats <- mler(fixed= predictedValue~1,
#                   random=~1,
#                   ginverse=NULL,
#                   group=list(id=14:112),
#                   data=VP, # building the same matrices for the other set
#                   yourLayers=myLayers, epochs_M=50,
#                   batch_size=50, returnMatrices=TRUE)
# # produce predictions
# Yhat <-  predict(model1$model, mats$W) %>%
#   data.frame() %>%
#   setNames(colnames(mats$yvar))
# preds <- (Yhat * mats$sds) + mats$means
# # look at correlation between real value and prediction
# plot(preds$predictedValue, VP$predictedValue)
# cor(preds$predictedValue, VP$predictedValue)
# 
# # compare with reml
# library(sommer)
# M <- TP[,14:112]
# model2 <- mmer(predictedValue~1,
#                   random=~vsr(list(M)),
#                   data=TP)
# summary(model2)$varcomp
# u <- as.matrix(VP[,14:112]) %*% model2$U$`u:M`$predictedValue
# plot(u, VP$predictedValue)
# cor(u, VP$predictedValue)

}

