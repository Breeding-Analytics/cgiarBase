\name{mlModel}
\alias{mlModel}

\title{
Machine learning model
}
\description{
This function forms the linear equations to fit a machine learning model based on back propagation using the keras and tensorflow libraries.
}
\usage{
mlModel(fixed, 
        random=NULL, ginverse=NULL,
        group, data,
        markerEffects=TRUE,
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
\item{markerEffects}{
A TRUE/FALSE statement to indicate if the model should use the group or the ginverse argument to create the incidence matrices.
}
\item{yourLayers}{
A data frame with as many rows as layers of neurons to be fitted. The data frame should contain 3 columns named 'units', 'activation' and 'rate' to specify the number of neurons in the layer, the activation function for the layer and the rate at which the neurons should be dropped.
}
\item{epochs_M}{
The number of epochs or iterations to train the model on. This only makes sense when the batch size is smaller than the number of records (which normally is).
}
\item{batch_size}{
A numeric value indicating the size of each batch to be used in the training process.
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
# myLayers <- data.frame(units=c(16,10), activation="sigmoid", rate=0.3);
# myLayers
# model1 <- mlModel(fixed= predictedValue~1, 
#                   random=~designation, ginverse=list(id=A),
#                   # group, 
#                   data=DT,
#                   markerEffects=FALSE,
#                   yourLayers=myLayers, epochs_M=50, 
#                   batcheSize=50, returnMatrices=FALSE)
# 
# #################################
# ## model with markers
# DT2 <- cbind(DT, GT[DT$designation,])
# ncol(DT); ncol(DT2)
# myLayers <- data.frame(units=c(16,10), activation="sigmoid", rate=0.3);
# myLayers
# model1 <- mlModel(fixed= predictedValue~1, 
#                   random=~designation, #ginverse=list(id=A),
#                   group=list(id=14:112),
#                   data=DT,
#                   markerEffects=FALSE,
#                   yourLayers=myLayers, epochs_M=50, 
#                   batcheSize=50, returnMatrices=FALSE)

}

