\name{geneticEvaluation3}
\alias{geneticEvaluation3}

\title{
Genetic evaluation
}
\description{
This function forms the mixed model equations developed by Henderson and solves the equations given the specified variance components.
}
\usage{
geneticEvaluation3(fixed=NULL,random=NULL,rcov=NULL,weights=NULL, 
                   AFGE=NULL, varComp=NULL, data=NULL, keep=NULL, 
                   tolParInv=1e-5)
}

\arguments{
  \item{fixed}{
Fixed formula for the desired model fitting.
}
\item{random}{
Random formula for the desired model fitting.
}
\item{rcov}{
Residual formula for the desired model fitting.
}
\item{weights}{
Name of the column in the dataset containing the weights vector.
}
\item{AFGE}{
Additive relationship matrix for the individual to predict. It should include the individuals in the dataset plus the individuals to be predicted.
}
\item{varComp}{
List of variance components for the random and residual factors.
}
\item{data}{
Input dataset.
}
\item{keep}{
A numeric value indicating which individuals should not be returned in the predictions.
}
\item{tolParInv}{
Numeric value for the inversion of the coefficient matrix.
}

}
\details{

This function forms the mixed model equations develpoed by Henderson and solves the equations given the specified variance components. It uses de mmec function from sommer to build the model matrices and then we form the coeficient matrix for solveing.

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{result}{a list of predicted values and standard errors for individuals.}
}
}

\references{

None

}

\examples{
####=========================================####

data(DT_example)

# extract data
result$status
mydataSub <- result$predictions
mydataSub <- mydataSub[mydataSub$analysisId == result$status$analysisId[2], ]
mydataSub <- mydataSub[mydataSub$trait == "Plant_Height_cm",]
mydataSub$w <- 1/(mydataSub$stdError^2)
nrow(mydataSub)
# make a relationship matrix
AGE <- grm2(result$data$geno)
# make variance component table
VarDf <- data.frame(Variance=c(0.5,1), VarComp=c("designation","residual")); VarDf
# make the genetic evalaution assuming some variance components
head(mydataSub)
final <- geneticEvaluation3(
  fixed =as.formula("predictedValue~1"),
  random = as.formula("~designation"),
  rcov=NULL,
  weights = "w",
  data = mydataSub,
  varComp=VarDf,
  keep=NULL,
  AFGE=AGE
)
str(final)

}

