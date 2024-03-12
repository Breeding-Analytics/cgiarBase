\name{asremlFormula}
\alias{asremlFormula}

\title{
Factor identification for model fitting
}
\description{
This function is the core of the single trial analysis pipeline function in the cgiarPipeline package. It is in charge of identifying which factors have enough levels to be fitted to come up with the optimal model to fit.
}
\usage{
asremlFormula(fixed, random, rcov, dat, minRandomLevels=NULL, minResidualLevels=NULL, 
              exchangeRandomEffects=NULL, exchangeResidualEffects=NULL,customRandomLevels=NULL, 
              customResidualLevels=NULL, xCoordinate= "ROWf",yCoordinate ="RANGEf", 
              doubleConstraintRandom=c("ROWf","RANGEf"), verbose=FALSE)
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
\item{dat}{
Input dataset.
}
\item{minRandomLevels}{
Number of levels to be considered a minimum to allow the factor to get in the random formula.
}
\item{minResidualLevels}{
Number of levels to be considered a minimum to allow the factor to get in the residual formula. This is provided by each factor (e.g. row=2, col=4, ...).
}
\item{exchangeRandomEffects}{
List of random effects that should be exchanged to meet the condition of minRandomLevels. For example, if user specifies that row and range should be exchnaged it means that both factors should meet the minRandomLevels criteria of each other.
}
\item{exchangeResidualEffects}{
List of residual effects that should be exchanged to meet the condition of minRandomLevels. For example, if user specifies that row and range should be exchnaged it means that both factors should meet the minResidualLevels criteria of each other.
}
\item{customRandomLevels}{
Other customized random levels.
}
\item{customResidualLevels}{
Other customized residual levels.
}
\item{xCoordinate}{
Name of the column that indicates the x-coordinate.
}
\item{yCoordinate}{
Name of the column that indicates the y-coordinate.
}
\item{doubleConstraintRandom}{
Pair of factors that should meet the constraint of each other.
}
\item{verbose}{
Should logs and warning messages be printed?
}

}
\details{

This function is in charge of identifying which factors have enough levels to be fitted to come up with the optimal model to fit. This aims to create the optimal formula for each field for either a one or two-stage approach.

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{result}{a list with the best formula and factors to be fitted per field.}
}
}

\references{

None

}

\examples{
####=========================================####

data(DT_example)

# extract data
mydataSub <- result$data$pheno
mydataSub[,"trait"] <- mydataSub$Yield_Mg_ha
# replace names
colnames(mydataSub) <- replaceValues(Source = colnames(mydataSub) , 
                                     Search = c("Env","Experiment","Replicate", "Block" ,"Range","Pass"),
                                     Replace = c("environment","trial","rep","iBlock","row","col") )
for(iFactor in c("environment","trial","rep","iBlock","row","col")){
  mydataSub[,paste0(iFactor,"F")] <- as.factor(mydataSub[,iFactor])
}

mde <- cgiarBase::asremlFormula(fixed=as.formula(paste("trait","~ 1")),
                                random=~ at(environmentF):rowF + at(environmentF):colF + at(environmentF):trialF + at(environmentF):repF + at(environmentF):iBlockF,
                                rcov=~at(environmentF):id(rowF):id(colF),
                                dat=droplevels(mydataSub[which(!is.na(mydataSub[,"trait"])),]),
                                
                                minRandomLevels=list(rowF= 3, colF=3, trialF=2,repF=2, blockF=4),
                                minResidualLevels=list(rowF=5, colF=5),
                                
                                exchangeRandomEffects=list(rowF="colF", colF="rowF"),
                                
                                exchangeResidualEffects=list(rowF="colF", colF="rowF"),
                                
                                customRandomLevels=NULL, customResidualLevels=NULL,
                                
                                xCoordinate= "rowF",yCoordinate ="colF",
                                doubleConstraintRandom=c("rowF","colF"), verbose=TRUE)

factorsFitted <- unlist(lapply(mde$used$environmentF,length))
factorsFitted
factorsFittedGreater <- which(factorsFitted > 0)
factorsFittedGreater

}

