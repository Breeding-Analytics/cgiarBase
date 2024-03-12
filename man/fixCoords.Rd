\name{fixCoords}
\alias{fixCoords}

\title{
Fix coordinates information for miss-specified field trials
}
\description{
This function identifiees if the same row and column coordinate has been assigned to two different plots in a field instance or environment and splits them into a different replicate to avoid issues in the model fitting.
}
\usage{
fixCoords(mydataSub, rowcoord, colcoord, rep)
}

\arguments{

\item{mydataSub}{
Data set to correct with new field coordinates.
}
\item{rowcoord}{
Column in the dataset containing the row coordinate information.
}
\item{colcoord}{
Column in the dataset containing the range coordinate information.
}
\item{rep}{
Column in the dataset containing the replicate information.
}

}
\details{

This function identifiees if the same row and column coordinate has been assigned to two different plots in a field instance or environment and splits them into a different replicate to avoid issues in the model fitting.

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{result}{a corrected data set with new coordinates.}
}
}

\references{

None

}

\examples{
####=========================================####
data(DT_example)

dt <- result$data$pheno
dt$environment <- paste(dt$Env, dt$Experiment)
dtList <- split(dt, dt$environment)
# field view
lapply(dtList, function(x){with(x, table(Range, Pass))})
# mock up some data to add noise
sampleDt1 <- sampleDt2 <- dtList[[1]]
sampleDt2$Replicate <- sampleDt2$Replicate+2
sampleDt <- rbind( sampleDt1, sampleDt2 ) # introduce noise
with(sampleDt, table(Range, Pass, Env))
# fix the cordinates
dt2 <- fixCoords(mydataSub=sampleDt,rowcoord="Range",colcoord="Pass", rep="Replicate")
# look at the before and after
# lattice::levelplot(Plant_Height_cm~Range*Pass|Replicate,data=sampleDt)
# lattice::levelplot(Plant_Height_cm~Range*Pass|Replicate,data=dt2)

}

