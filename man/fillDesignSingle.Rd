\name{fillDesignSingle}
\alias{fillDesignSingle}

\title{
Fill field-gaps in row and column
}
\description{
This function adds the missing rows and column records for a field/environment to be able to perform and AR1 model.
}
\usage{
fillDesignSingle(x,rows="ROW",ranges="RANGE", by, extra)
}

\arguments{

\item{x}{
Data set to enhance with new field coordinates.
}
\item{rows}{
Column in the dataset containing the row coordinate information.
}
\item{ranges}{
Column in the dataset containing the range coordinate information.
}
\item{by}{
Column in the dataset containing the environment or field instance information.
}
\item{extra}{
Column(s) in the dataset containing extra nested fields information.
}

}
\details{

This function adds the missing rows and column records for a field/environment to be able to perform and AR1 model.

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{result}{an enhanced data set with missing coordinates.}
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
sampleDt <- dtList[[1]]
# lattice::levelplot(Plant_Height_cm~Range*Pass,data=sampleDt)
head(sampleDt)
# create some gaps
gaps <- c(4,6,8)
deleted <- sampleDt[gaps, c("Replicate", "Block", "Plot", "Range", "Pass")];deleted
sampleDt <- sampleDt[-gaps,]
# fill design
dt2 <- fillDesignSingle(sampleDt,rows="Range",ranges="Pass", by="Env")
# gaps are filled
dt2[gaps, c("Replicate", "Block", "Plot", "Range", "Pass")]

}

