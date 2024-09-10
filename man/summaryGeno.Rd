\name{summaryGeno}
\alias{summaryGeno}

\title{
Summary of a modifications to a marker matrix
}
\description{
This function takes a table of modifications as established by the OneCGIAR Biometrical Genetics team and summarizes such changes to a marker matrix.
}
\usage{
summaryGeno(object, analysisId)
}

\arguments{
  \item{object}{
The data object structure produced from bioflow which contains data, metadata, modifications and tables of results (metrics, predictions, modeling, status), in order to pick the necessary information to calculate an index.
}
\item{analysisId}{
A numeric id to link to the predictions to be used for the analysis.
}

}
\details{

This function takes the modifications table and gives some summary statistics (imputation, deletion, etc.).

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{result}{a summary table.}
}
}

\references{

None

}

\examples{
####=========================================####
data(DT_example)
res <- summaryGeno(result, result$status$analysisId[1])
res

}

