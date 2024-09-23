\name{summaryData}
\alias{summaryData}

\title{
Summary of data records in bioflow object
}
\description{
This function takes a table of number of records and features for each data element in the bioflow data structure.
}
\usage{
summaryData(object)
}

\arguments{
  \item{object}{
The data object structure produced from bioflow which contains data, metadata, modifications and tables of results (metrics, predictions, modeling, status), in order to pick the necessary information to calculate an index.
}

}
\details{

This function takes the bioflow object and summarizes which data types are available.

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{res}{a summary table.}
}
}

\references{

None

}

\examples{
####=========================================####
data(DT_example)
res <- summaryData(result)
res

}

