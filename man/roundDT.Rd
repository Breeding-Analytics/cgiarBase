\name{roundDT}
\alias{roundDT}

\title{
Rounding of numbers for data tables
}
\description{
This function rounds numbers in a data table with a specified number of singificant digits.
}
\usage{
roundDT(dataTable, numRound)
}

\arguments{
  \item{dataTable}{A data table.}
  \item{numRound}{Number of digits for rounding.}
}
\details{
It ensures that all values have at least \code{numRound} significant digits.
}
\value{
It returns a data table with the rounded numbers.
}

\references{
None
}

\examples{
####=========================================####
data(DT_example)

dataTable <- DT::datatable(head(result$data$pheno))

roundDT(dataTable, 3)

}