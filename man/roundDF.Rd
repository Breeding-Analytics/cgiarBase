\name{roundDF}
\alias{roundDF}

\title{
Rounding of numbers for data frames
}
\description{
This function rounds numbers in a data frame with a specified number of singificant digits.
}
\usage{
roundDF(dataFrame, numRound)
}

\arguments{
  \item{dataFrame}{A data frame.}
  \item{numRound}{Number of digits for rounding.}
}

\details{
It ensures that the smalles number in each column has at least \code{numRound} significant digits.
}

\value{
It returns a data frame with the rounded numbers.
}

\references{
None
}

\examples{
####=========================================####
data(DT_example)

roundDF(result$data$pheno, 3)

}

