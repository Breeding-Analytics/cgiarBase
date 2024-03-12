\name{matToTab}
\alias{matToTab}

\title{
Matrix to Table 
}
\description{
Converts a matrix to a table format.
}
\usage{
matToTab(A, symmetric=TRUE)
}

\arguments{
  \item{A}{
A numeric matrix.
}
\item{symmetric}{
Whether a symmetric matrix should be assumed or not
}

}
\details{

A simple function to convert a matrix into a table for data management purposes to return coordinates of data information in numeric and character format.

}
\value{
The function returns a data frame with coordinates of data information.
}

\references{
None
}

\examples{

A <- matrix(1:9,3,3);A # original
B <- matToTab(A, symmetric = FALSE);B # move to data frame format
A2 <- tabToMat(B, symmetric = FALSE);A2 # move back to matrix

}