\name{tabToMat}
\alias{tabToMat}

\title{
Table to Matrix 
}
\description{
Converts a table a matrix format.
}
\usage{
tabToMat(Adf, symmetric=TRUE, V1="Var1n", V2="Var2n", Freq="Freq")
}

\arguments{
  \item{Adf}{
A table with coordinates to transform to a matrix.
}
\item{symmetric}{
Whether a symmetric matrix should be assumed or not.
}
\item{V1}{
Name of the numerical variable indicating rows of the original matrix.
}
\item{V2}{
Name of the numerical variable indicating column of the original matrix.
}
\item{Freq}{
Name of the numerical variable indicating value of the cell in the original matrix.
}

}
\details{

A simple function to a table into a matrix for data management purposes to return a numeric matrix.

}
\value{
The function returns a matrix.
}

\references{
None
}

\examples{

A <- matrix(1:9,3,3);A # original
B <- matToTab(A, symmetric = FALSE);B # move to data frame format
A2 <- tabToMat(B, symmetric = FALSE);A2 # move back to matrix

}