\name{grm2}
\alias{grm2}

\title{
Genomic Relationship Matrix
}
\description{
Wrapper to calculate the genomic relationship matrix while accounting for different ploidies and missing data present in the marker set.
}
\usage{
grm2( M= NULL, verbose=FALSE )
}

\arguments{

  \item{M}{
Marker matrix numerically coded for allele dosage and with dimension n rows (individuals) and m columns (markers).
}
\item{verbose}{
A TRUE/FALSE statement indicating if logs and warning messages should be printed or not.
}

}
\details{

Wrapper to calculate the genomic relationship matrix while accounting for different ploidies and missing data present in the marker set

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{result}{a genomic relationship matrix.}
}
}

\references{

None

}

\examples{
####=========================================####
data(DT_example)
# make a relationship matrix
AGE <- grm2(result$data$geno)
AGE[1:4,1:4]
}

