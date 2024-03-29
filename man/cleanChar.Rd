\name{cleanChar}
\alias{cleanChar}

\title{
Clean names from special characters
}
\description{
This function simply removes special characters from a name.
}
\usage{
cleanChar(x)
}

\arguments{

  \item{x}{
Character vector with names to be cleaned from special characters.
}

}
\details{

This function simply removes special characters from a name. Intended to be used to clean names of the designation labels.

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{result}{a vector of names without special characters.}
}
}

\references{

None

}

\examples{
####=========================================####
data(DT_example)

hybNames <- unique(result$data$pheno$Hybrid)
# example 1
cleanChar(hybNames)
# example 2
cleanChar("weIRDo?:h")

}

