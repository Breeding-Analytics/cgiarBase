\name{applyGenoModifications}
\alias{applyGenoModifications}

\title{
Modifications to a marker matrix
}
\description{
This function takes a table of modifications as established by the OneCGIAR Biometrical Genetics team and apply such changes to a marker matrix.
}
\usage{
applyGenoModifications(M, modifications)
}

\arguments{
  \item{M}{
Marker matrix numerically coded for allele dosage and with dimension n rows (individuals) and m columns (markers).
}
\item{modifications}{
A table of modifications as established by the OneCGIAR Biometrical Genetics team with columns: module, analysisId, reason, row, col, value. Please see example section.
}

}
\details{

This function takes the modifications table and goes row by row assigning new values to the marker matrix depending on the modification identified (imputation, deletion, etc.).

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{result}{a marker matrix modified.}
}
}

\references{

None

}

\examples{
####=========================================####
data(DT_example)
mods <- result$modifications$geno
M <- result$data$geno
M[12:15,20:24]
# apply modifications
M2 <- applyGenoModifications(M=M, modifications=mods)
M2[12:15,20:24] # check result

}

