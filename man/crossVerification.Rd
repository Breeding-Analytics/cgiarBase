\name{crossVerification}
\alias{crossVerification}

\title{
Verification of a cross using markers
}
\description{
This function takes a matrix of genotypes for males, another for females, and another for progeny and calculates the probability of matching the expected genotype. Additional metrics like hetrozygosity and deviation from expectations are returned.
}
\usage{
crossVerification(Mf,Mm,Mp, Mexp=NULL, ploidy=2)
}

\arguments{
  \item{Mf}{
Marker matrix numerically coded for allele dosage (0,...,n where n is the ploidy of the species) and with dimension n rows (mother) and m columns (markers).
}
  \item{Mm}{
Marker matrix numerically coded for allele dosage (0,...,n where n is the ploidy of the species) and with dimension n rows (fathers) and m columns (markers).
}
  \item{Mp}{
Marker matrix numerically coded for allele dosage (0,...,n where n is the ploidy of the species) and with dimension n rows (progeny) and m columns (markers).
}
  \item{Mexp}{
Marker matrix numerically coded for allele dosage (0,...,n where n is the ploidy of the species) and with dimension n rows (individuals) and m columns (markers) of the expected genotype. If NULL it will be assumed that the expected genotype is the cross between males and females:

Mexp = (Mf + Mm)/2

which is the F1 cross of two parents. Alternatively the user can provide a different expected genotype. For example, in back crossing applications this could be the recurrent parent.

}
\item{ploidy}{
A numeric value indicating the ploidy of the individuals.
}

}
\details{

This function takes the markers matrices and performs the cross verification.

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{metrics}{a table of individuals' metrics.}
\item{matchMat}{the matrix of probabilities.}
}
}

\references{

None

}

\examples{
####=========================================####
data(DT_example)

M <- result$data$geno
Mf <- M[1,,drop=FALSE]; Mm <- M[2,,drop=FALSE]; Mp <- M[3,,drop=FALSE];
Mf <- M[sample(1:nrow(M)),]; Mm <- M[sample(1:nrow(M)),]; Mp <- M[sample(1:nrow(M)),]
ss <- crossVerification(Mf,Mm,Mp, Mexp=NULL, ploidy=2)
head(ss$metrics)
Matrix::image(ss$matchMat[1:40,1:40])

}

