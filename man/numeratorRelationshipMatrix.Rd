\name{nrm2}
\alias{nrm2}

\title{
Numerator Relationship Matrix
}
\description{
Calculates a numerator relationship matrix for a data with 3 columns indicating sire, dam and progeny.
}
\usage{
nrm2(
    pedData= NULL,
    verbose=FALSE,
    indivCol=NULL, damCol=NULL, sireCol=NULL,
    returnMatrix=TRUE
)
}

\arguments{

  \item{pedData}{
Data frame with 3 columns indicating father, mother and progeny.
}
\item{verbose}{
A TRUE/FALSE statement indicating if logs and warning messages should be printed or not.
}
\item{indivCol}{
A character vector with one value indicating the column where progeny names can be found in the dataset.
}
\item{damCol}{
A character vector with one value indicating the column where the mother names can be found in the dataset.
}
\item{sireCol}{
A character vector with one value indicating the column where the father names can be found in the dataset.
}
\item{returnMatrix}{
A TRUE/FALSE statement indicating if the relationship matrix should be returned. If FALSE the curated pedigree is returned instead.
}

}
\details{

Wrapper of pedigreem package functions to clean and shape a pedigree file to posteriorly calculate the numerator relationship matrix using efficient methods.

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{pedigree}{a numerator relationship matrix or curated pedigree.}
}
}

\references{

None

}

\examples{
####=========================================####
data(DT_example)

dt <- result$data$pedigree
head(dt)
N <- nrm2(
    pedData= dt,
    verbose=FALSE,
    indivCol="Hybrid", damCol="Hybrid_Parent1", sireCol="Hybrid_Parent2",
    returnMatrix=TRUE
)
N[1:4,1:4]
# Matrix::image(as(N, Class = "dgCMatrix"))

}

