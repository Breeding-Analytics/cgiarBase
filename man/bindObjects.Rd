\name{bindObjects}
\alias{bindObjects}

\title{
Binding bioflow objects
}
\description{
This function takes two bioflow objects as established by the OneCGIAR Biometrical Genetics team and merges the data in the best possible way.
}
\usage{
bindObjects(
    object1= NULL,
    object2= NULL,
    verbose=FALSE
)
}

\arguments{
  \item{object1}{
List object with the different data slots specified by the OneCGIAR Biometrical Genetics team: status, predictions, modeling, metrics, data, metadata and modifications.
}
  \item{object2}{
List object with the different data slots specified by the OneCGIAR Biometrical Genetics team: status, predictions, modeling, metrics, data, metadata and modifications.
}
\item{verbose}{
A TRUE/FALSE statement indicating if logs and warning messages should be printed or not.
}

}
\details{

This function takes two bioflow objects and merges them. In the case of phenotype data it will row bind the columns in common and keep separated the different columns. In the case of genotype data it will keep all data from both sets (augmented) and in case individuals are repeated it will keep the one with the highest marker density.

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{result}{a new bioflow object merging the two provided.}
}
}

\references{

None

}

\examples{
####=========================================####
data(DT_example)

obj1 <- result
obj2 <- result

obj3 <- bindObjects(
    object1= obj1,
    object2= obj2,
    verbose=FALSE
)
head(obj3$status)

}

