\name{computeGenoModifications}
\alias{computeGenoModifications}

\title{
Marker modifications
}
\description{
This function computes a table of genetic marker modifications given the specified parameters.
}
\usage{
computeGenoModifications(
    M = NULL, 
    propNaUpperThreshForMarker=.3,
    propNaUpperThreshForInds=.3,
    propHetUpperThreshForMarker = 1, 
    propFisUpperThreshForMarker = 1,
    maf=0.05, ploidy=2,
    imputationMethod="median"
)
}

\arguments{
  \item{M}{
Marker matrix numerically coded for allele dosage and with dimension n rows (individuals) and m columns (markers).
}
\item{propNaUpperThreshForMarker}{
This sets a threshold for how much missing data in a marker is allowed. Any marker with more than this value will be marked as a column to be removed in posterior analyses. Value between 0 and 1.
}
\item{propNaUpperThreshForInds}{
This sets a threshold for how much missing data in an individual is allowed. Any individual with more than this value it will be marked as a row to be removed in posterior analyses. Value between 0 and 1.
}
\item{propHetUpperThreshForMarker}{
This sets an upper threshold for what is the maximum level of heterozygosity allowed in the markers. If greater than this value it will be marked as a column to be removed in posterior analyses. Value between 0 and 1. For example, a line dataset should not have markers with high heterozigosity.
}
\item{propFisUpperThreshForMarker}{
This sets an upper threshold for what is the maximum level of inbreeding allowed in the markers. If lower than this value it will be marked as a column to be removed in posterior analyses. Value between 0 and 1.
}
\item{maf}{
This sets a lower threshold for what is the minimum allele frequency allowed in the dataset. If lower than this value it will be marked as a column to be removed in posterior analyses. Value between 0 and 1.
}
\item{ploidy}{
Number of chromosome copies. This value is important to compute some of the paramters. Default is 2 or diploid.
}
\item{imputationMethod}{
Method to impute missing cells. Median is the only method currently available.
}

}
\details{

When genetic evaluation is carried using genomic data, we need to ensure the quality of genetic markers. This option aims to allow users to identify bad markers or individuals given certain QA parameters.

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{result}{a table of potential modifications given the specified thresholds for the parameters.}
}
}

\references{

None

}

\examples{
####=========================================####

data(DT_example)

M <- result$data$geno

mods <- computeGenoModifications(
  M = M, 
  propNaUpperThreshForMarker=.3,
  propNaUpperThreshForInds=.3,
  propHetUpperThreshForMarker = 1, 
  propFisUpperThreshForMarker = 1,
  maf=0.05, ploidy=2,
  imputationMethod="median"
)
head(mods)

}

