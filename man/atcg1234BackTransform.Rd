\name{markerBackTransform}
\alias{markerBackTransform}
\title{
Letter to number converter
}
\description{
This function was designed to help users back transform the numeric marker matrices from the function atcg1234 into letters. 

}
\usage{
markerBackTransform(marks, refs)
}
\arguments{
  \item{marks}{a centered marker matrix coming from atcg1234 (not centered).}
  \item{refs}{a 2 x m matrix for m markers (columns) and 2 rows where the reference and alternate alleles for each marker are indicated.}
  
}
\value{
\describe{

\item{markers}{ a new marker matrix leter coded according to the reference allele matrix.}

}
}

\author{
Giovanny Covarrubias-Pazaran
}

\examples{

data(DT_example)

refs <- t(result$metadata$geno[,c("refAllele", "altAllele")])
rownames(refs) <- c("Ref","Alt")
# letter coded
numob <- markerBackTransform(marks =  result$data$geno+1, refs = refs   )
numob[1:4,1:4]
# numeric coded
bad <- apply(numob,2,function(x){length(which(x==""))/length(x)})
numoc <- atcg1234(numob[,which(bad < 1)], maf = -1  )
numoc$M[1:4,1:4]


}