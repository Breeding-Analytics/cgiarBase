\name{gdiv}
\alias{gdiv}

\title{
  gdiv
}
\description{
  For calculate the Fst among groups and within each group.
}
\usage{
 gdiv(datos,datos1,selv,quitomono,mdrMAT)
}

\arguments{
  \item{datos}{A dataframe in which columns are genotype information and in row marker information.}
  \item{datos1}{A dataframe with information about 3 first factors from the multidimensional scaling analysis for each designation and the group to which it belongs in the cluster or the passport information.}
  \item{selv}{A character with the name of the column that identifies which grouping we will analyze.}
  \item{quitomono}{TRUE or FALSE. If TRUE for each group delete the monomorphic markers.}
  \item{mdrMAT}{Object of class "dist", represent the distance matrix}
}

\details{
\bold{Diversity within groups}

\deqn{H_{sl}=1-\sum p_s_i^2 }

where \eqn{p_s_i} is the frequency in the \eqn{ith} allele in one locus in the \eqn{sth} subpopulation and if you calculate the mean of this you obtain:

\deqn{H_s=\frac{1}{L} \sum_{j=1}^L{H_{slj}} }

\bold{Diversity among groups (Wright's statistics)}

\deqn{F_{ST}=\frac{D_{ST}}{H_e}}

where \eqn{D_{ST}=H_e-H_s} is the diversity among individuals within subpopulation and \eqn{H_e} is the expected heterozygosity, so, \eqn{F_{ST}} measure the degree of genetic differentiation among populations, depending of allele frequencies.

Then we say that if \eqn{F_{ST}} is 0-0.05 the genetic differentiation is small.

Then we say that if \eqn{F_{ST}} is 0.05-0.15 the genetic differentiation is middle.

Then we say that if \eqn{F_{ST}} is 0.15-0.25 the genetic differentiation is big.

Then we say that if \eqn{F_{ST}} is >0.25 the genetic differentiation is very big.

}

\value{
  The function returns a list
  \describe{
    \item{$div2}{A dataframe with the values for the Fst among groups and within each group that you was analyzed, also give us the total number of markers which was used for do the calculus.}
    \item{$writemono}{If quitomono is TRUE, output a dataframe with information about which markers in ecah group was deleted because are monomorphics.}
  }
}

\references{
de Vicente, M.C., Lopez, C. y Fulton, T. (eds.). 2004. Analisis de la Diversidad Genetica Utilizando
Datos de Marcadores Moleculares: Modulo de Aprendizaje. Instituto Internacional
de Recursos Fitogeneticos (IPGRI), Roma, Italia.
}

\examples{

}
