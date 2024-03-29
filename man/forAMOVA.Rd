\name{forAMOVA}
\alias{forAMOVA}

\title{
  forAMOVA
}
\description{
  For calculate the AMOVA.
}
\usage{
  forAMOVA(agc.dist,agc.env)
}

\arguments{
  \item{agc.dist}{Object of class "dist", represent the distance matrix.}
  \item{agc.env}{A dataframe of type of factor that contain the selected groups.}
  }

\details{
adonis is a function for the analysis and partitioning sums of squares using semimetric and metric distance matrices. Insofar as it partitions sums of squares of a multivariate data set, it is directly analogous to MANOVA (multivariate analysis of variance). M.J. Anderson (McArdle and Anderson 2001, Anderson 2001) refers to the method as “permutational manova” (formerly “nonparametric manova”). Further, as its inputs are linear predictors, and a response matrix of an arbitrary number of columns (2 to millions), it is a robust alternative to both parametric MANOVA and to ordination methods for describing how variation is attributed to different experimental treatments or uncontrolled covariates. It is also analogous to redundancy analysis (Legendre and Anderson 1999).

Typical uses of adonis include analysis of ecological community data (samples X species matrices) or genetic data where we might have a limited number of samples of individuals and thousands or millions of columns of gene expression data (e.g. Zapala and Schork 2006).

adonis is an alternative to AMOVA (nested analysis of molecular variance, Excoffier, Smouse, and Quattro, 1992; amova in the ade4 package) for both crossed and nested factors.
}

\value{
  The function returns a dataframe "variance" that is the AMOVA table where can find for each source of variance (Between Groups, Within Groups, Total)
  \describe{
    \item{$Df}{degree of freedom}
    \item{$SumsOfSqs}{Sums of squares.}
    \item{$Pr..F.}{p-value}
    \item{$F.Model}{pseudo F-statistic.}
    \item{$R2}{partial R2.}
    \item{$Sigma}{Components of variance.}
    \item{$PercVar}{Percent of variance explained.}
    \item{$Phi}{Provides the population differentiations statistics. We would expect a higher Phi value to represent a higher amount of differentiation.}
  }
}

\references{
Anderson, M.J. 2001. A new method for non-parametric multivariate analysis of variance. Austral Ecology, 26: 32–46.

Crawley, M.J. 2002. Statistical Computing: An Introduction to Data Analysis Using S-PLUS

Excoffier, L., P.E. Smouse, and J.M. Quattro. 1992. Analysis of molecular variance inferred from metric distances among DNA haplotypes: Application to human mitochondrial DNA restriction data. Genetics, 131:479–491.

Legendre, P. and M.J. Anderson. 1999. Distance-based redundancy analysis: Testing multispecies responses in multifactorial ecological experiments. Ecological Monographs, 69:1–24.

McArdle, B.H. and M.J. Anderson. 2001. Fitting multivariate models to community data: A comment on distance-based redundancy analysis. Ecology, 82: 290–297.

Warton, D.I., Wright, T.W., Wang, Y. 2012. Distance-based multivariate analyses confound location and dispersion effects. Methods in Ecology and Evolution, 3, 89–101.

Zapala, M.A. and N.J. Schork. 2006. Multivariate regression analysis of distance matrices for testing associations between gene expression patterns and related variables. Proceedings of the National Academy of Sciences, USA, 103:19430–19435.
}

\examples{

}
