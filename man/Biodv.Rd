\name{Biodv}
\alias{Biodv}

\title{
  Biodv
}
\description{
  Calculate heterozygosity (observed and expected), shannon index, number of effective allele, specificity, rareness,  percent of missing values, percent of polymorphic loci, Rogers distance, Nei distance, cluster analysis and multidimensional scaling.
}
\usage{
  Biodv(file_name,datos,nall,distk,mayorque,menorque,missval)
}

\arguments{
  \item{file_name}{
    The name for the output files (character).
  }
  \item{datos}{
    A dataframe in which columns are genotype information and in row marker information.
  }
  \item{nall}{
    Numeric value for identify the number of alleles, by default is 2.
  }
  \item{distk}{
    "Rogers" or "Nei" for identify which distance will be calculated.
  }
  \item{mayorque}{
    A number between 0-1, upper limit for MAF.
  }
  \item{menorque}{
    A number between 0-1, lower limit for MAF.
  }
  \item{missval}{
    A number between 0-1, delete genotypes that have more percent of missing values that the threshold
  }
}

\details{

\bold{Percent of polymorphic loci}

A gene is defined as polymorphic if the frequency one of its alleles is less than or equal to 0.95 or 0.99
\eqn{P_j=q\le 0.95} or \eqn{P_j=q\le 0.99} where \eqn{P_j} is the polymorphic rate and \eqn{q} is the frequency allele. This measure provides the criteria to determine whether a gene has variation.

From here we can calculate the proportion of polymorphic loci:
\deqn{P= \frac{n_{P_j}} {n_{total}}}

where \eqn{n_{P_j}} is the number of polymorphic loci and \eqn{n_{total}} is the total number of loci. Expresses the percentage of loci in a population variables.

\bold{Number of effective allele}

The number of alleles that may be present in a population

\deqn{A_e=\sum_{l=1}^L\frac{1}{1-h_l}=\sum_{l=1}^L\frac{1}{\sum p_i^2}}

Where \eqn{p_i} is the frequency in the \eqn{ith} allele in one locus and \eqn{h_l=1-\sum p_i^2} is the heterozygosity in a locus. This measure of diversity can provide useful information for establishing collection strategies. For example, we estimate the number effective alleles in a sample. Then the check in a or different shows throughout the collection. If the figure obtained the second time is less than the first, this could mean that our strategy collection needs revision.

\bold{Expected Heterozygosity}

Genetic diversity of Nei is the probability that, in a single locus, any pair of alleles, chosen at random population, different from each other. It can be calculated expected heterozygosity, which is the average of all loci is a estimate the degree of genetic variability in
the population

\deqn{H_e=\frac{1}{L} \sum_{j=1}^L{h_j} }

where \eqn{h_j} is the heterozygosity per locus and \eqn{L} is the total number of loci. The values for \eqn{H_e} are between 0 and 1 and a minimum of 30 loci to be analyzed in 20 individuals per population, to reduce the risk of statistical bias.

\bold{Shannon's Index}

Shannon's index accounts for both abundance and evenness of the species present.
\deqn{SH=-\sum_{a=1}^A \hat{p}_a log(\hat{p}_a) }

where \eqn{\hat{p}_a} is the estimated frequency of the allele \eqn{a} on the whole sample and \eqn{A} is the total number of alleles in the sample.

\bold{Genetic distance}

\bold{Modified Rogers distance}
\deqn{MR_{xy}=\sqrt{\frac{\sum_{l=1}^L \sum_{a=1}^{n_l} (p_{lax}-p_{lay})^2}{2L}} }

where p_{lax} is the estimated frequency of the allele a, within the locus l, at the genotype x; L the number of loci, and n_l the number of alleles within the lth locus.

\bold{Nei distance}

\deqn{NeiD_{xy}=-ln\left( \frac{\sum_{l=1}^L \sum_{a=1}^{n_l} p_{lax}p_{lay}}{\sum_{l=1}^L \sum_{a=1}^{n_l} p_{lax}^2 \sum_{l=1}^L \sum_{a=1}^{n_l} p_{lay}^2} \right) }

where p_{lax} is the estimated frequency of the allele a, within the locus l, at the genotype x; L the number of loci, and n_l the number of alleles within the lth locus.

\bold{Specificity of a marker in each allele}
The so-defined mutual information equals the average allele specificity, defined in this context as the information gained about an accession's identity, by the random extraction and identification of a the allele.

\deqn{S_i=\sum_{j=1}^N {\frac{p_{ij}}{Np_i}log_2 \left(\frac{p_{ij}}{p_i} \right)} }

where:

\itemize{
\item N: is the total number of accessions.
\item p_{ij}: the allele frequency of i-th allele within the accession j.
\item p_i=\frac{1}{N}\sum_{j=1}^N p_{ij}: the average frequency of the i-th allele across accessions.
}

\bold{Rarity of an accession}

The rarity of an accession is defined as the average specificity of the alleles it contains.

\deqn{R_j=\sum_{i=1}^k p_{ij}S_i}

where:

\itemize{
\item k: is the total number of alleles, in our case always is 2.
}


\bold{Multidimensional scaling analysis (MDS)}

After calculating the distance matrix proceeds to represent the distances among the objects in a parsimonious (and visual) way (i.e., a lower k-dimensional space). The goal of an MDS analysis is to find a spatial configuration of objects when all that is known is some measure of their general (dis)similarity. The spatial configuration should provide some insight into how the subject(s) evaluate the stimuli in terms of a (small) number of potentially unknown dimensions.

Classical MDS algorithms typically involve some linear algebra. The classical MDS algorithm rests on the fact that the coordinate matrix X can be derived by eigenvalue decomposition from the scalar product matrix B = XX'.

The problem of constructing B from the proximity matrix P is solved by multiplying the squared proximities with the matrix J=I-n^{-1}11'. This procedure is called double centering. The following steps summarize the algorithm of classical MDS:

\itemize{
\item Set up the matrix of squared proximities P^{(2)}.
\item Apply the double centering: B=\frac{-1}{2}JP^{(2)}J using the matrix J=I-n^{-1}11', where n is the number of objects.
\item Extract the m largest positive eigenvalues \lambda_1 \cdots \lambda_m of B and the corresponding m eigenvectors e_1 \cdots e_m.
\item A m-dimensional spatial configuration of the n objects is derived from the coordinate matrix X=E_m\Lambda^{1/2}_m , where E_m is the matrix of m eigenvectors and \Lambda_m is the diagonal matrix of m eigenvalues of B, respectively.
}

\bold{Cluster Analysis}

The best numerical classification strategy is the one that produces the most compact and well separated groups, i.e.,minimum variability
within each group and maximum variability among groups.

\bold{Ward's Method}
In Ward's minimum variance method, the distance between two clusters is the ANOVA sum of squares between the two clusters added up over all the variables. At each generation, the within-cluster sum of squares is minimized over all partitions obtainable by merging two clusters from the previous generation. The sums of squares are easier to interpret when they are divided by the total sum of squares to give the proportions of variance (squared semipartial correlations).

Ward's method joins clusters to maximize the likelihood at each level of the hierarchy under the assumptions of multivariate normal mixtures, spherical covariance matrices, and equal sampling probabilities.

Ward's method tends to join clusters with a small number of observations and is strongly biased toward producing clusters with approximately the same number of observations. It is also very sensitive to outliers.

Distance for Ward's method is:

\deqn{D_{KL}=\frac{\parallel \bar{x}_K-\bar{x}_L \parallel ^2}{\frac{1}{N_K}+\frac{1}{N_L}} }

where \bar{x}_K and \bar{x}_L is the mean vector for the Kth and Lth cluster respectively; N_K and N_L is the number of observation in the Kth and Lth cluster respectively and \parallel x \parallel is the square root of the sum of the squares of the elements of x.
}

\value{
The function returns a list containing
\describe{
\item{$div}{A dataframe with information about percent of polymorphic loci, expected heterozygosity, standar desviation for HE, observed heterozygosity, standar desviation for HO, number of effective allele, standar desviation for Ae, shannon diversity index, standar desviation for shannon index"}
\item{$coord2}{A dataframe with information about 3 first factors from the multidimensional scaling analysis}
\item{$dirsave}{The work directory where we are working}
\item{$clust}{Object of class "agnes", represent an agglomerative hierarchical clustering}
\item{$datos}{A dataframe in which columns are genotype information and in row marker information.}
\item{$mrdMAT}{Object of class "dist", represent the distance matrix}
\item{$perctCP12}{Percent of variability explained by the first two factors from the multidimensional scaling analysis}
\item{$OtherResults$Repeated_Markers}{A vector with the names of markers that the sequence of values are the same with others}
\item{$OtherResults$Repeated_Genotypes}{A vector with the names of genotypes that the sequence of values are the same with others}
\item{$OtherResults$MarkOutFilterPoly}{A vector with the names of markers that results polimorphics}
\item{$OtherResults$CalculusPerMarker}{A dataframe with information about heterozygosity (expected and observed), number of effective allele, shannon index, specificity and percent of missing values for each marker.}
\item{$OtherResults$CalculusPerGenotype}{A dataframe with information about heterozygosity (expected and observed), number of effective allele, shannon index, rareness and percent of missing values for each genotype.}
}
}

\references{
de Vicente, M.C., Lopez, C. y Fulton, T. (eds.). 2004. Analisis de la Diversidad Genetica Utilizando
Datos de Marcadores Moleculares: Modulo de Aprendizaje. Instituto Internacional
de Recursos Fitogeneticos (IPGRI), Roma, Italia.
}

\examples{



}
