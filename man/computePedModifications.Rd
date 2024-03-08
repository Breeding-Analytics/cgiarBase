\name{computePedModifications}
\alias{computePedModifications}

\title{
  computePedModifications
}
\description{
  For do the quality control using the pedigree and genetic information.
}
\usage{
  computePedModifications(M = NULL, Tb1=NULL, Tb2=NULL, Limp1=14000,
  Limp2=17000, QCc1=1.2, QCc2=2, GmatFL1=0.2,  GmatFL2=1.6, GmatFL3=0.7,
  GmatML1=0.2, GmatML2=0.7, GmatD_FL1=0.2, HetL1=5, HetF1=5, HetM1=5,
  Fimp1=16000, Mimp1=16000, QcTest=1000, QcTestG1=1400, QcTestG2=1200)
}

\arguments{
  \item{M}{
    Object of class "Matrix", genetic information in numeric format.
  }
  \item{Tb1}{
    A dataframe with 3 columns (Line, Female and Male) that correspond to pedigree information.
  }
  \item{Tb2}{
    A character url of hapmap file.
  }
  \item{Limp1}{A numeric value for assign the lower threshold for designation imputed. The default value is 14000.}
  \item{Limp2}{A numeric value for assign the upper threshold for designation imputed. The default value is 17000.}
  \item{Fimp1}{A numeric value for assign the threshold for Father imputed. The default value is 16000.}
  \item{Mimp1}{A numeric value for assign the threshold for Mother imputed. The default value is 16000.}
  \item{QCc1}{A numeric value for assign the lower threshold for allele quality control (QC pass). The default value is 1.2.}
  \item{QCc2}{A numeric value for assign the upper threshold for allele quality control (QC fail). The default value is 2.}
  \item{QcTest}{A numeric value for assign the lower threshold for allele quality control (QC pass).The default value is 1000.}
  \item{QcTestG1}{A numeric value for define group 1 with the condition that designation imputed > Qc G1 (No assessment-designation sample too poor).The default value is 1400.}
  \item{QcTestG2}{A numeric value for define group 2 with next conditions; designation imputed <= Qc G1 & Designation mother>Qc G2 & Designation father>Qc G2 (No assessment mother and father samples too poor), Mother imputed>Father imputed (No QC test, no mother assessment - mother sample too poor),Mother imputed<Father imputed (No QC test, no father assessment - father sample too poor).The default value is 1200.}
  \item{GmatFL1}{A numeric value for assign the lower threshold for G matrix value comparison between designation and mother, and conclude if the mother is wrong. The default value is 0.2.}
  \item{GmatFL2}{A numeric value for assign the upper threshold for G matrix value comparison between designation and mother, and conclude if the mother is correct. The default value is 1.6.}
  \item{GmatFL3}{A numeric value for assign the intermediate threshold for G matrix value comparison between designation and mother, and conclude if the mother is in the borderline. The default value is 0.7.}
  \item{GmatML1}{A numeric value for assign the lower threshold for G matrix value comparison between designation and father, and conclude if the father is wrong. The default value is 0.2.}
  \item{GmatML2}{A numeric value for assign the intermediate threshold for G matrix value comparison between designation and father, and conclude if the father is in the borderline. The default value is 0.7.}
  \item{GmatD_FL1}{A numeric value for assign the lower threshold for G matrix value comparison between diagonal designation and designation-mother, and conclude if the mother is correct or probable self. The default value is 0.2.}
  \item{HetL1}{A numeric value for assign the threshold for heterozygosity in designation.The default value is 5.}
  \item{HetF1}{A numeric value for assign the threshold for heterozygosity in mother.The default value is 5.}
  \item{HetM1}{A numeric value for assign the threshold for heterozygosity in father.The default value is 5.}
}

\details{
}

\value{
The function returns a dataframe final1
\describe{
\item{$Line}{Name of designation.}
\item{$Female}{Name of mother.}
\item{$Male}{Name of father.}
\item{$wrong}{Number of wrong triplets found it.}
\item{$QCPer_wrong}{Number of	testable triplets found it, percent of wrong triplets found it.}
\item{$QCPer_wrong_class}{Classification assigned depend of the percent of wrong triplets found it.}
\item{$LineImputed}{number of missing values in designation information.}
\item{$Line_imputed_class}{The classification assigned depend of number of missing values in designation information.} \item{$G_matrix_Line_Female}{G corresponding value to the pair of designation and mother.}
\item{$G_matrix_Line_Female_Class}{The classification assigned to the mother depend of G value obtained in G_matrix_Line_Female.}
\item{$G_matrix_Line_Male}{G corresponding value to the pair of designation and father.}
\item{$G_matrix_Line_Male_Class}{The classification assigned to the father depend of G value obtained in G_matrix_Line_Male.}
\item{$Hets_Line}{Heterozygosity calculation value for designation.}
\item{$Hets_Line_class}{The classification assigned depend of value of heterozygosity for designation}
\item{$Hets_Female}{Heterozygosity calculation value for mother.}
\item{$Hets_Female_class}{The classification assigned depend of value of heterozygosity for mother}
\item{$Hets_Male}{Heterozygosity calculation value for father.}
\item{$Hets_Male_class}{The classification assigned depend of value of heterozygosity for father}
\item{$Diagonal_minus_Female_G_matrix}{G corresponding value to the designation minus G corresponding value to the pair of designation and mother.}
\item{$Diagonal_Female_G_matrix_class}{The classification assigned depend of G value obtained in Diagonal_minus_Female_G_matrix.}
\item{$Imputed_Female}{Number of missing values in mother information.}
\item{$Imputed_Female_class}{The classification assigned depend of number of missing values in mother information.}
\item{$Imputed_Male}{Number of missing values in father information.}
\item{$Imputed_Male_class}{The classification assigned depend of number of missing values in father information.}
\item{$Final_Classification}{Considering all the previous classifications, a final classification is given.}
\item{$KEEP_0}{O is assigned if the final classification is positive and 1 otherwise.}
}
}

\references{

}

\examples{

}
