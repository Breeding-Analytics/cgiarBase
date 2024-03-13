\name{printTable}
\alias{printTable}

\title{
Printing a table 
}
\description{
Basic function to print a table.
}
\usage{
printTable(DT, pageLength = 5, 
          numericColNames = NULL, numRound = 3, 
          scrollXOpt = TRUE, ...)
}

\arguments{
  \item{DT}{
Data frame to be printed as table.
}
\item{pageLength}{
The number of items that should be presented by page.
}
\item{numericColNames}{
Character vector specifying which columns should be considered numeric.
}
\item{numRound}{
The number of decimals to round up to.
}
\item{scrollXOpt}{
Should the scrolling option should be activated.
}
\item{...}{
Additional parameters to be passed to the data.table function.
}

}
\details{

Simple function to print a table.

}
\value{

An interactive table. 

}
}

\references{
None
}

\examples{

data(DT_example)

dt <- result$data$pedigree

printTable(DT=dt, pageLength = 5, 
           numericColNames = NULL, numRound = 3, 
           scrollXOpt = TRUE)

}