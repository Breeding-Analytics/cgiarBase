# Function to give format to a data.frame
roundDF <- function(dataFrame, numRound) {
  data.frame(lapply(X = dataFrame,
                    FUN = function(x) {
                      if(is.numeric(x)){
                        n1 <- lapply(x, function(x) nchar(strsplit(as.character(x), "\\.")[[1]][1]))
                        n1 <- min(unlist(n1), na.rm = TRUE)
                        if (n1 >= numRound) {
                          round(x, 0)
                        } else {
                          if (min(x, na.rm = TRUE) >= 1) {
                            round(x, numRound - n1)
                          } else {
                            tmp <- signif(x, numRound)
                            n2 <- lapply(tmp, function(x) nchar(strsplit(as.character(x), "\\.")[[1]][2]))
                            n2 <- max(unlist(n2), na.rm = TRUE)
                            round(x, n2)
                          }
                        }
                      } else {
                        x
                      }
                    }))
}

# Function to give format to a datatable
roundDT <- function(dataTable, numRound) {
  x <- dataTable$x$data
  numericColNames <- unlist(lapply(x, is.numeric))
  for (i in 1:length(numericColNames)) {
    if (numericColNames[i]) {
      n1 <- lapply(x[, i], function(x) nchar(strsplit(as.character(x), "\\.")[[1]][1]))
      n1 <- unlist(n1)
      for (j in 1:dim(x)[1]) {
        if (n1[j] >= numRound)
          dataTable$x$data[j, names(numericColNames[i])] <- round(x[j, names(numericColNames[i])], 0)
        if (n1[j] < numRound)
          dataTable$x$data[j, names(numericColNames[i])] <- signif(x[j, names(numericColNames[i])], numRound)
      }
    }
  }
  dataTable
}
