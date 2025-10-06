replaceValues <- function(Source, Search, Replace){
  if (length(Search) != length(Replace))
    stop("Search and Replace Must Have Equal Number of Items\n")
  Changed <- as.character(Source)

  for (i in 1:length(Search)) {
    if(!is.na(Replace[i])){
      if((Replace[i] %in% Changed) & (Replace[i] != Search[i])){
        Changed <- replace(Changed, Changed == Replace[i], paste0(Replace[i],"_dup"))
      }
    }
  }

  for (j in 1:length(Search)) {
    Changed <- replace(Changed, Changed == Search[j], Replace[j])
  }
  # silly comment

  if(is.numeric(Replace))
    Changed <- as.numeric(Changed)

  return(Changed)
}
