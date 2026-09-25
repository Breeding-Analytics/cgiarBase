## Resolve a genotype QA/QC stamp to position(s) in a geno_imp list.
##
## The genotype QA/QC module records the same analysisId in two different forms:
##
##   status$analysisId          as.numeric(Sys.time())            full-precision double
##   names(data$geno_imp)       as.character(round(analysisId))   rounded to whole seconds
##
## Every consumer builds its version dropdown from status$analysisId, so the value
## handed to the matching code is the unrounded one. Comparing it as a string against
## the rounded list key therefore never matches, and the consumer either errors on
## geno_imp[[integer(0)]] or degrades silently to an empty marker frame. Objects written
## by older app versions (and by the EBS import scripts) key on the unrounded value and
## do match exactly, so both generations are in circulation and both must resolve.
##
## Matching is done on set membership rather than equality so that dropdowns declared
## with multiple = TRUE work: comparing with `==` against a multi-element selection
## compares position-wise and silently returns the wrong answer.
##
## Returns an integer vector of positions in `genoImp`. Callers are expected to check
## its length: zero means the stamp did not resolve, more than one means the stamp is
## ambiguous (two QA runs collapsing onto the same rounded key). Both cases are left to
## the caller because UI code needs to skip quietly while pipeline code must stop.
resolveGenoStamp <- function(genoImp, stamp){

  nms <- names(genoImp)
  if(is.null(nms) || length(nms) == 0){return(integer(0))}

  # Drop the placeholder values the various dropdowns use for "nothing selected":
  # "" (unset), "0" (Population Structure / Core Set / OCS) and the literal
  # "No data available" (MTA-ASReml). None can be a real Unix-timestamp key.
  stampChar <- trimws(as.character(stamp))
  stampChar <- stampChar[!is.na(stampChar)]
  stampChar <- stampChar[!stampChar %in% c("", "0", "No data available")]
  if(length(stampChar) == 0){return(integer(0))}

  # Exact key match: objects keyed on the unrounded analysisId.
  idx <- which(nms %in% stampChar)

  # Rounded match: objects keyed on as.character(round(analysisId)).
  if(length(idx) == 0){
    nmsNum   <- suppressWarnings(as.numeric(nms))
    stampNum <- suppressWarnings(as.numeric(stampChar))
    stampNum <- stampNum[!is.na(stampNum)]
    if(length(stampNum) > 0){
      idx <- which(!is.na(nmsNum) & round(nmsNum) %in% round(stampNum))
    }
  }

  return(idx)
}
