selfEvidenceStatus <- function(status, selfEvidence, nDiffCross = NA){

  # Final F1 QA/QC status from the probMatch status and the self evidence of crossVerification().
  # The evidence can only make a status stricter:
  #   self-supported (1, 2, 3) -> BAD;  unresolved (0) -> GOOD becomes QUESTIONABLE;  F1-supported (-1) -> unchanged.
  # PARENT FAIL / NO GENO DATA, and samples without evidence (NA), are left unchanged.
  # Returns data.frame(status, Warning); Warning is "" when there is nothing to report.
  n <- length(status)
  selfEvidence <- rep_len(as.numeric(selfEvidence), n)
  nDiffCross <- rep_len(as.numeric(nDiffCross), n)
  act <- status %in% c("GOOD", "QUESTIONABLE", "BAD") & !is.na(selfEvidence)

  new <- status
  new[act & selfEvidence %in% 1:3] <- "BAD"
  new[act & selfEvidence %in% 0 & status %in% "GOOD"] <- "QUESTIONABLE"

  w1 <- rep("", n)
  w1[act & new %in% "QUESTIONABLE" & selfEvidence %in% 0] <- "Cannot distinguish F1 from self"
  w1[act & new %in% "BAD" & selfEvidence %in% 1] <- "Matches mother (likely self)"
  w1[act & new %in% "BAD" & selfEvidence %in% 2] <- "Matches father (likely self)"
  w1[act & new %in% "BAD" & selfEvidence %in% 3] <- "Not F1 (mixed parental calls)"
  w2 <- rep("", n)
  few <- act & nDiffCross %in% 1:2
  w2[few] <- sprintf("Parents differ at only %d marker%s", nDiffCross[few], ifelse(nDiffCross[few] == 1, "", "s"))

  Warning <- ifelse(w1 != "" & w2 != "", paste(w1, w2, sep = "; "), paste0(w1, w2))
  data.frame(status = new, Warning = Warning, stringsAsFactors = FALSE)
}
