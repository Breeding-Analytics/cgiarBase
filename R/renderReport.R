renderReport <- function(object){
  src <- normalizePath(system.file("rmd",paste0("report",object,".Rmd"),package="bioflow"))
  src2 <- normalizePath(paste0('data/result',object,'.RData'))

  # temporarily switch to the temp dir, in case you do not have write
  # permission to the current working directory
  owd <- setwd(tempdir())
  on.exit(setwd(owd))

  file.copy(src, 'report.Rmd', overwrite = TRUE)
  file.copy(src2, paste0('result',object,'.RData'), overwrite = TRUE)

  outReport <- rmarkdown::render('report.Rmd', params = list(toDownload=TRUE),
                                 switch("HTML",HTML = rmdformats::robobook(toc_depth = 4)))

  return(outReport)
}
