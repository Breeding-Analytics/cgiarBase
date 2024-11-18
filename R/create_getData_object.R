create_getData_object <- function() {
  obj <- list()

  obj$data <- list()
  obj$metadata <- list()
  obj$modifications <- list()

  obj$data$pheno    <- data.frame()
  obj$data$geno     <- matrix()
  obj$data$weather  <- data.frame()
  obj$data$pedigree <- data.frame(designation = character(),
                                  mother = character(),
                                  father = character(),
                                  yearOfOrigin = integer())
  obj$data$qtl     <- matrix()

  obj$metadata$pheno <- data.frame(parameter = character(),
                                   value = character())

  obj$metadata$geno <- data.frame(marker = character(),
                                  chr = character(),
                                  pos = double(),
                                  refAllele = character(),
                                  altAllele = character())

  obj$metadata$pedigree <- data.frame(parameter = character(),
                                      value = double())

  obj$metadata$weather <- data.frame(parameter = character(),
                                     value = double())

  obj$metadata$qtl <- data.frame(parameter = character(),
                                 value = double())

  obj$modifications$pheno <- data.frame(module = character(),
                                        analysisId = character(),
                                        trait = character(),
                                        reason = character(),
                                        row = integer(),
                                        value = double())

  obj$modifications$geno <- data.frame(module = character(),
                                       analysisId = character(),
                                       reason = character(),
                                       row = integer(),
                                       col = integer(),
                                       value = double())

  obj$modifications$pedigree <- data.frame(module = character(),
                                           analysisId = character(),
                                           reason = character(),
                                           row = integer())

  obj$modifications$weather <- data.frame()

  obj$modifications$qtl <- data.frame()

  obj$predictions <- data.frame(module = character(),
                                analysisId = character(),
                                pipeline = character(),
                                trait = character(),
                                gid = character(),
                                designation = character(),
                                mother = character(),
                                father = character(),
                                effectType = character(),
                                entryType = character(),
                                environment = character(),
                                predictedValue = double(),
                                stdError = double(),
                                reliability = double())

  obj$metrics <- data.frame(module = character(),
                            analysisId = character(),
                            trait = character(),
                            environment = character(),
                            parameter = character(),
                            method = character(),
                            value = double(),
                            stdError = double())

  obj$modeling <- data.frame(module = character(),
                             analysisId = character(),
                             trait = character(),
                             environment = character(),
                             parameter = character(),
                             value = character())

  obj$status <- data.frame(module = character(),
                           analysisId = character(),
                           analysisIdName=character())

  return(obj)
}
