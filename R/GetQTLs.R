#' Get QTLs out of given csv files
#'
#' @return
#' @export
#'
#' @examples
GetQTLs <- function(){

  # Read given files
  cis_eQTL_file <- system.file("extdata", "BRCA_tumor.cis_eQTL.csv", package = "omicsimulator2.0")
  trans_eQTL_file <- system.file("extdata", "BRCA_tumor.trans_eQTL.csv", package = "omicsimulator2.0")

  cis_eQTL <- read.table(cis_eQTL_file, header = TRUE)
  trans_eQTL <- read.table(trans_eQTL_file, header = TRUE)

  # Select random 70% of eQTLs
  cis_eQTL_selection <- sample(nrow(cis_eQTL), 0.7 * nrow(cis_eQTL))
  trans_eQTL_selection <- sample(nrow(trans_eQTL), 0.7 * nrow(trans_eQTL))

  cis_eQTL <- cis_eQTL[cis_eQTL_selection, ]
  trans_eQTL <- trans_eQTL[trans_eQTL_selection, ]

  # Save eQTL in maf format


  print(cis_eQTL)

}
