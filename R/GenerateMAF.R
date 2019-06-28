#' Generate MAF file out of QTL data (from csv)
#'
#' @return list; Influenced Genes
#' @export
#'
#' @examples
GenerateMAF <- function(){

  # Read given files
  cis_eQTL_file <- system.file("extdata", "BRCA_tumor.cis_eQTL.csv", package = "omicsimulator2.0")
  trans_eQTL_file <- system.file("extdata", "BRCA_tumor.trans_eQTL.csv", package = "omicsimulator2.0")

  cis_eQTLs <- read.table(cis_eQTL_file, sep=";", header = TRUE)
  trans_eQTLs <- read.table(trans_eQTL_file, sep=";", header = TRUE)

  # Combine eQTLs
  names(cis_eQTL) <- c("eQTLs", "SNP_chr", "SNP_position", "SNP_alleles", "egenes", "gene_position")
  names(trans_eQTL) <- c("eQTLs", "SNP_chr", "SNP_position", "SNP_alleles", "egenes", "gene_position")
  eQTL <- rbind(cis_eQTL, trans_eQTL)

  # Select random 70% of eQTLs
  eQTL_selection <- sample(nrow(eQTL), 0.7 * nrow(eQTL))
  eQTL <- eQTL[eQTL_selection, ]

  ######TEST###########
  eQTL <-  eQTL[1:10, ]
  ######TEST###########

  # Transfer eQTLs to maf format
  chrom <- eQTL[,2]
  start = eQTL[,3]
  end = eQTL[,3]
  ref = sub('\\/.', '', eQTL[,4])
  alt = sub('.*\\/', '', eQTL[,4])

  sample.var = data.frame(chromsome = chrom, start = start, end = end, ref = ref, alt = alt)
  write.table(sample.var, 'sampleVars.txt', sep='\t',quote = FALSE, row.names = FALSE)
  var.maf <- maftools::oncotate(maflite = 'sampleVars.txt', header = TRUE)

  print(var.maf[, 1:16])

  # Return influenced genes
  return(eQTL[5])

}
