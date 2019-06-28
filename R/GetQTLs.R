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
  # Minimum: Chromosome	Start_Position	Reference_Allele	Tumor_Seq_Allele2	Tumor_Sample_Barcode

  # Annotate avinput format R data and file using ANNOVAR
  annovar.dir <- "/opt/bin/annovar"
  database.dir <- "/opt/bin/annovar/humandb"
  chr = "chr1"
  start = "123"
  end = "123"
  ref = "A"
  alt = "C"
  dat <- data.table(chr, start, end, ref, alt)
  tmpfn <- tempfile()
  write.table(dat, fn, row.names = FALSE, quote = FALSE, sep = "\t", col.names = FALSE)
  x <- annotation(dat, "perl_annovar_refGene", annovar.dir = "/opt/bin/annovar",
                  database.dir = database.dir)
  x <- annotation(input.file = tmpfn, "perl_annovar_refGene", annovar.dir = "/opt/bin/annovar",
                  database.dir = database.dir)


  maf <- maftools::maf2maf(cis_eQTL)

  print(maf)

  sample.var = data.frame(chr = c('chr4', 'chr15'), start = c(55589774, 41961117),
                          end = c(55589774, 41961117), ref_allele = c('A', 'TGGCTAA'), alt_allele = c('G', '-'))
  write.table(sample.var, file='sampleVars.tsv', quote=FALSE, sep='\t', row.names = FALSE)
  #write.table(sample.var, 'sampleVars.txt', sep='\t',quote = FALSE, row.names = FALSE)
  var.maf <- oncotate(maflite = 'sampleVars.tsv', header = TRUE)

  var.maf <- icgcSimpleMutationToMAF(icgc = 'sampleVars.tsv', addHugoSymbol = TRUE)

  var.file = system.file("extdata", "sampleVars.tsv", package = "omicsimulator2.0")
  print(var.file)
  var.maf = oncotate(maflite = var.file, header = TRUE)


  sample.var = data.frame(chromsome = c('chr4', 'chr15'), Start = c(55589774, 41961117),
                          end = c(55589774, 41961117), ref = c('A', 'TGGCTAA'), alt = c('G', '-'),
                          Tumor_Sample_Barcode = c('fake_1', 'fake2'))
  write.table(sample.var, 'sampleVars.txt', sep='\t',quote = FALSE, row.names = FALSE)
  var.maf <- oncotate(maflite = 'sampleVars.txt', header = TRUE)

  print(cis_eQTL)

}
