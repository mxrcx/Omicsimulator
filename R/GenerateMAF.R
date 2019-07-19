#' Generate MAF file out of QTL data (from csv)
#'
#' @param threshold_eQTls Int; Threshold of the number of eQTLs used
#'
#' @return list; Influenced Genes
#' @export
#'
#' @examples
GenerateMAF <- function(threshold_eQTls){

  # Read given files
  cis_eQTL_file <- system.file("extdata", "BRCA_tumor.cis_eQTL.csv", package = "omicsimulator2.0")
  trans_eQTL_file <- system.file("extdata", "BRCA_tumor.trans_eQTL.csv", package = "omicsimulator2.0")

  cis_eQTLs <- read.table(cis_eQTL_file, sep=";", header = TRUE)
  trans_eQTLs <- read.table(trans_eQTL_file, sep=";", header = TRUE)

  # Combine eQTLs
  names(cis_eQTL) <- c("eQTLs", "SNP_chr", "SNP_position", "SNP_alleles", "egenes", "gene_position")
  names(trans_eQTL) <- c("eQTLs", "SNP_chr", "SNP_position", "SNP_alleles", "egenes", "gene_position")
  eQTL <- rbind(cis_eQTL, trans_eQTL)

  # Select random x % of eQTLs (depending on threshold value)
  eQTL_selection <- sample(nrow(eQTL), threshold_eQTls * nrow(eQTL), replace = F)

  eQTL_no_impact <- eQTL[-(eQTL_selection), ] # eQTLS with no impact

  eQTL <- eQTL[eQTL_selection, ]

  # 30% of eQTLS with no impact
  eQTL_no_impact_selection <- sample(nrow(eQTL_no_impact), 0.3 * nrow(eQTL_no_impact), replace = F)
  eQTL_no_impact <- eQTL_no_impact[eQTL_no_impact_selection, ]

  # Combine eQTLs with and without impact
  eQTL <- rbind(eQTL, eQTL_no_impact)

  ####TEST############
  eQTL <- eQTL[1, ]
  ####################

  # Transfer eQTLs to maf format
  chrom <- eQTL[,2]
  start = eQTL[,3]
  end = eQTL[,3]
  ref = sub('\\/.', '', eQTL[,4])
  alt = sub('.*\\/', '', eQTL[,4])

  sample.var = data.frame(chromsome = chrom, start = start, end = end, ref = ref, alt = alt)
  write.table(sample.var, 'sampleVars.txt', sep='\t',quote = FALSE, row.names = FALSE)
  var.maf <- maftools::oncotate(maflite = 'sampleVars.txt', header = TRUE)

  # Manipulate the impact field
  for (eQTL_entry in 1:length(eQTL_selection)){
    impact <- sample(1:18, 1)

    if (impact == (1 || 7 || 8 || 9 || 10)){
      # SIFT --> deleterious(0.3)
      var.maf["SIFT", eQTL_entry] <- "deterious(0.3)"
    }
    if (impact == (2 || 11 || 12 || 13 || 14)){
      # SIFT --> deleterious_low_confidence(0.3)
      var.maf["SIFT", eQTL_entry] <- "deterious_low_confidence(0.3)"
    }
    if (impact == (3 || 7 || 11 || 15 || 16)){
      # PolyPhen --> probably damaging (PR)
      var.maf["PolyPhen", eQTL_entry] <- "PR(0.3)"
    }
    if (impact == (4 || 8 || 12 || 17 || 18)){
      # PolyPhen --> possibly damaging (PO)
      var.maf["PolyPhen", eQTL_entry] <- "PO(0.3)"
    }
    if (impact == (5 || 9 || 13 || 15 || 17)){
      # IMPACT (VEP) --> HIGH (H)
      var.maf["IMPACT", eQTL_entry] <- "H(0.3)"
    }
    if (impact == (6|| 10 || 14 || 16 || 18)){
      # IMPACT (VEP) --> MODERATE (M)
      var.maf["IMPACT", eQTL_entry] <- "M(0.3)"
    }
  }


  # Manipulate impact field from eQTLs with no impact
  for(eQTL_entry in (length(eQTL_selection)+1):nrow(var.maf)){
    impact <- sample(1:18, 1)

    if (impact == (1 || 7 || 8 || 9 || 10)){
      # SIFT --> tolerated(0.3)
      var.maf["SIFT", eQTL_entry] <- "tolerated(0.3)"
    }
    if (impact == (2 || 11 || 12 || 13 || 14)){
      # SIFT --> tolerated_low_confidence(0.3)
      var.maf["SIFT", eQTL_entry] <- "tolerated_low_confidence(0.3)"
    }
    if (impact == (3 || 7 || 11 || 15 || 16)){
      # PolyPhen --> unknown (UN)
      var.maf["PolyPhen", eQTL_entry] <- "UN"
    }
    if (impact == (4 || 8 || 12 || 17 || 18)){
      # PolyPhen --> benign (BE)
      var.maf["PolyPhen", eQTL_entry] <- "BE(0.3)"
    }
    if (impact == (5 || 9 || 13 || 15 || 17)){
      # IMPACT (VEP) --> MODIFIER (MO)
      var.maf["IMPACT", eQTL_entry] <- "MO(0.3)"
    }
    if (impact == (6|| 10 || 14 || 16 || 18)){
      # IMPACT (VEP) --> Low (L)
      var.maf["IMPACT", eQTL_entry] <- "L(0.3)"
    }
  }

  # Return influenced genes
  return(c(var.maf, as.vector(eQTL[,5])))

}
