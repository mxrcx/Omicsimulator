#' Generate MAF file out of QTL data (from csv)
#'
#' @param threshold_eQTls Int; Threshold of the number of eQTLs used
#' @param tumor_sample_barcodes
#'
#' @return list; Influenced Genes
#' @export
#'
#' @examples
GenerateMAF <- function(threshold_eQTls, tumor_sample_barcodes, output_directory, file_prefix, disease){

  # Calculate all standard derivations in advance

  sd_list <- sd(tcga_matrix_normal[1, 1:ncol(tcga_matrix_normal)])

  for (row in 2:nrow(tcga_matrix_normal)){

    # Calculate Standard Derivation for gene (from normal matrix)
    sd_list <- c(sd_list, sd(tcga_matrix_normal[row, 1:ncol(tcga_matrix_normal)]))

  }

  # Create maf
  maf_file <- NULL

  # Read given files
  cis_eQTL_file <- system.file("extdata", "BRCA_tumor.cis_eQTL.csv", package = "omicsimulator2.0")
  cis_eQTLs <- read.table(cis_eQTL_file, sep=";", header = TRUE)
  names(cis_eQTL) <- c("eQTLs", "SNP_chr", "SNP_position", "SNP_alleles", "egenes", "gene_position")

  trans_eQTL_file <- system.file("extdata", "BRCA_tumor.trans_eQTL.csv", package = "omicsimulator2.0")
  trans_eQTLs <- read.table(trans_eQTL_file, sep=";", header = TRUE)
  names(trans_eQTL) <- c("eQTLs", "SNP_chr", "SNP_position", "SNP_alleles", "egenes", "gene_position")

  # Combine eQTLs
  eQTL <- rbind(cis_eQTL, trans_eQTL)

  # For each sample select x % of eQTLS which have impact (depending on threshold value)
  for (sample in tumor_sample_barcodes){

    cat(sample, "\n")

    # Randomly select eQTLs with impact
    eQTL_selection <- sample(nrow(eQTL), threshold_eQTls * nrow(eQTL), replace = F)
    eQTL_with_impact <- eQTL[eQTL_selection, ]

    # Randomly selct eQTLS with no impact
    eQTL_no_impact <- eQTL[-(eQTL_selection), ]
    eQTL_no_impact_selection <- sample(nrow(eQTL_no_impact), 0.3 * nrow(eQTL_no_impact), replace = F)
    eQTL_no_impact <- eQTL_no_impact[eQTL_no_impact_selection, ]

    # Combine eQTLs with and without impact
    eQTL_current <- rbind(eQTL_with_impact, eQTL_no_impact)

    ####TEST############
    #eQTL_current <- eQTL_current[1:3, ]
    ####################

    # Transfer eQTLs to maf format
    chrom <- eQTL_current[,2]
    start = eQTL_current[,3]
    end = eQTL_current[,3]
    ref = sub('\\/.', '', eQTL_current[,4])
    alt = sub('.*\\/', '', eQTL_current[,4])
    genes = eQTL_current[,5]

    # Randomly select SOMATIC
    somatic <- sample(0:1, length(chrom), replace = TRUE)

    # Randomly select impact values
    sift <- NULL
    polyphen <- NULL
    impact <- NULL

    # Create progress bar for random value generation
    progress_bar_random_values <- txtProgressBar(min = 0, max = nrow(eQTL_current), style = 3)

    for (eQTL_entry in 1:nrow(eQTL_current)){

      setTxtProgressBar(progress_bar_random_values, eQTL_entry)

      repeat{

        if(eQTL_entry <= length(eQTL_selection)){         # with impact

          random_sift <- sample(1:4, 1)
          random_polyphen <- sample(1:4, 1)
          random_impact <- sample(1:4, 1)

          # At least one field has an impact
          if((random_sift + random_polyphen + random_impact) < 9){
            break
          }
        }
        else {                                            # no impact

          random_sift <- sample(3:4, 1)
          random_polyphen <- sample(3:4, 1)
          random_impact <- sample(3:4, 1)

          # All fields no impact
          if((random_sift > 2) && (random_polyphen > 2) && (random_impact > 2)){
            break
          }
        }
      }

      ######### SIFT #########
      if (random_sift == 1){
        sift <- c(sift, "deterious(0.3)")
      }
      if (random_sift == 2){
        sift <- c(sift, "deterious_low_confidence(0.3)")
      }
      if (random_sift == 3){
        sift <- c(sift, "tolerated_low_confidence(0.3)")
      }
      if (random_sift == 4){
        sift <- c(sift, "tolerated(0.3)")
      }

      ######### PolyPhen #########
      if (random_polyphen == 1){
        polyphen <- c(polyphen, "probably damaging(0.3)")
      }
      if (random_polyphen == 2){
        polyphen <- c(polyphen, "possibly damaging(0.3)")
      }
      if (random_polyphen == 3){
        polyphen <- c(polyphen, "benign(0.3)")
      }
      if (random_polyphen == 4){
        polyphen <- c(polyphen, "unknown(0.3)")
      }

      ######### IMPACT #########
      if (random_impact == 1){
        impact <- c(impact, "HIGH(0.3)")
      }
      if (random_impact == 2){
        impact <- c(impact, "MODERATE(0.3)")
      }
      if (random_impact == 3){
        impact <- c(impact, "LOW(0.3)")
      }
      if (random_impact == 4){
        impact <- c(impact, "MODIFIER(0.3)")
      }
    }

    close(progress_bar_random_values)


    # Build sample_maf

    sample_maf = data.frame(Chromsome = chrom, Start_Position = start, End_Position = end,
                            Reference_Allele = ref, Tumor_Seq_Allele1 = ref, Tumor_Seq_Allele2 = alt,
                            Tumor_Sample_Barcode = sample, Gene = genes, SIFT = sift, PolyPhen = polyphen,
                            SOMATIC = somatic, IMPACT = impact)

    maf_file <- rbind(maf_file, sample_maf)


    # Simulate influenced Genes per sample

    influenced_genes <- sample_maf[,8]

    cat("Size: ", length(influenced_genes))

    influenced_genes <- unique(influenced_genes)

    cat("Size: ", length(influenced_genes))

    genes_dictionary_from_eQTL <- hash::hash(influenced_genes, rep(1, length(influenced_genes)))


    # Simulate gene expression values of influenced genes

    simulated_matrix <- tcga_matrix_normal

    # Create progress bar
    progress_bar <- txtProgressBar(min = 0, max = nrow(simulated_matrix), style = 3)

    for (row in 1:nrow(simulated_matrix)){

      setTxtProgressBar(progress_bar, row)

      # Manipulate expression values
      if(rownames(simulated_matrix)[row] %in% ls(genes_dictionary_from_eQTL)){

        # Get Standard Derivation for gene
        sd <- sd_list[row]

        # Increasing or decreasing Influence
        influence <- get(toString(rownames(simulated_matrix)[row]), genes_dictionary_from_eQTL)

        # Increase or decrease expression value
        simulated_matrix[row, which(tumor_sample_barcodes == sample)] <- simulated_matrix[row, which(tumor_sample_barcodes == sample)] + sd * influence * 20

      }

    }

    close(progress_bar)

  }

  write.table(maf_file, file = file.path(output_directory, disease, paste(file_prefix, "_MAF_file.maf", sep="")), quote=FALSE, sep ="\t", row.names = TRUE)

  print(maf_file)

  cat("MAF-File Einträge: ", nrow(maf_file))

  # Return influenced genes
  return(maf_file)

}
