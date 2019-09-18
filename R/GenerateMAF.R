#' Generate MAF file out of QTL data (from csv)
#'
#' @param threshold_eQTls Int; Threshold of the number of eQTLs used
#' @param tumor_sample_barcodes
#'
#' @return Matrix; Simulated Matrix
#' @export
#'
#' @examples
GenerateMAF <- function(threshold_eQTls, tumor_sample_barcodes, output_directory, file_prefix, disease, tcga_matrix_normal, tcga_matrix_tumor){

  simulated_matrix <- tcga_matrix_normal

  cat("Calculate all standard derivations in advance...\n")
  # Calculate all standard derivations in advance

  sd_list_tumor <- sd(tcga_matrix_tumor[1, 1:ncol(tcga_matrix_tumor)])
  for (row in 2:nrow(tcga_matrix_tumor)){
    # Calculate Standard Derivation for gene (from tumor matrix)
    sd_list_tumor <- c(sd_list_tumor, sd(tcga_matrix_tumor[row, 1:ncol(tcga_matrix_tumor)]))
  }

  mean_list_tumor <- mean(tcga_matrix_tumor[1, 1:ncol(tcga_matrix_tumor)])
  for (row in 2:nrow(tcga_matrix_tumor)){
    # Calculate Mean for gene (from tunor matrix)
    mean_list_tumor <- c(mean_list_tumor, mean(tcga_matrix_tumor[row, 1:ncol(tcga_matrix_tumor)]))
  }

  sd_list_normal <- sd(tcga_matrix_normal[1, 1:ncol(tcga_matrix_normal)])
  for (row in 2:nrow(tcga_matrix_normal)){
    # Calculate Standard Derivation for gene (from normal matrix)
    sd_list_normal <- c(sd_list_normal, sd(tcga_matrix_normal[row, 1:ncol(tcga_matrix_normal)]))
  }

  mean_list_normal <- mean(tcga_matrix_normal[1, 1:ncol(tcga_matrix_normal)])
  for (row in 2:nrow(tcga_matrix_normal)){
    # Calculate Mean for gene (from normal matrix)
    mean_list_normal <- c(mean_list_normal, mean(tcga_matrix_normal[row, 1:ncol(tcga_matrix_normal)]))
  }

  cat("Prepare MAF-Files...\n")

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

    cat(sample, ":\n", "Generate random values...\n")

    # Randomly select eQTLs with impact
    eQTL_selection <- sample(nrow(eQTL), threshold_eQTls * nrow(eQTL), replace = F)
    eQTL_with_impact <- eQTL[eQTL_selection, ]

    # Randomly selct eQTLS with no impact
    eQTL_no_impact <- eQTL[-(eQTL_selection), ]
    no_impact_threshold <- (0.3 * threshold_eQTls) / 0.7
    eQTL_no_impact_selection <- sample(nrow(eQTL_no_impact), no_impact_threshold * nrow(eQTL_no_impact), replace = F)
    eQTL_no_impact <- eQTL_no_impact[eQTL_no_impact_selection, ]

    # Combine eQTLs with and without impact
    eQTL_current <- rbind(eQTL_with_impact, eQTL_no_impact)

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
    sift <- vector(length = nrow(eQTL_current))
    polyphen <- vector(length = nrow(eQTL_current))
    impact <- vector(length = nrow(eQTL_current))

    # Random value generation
    for(eQTL_entry in 1:nrow(eQTL_current)){

      if(eQTL_entry <= length(eQTL_selection)){         # with impact
        # At least one of the three fields has an impact (Level 1 or 2)
        random_value_influence <- sample(1:2, 1, replace = TRUE)
        random_values_maybe_influence <- sample(1:4, 2, replace = TRUE)

        #Randomize vector
        random_values <- sample(c(random_value_influence, random_values_maybe_influence))
      }
      else {                                            # no impact
        # All fields no impact (Level 3 or 4)
        random_values <- sample(3:4, 3, replace = TRUE)
      }

      sift[eQTL_entry] <- random_values[1]
      polyphen[eQTL_entry] <- random_values[2]
      impact[eQTL_entry] <- random_values[3]

    }

    # SIFT
    sift <- replace(sift, sift == 1, "deterious(0.3)")
    sift <- replace(sift, sift == 2, "deterious_low_confidence(0.3)")
    sift <- replace(sift, sift == 3, "tolerated_low_confidence(0.3)")
    sift <- replace(sift, sift == 4, "tolerated(0.3)")
    # PolyPhen
    polyphen <- replace(polyphen, polyphen == 1, "probably damaging(0.3)")
    polyphen <- replace(polyphen, polyphen == 2, "possibly damaging(0.3)")
    polyphen <- replace(polyphen, polyphen == 3, "benign(0.3)")
    polyphen <- replace(polyphen, polyphen == 4, "unknown(0.3)")
    # IMPACT
    impact <- replace(impact, impact == 1, "HIGH(0.3)")
    impact <- replace(impact, impact == 2, "MODERATE(0.3)")
    impact <- replace(impact, impact == 3, "LOW(0.3)")
    impact <- replace(impact, impact == 4, "MODIFIER(0.3)")

    # Build sample_maf
    sample_maf = data.frame(Chromsome = chrom, Start_Position = start, End_Position = end,
                            Reference_Allele = ref, Tumor_Seq_Allele1 = ref, Tumor_Seq_Allele2 = alt,
                            Tumor_Sample_Barcode = sample, Gene = genes, SIFT = sift, PolyPhen = polyphen,
                            SOMATIC = somatic, IMPACT = impact)

    maf_file <- rbind(maf_file, sample_maf)


    # Simulate influenced Genes per sample
    influenced_genes <- eQTL_with_impact[,5]
    influenced_genes <- unique(influenced_genes)
    genes_dictionary_from_eQTL <- hash::hash(influenced_genes, rep(1, length(influenced_genes)))


    # Simulate gene expression values of influenced genes

    # Create progress bar
    #progress_bar <- txtProgressBar(min = 0, max = nrow(simulated_matrix), style = 3)

    ####### NEWWWWWW
    wanted_rows = which(rownames(simulated_matrix) %in% ls(genes_dictionary_from_eQTL))

    print(system.time({
      simulated_matrix[wanted_rows, which(tumor_sample_barcodes == sample)] <- ((simulated_matrix[wanted_rows, which(tumor_sample_barcodes == sample)] - mean_list_normal[wanted_rows])/sd_list_normal[wanted_rows]*sd_list_tumor[wanted_rows])+mean_list_tumor[wanted_rows]
    }))

    #print(simulated_matrix)

    #simulated_matrix <- tcga_matrix_normal


    simulate_gene_expression <- function (row){

      #setTxtProgressBar(progress_bar, row)

      # Manipulate expression values
      if(rownames(simulated_matrix)[row] %in% ls(genes_dictionary_from_eQTL)){

        # Get Standard Derivation for genes
        sd_normal <- sd_list_normal[row]
        sd_tumor <- sd_list_tumor[row]

        # Get Means
        mean_normal <- mean_list_normal[row]
        mean_tumor <- mean_list_tumor[row]

        # Increasing or decreasing Influence
        influence <- get(toString(rownames(simulated_matrix)[row]), genes_dictionary_from_eQTL)

        # Increase or decrease expression value
        simulated_matrix[row, which(tumor_sample_barcodes == sample)] <- ((simulated_matrix[row, which(tumor_sample_barcodes == sample)] - mean_normal)/sd_normal*sd_tumor)+mean_tumor

      }

    }

    #library(parallel)

    #print(system.time({
     #results <- parallel::mclapply(1:nrow(simulated_matrix), simulate_gene_expression, mc.cores = detectCores())
    #}))

    #close(progress_bar)

  }

  write.table(maf_file, file = file.path(output_directory, disease, paste(file_prefix, "_MAF_file.maf", sep="")), quote=FALSE, sep ="\t", row.names = TRUE)

  cat("MAF-File Einträge: ", nrow(maf_file), "\n")

  # Return simulated matrix
  return(simulated_matrix)

}
