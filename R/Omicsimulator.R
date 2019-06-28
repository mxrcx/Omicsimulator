#' Main function which calls all relevant to simulate the influence of variations on gene expression
#'
#' Main function which calls all relevant to simulate the influence of variations on gene expression
#'
#' @return List; List containing the DEA results of both normal-tumor and normal-simulated
#' @export
Omicsimulator <-
function(disease, sample_number, top_DEG_number, output_directory, file_prefix)
{ ... }
#'
#' @param disease String; name of a disease to get a KEGG pathway for
#' @param output_directory String; directory of the output files
#' @param top_DEG_number Integer; number of top differential expressed genes
#' @param sample_number Integer; the number of samples which are simulated
#' @param file_prefix String; name of the results output file
#'
#' @examples
#' Omicsimulator(disease = "Breast cancer", sample_number = 10, top_DEG_number = 100)
Omicsimulator <- function(disease, sample_number, top_DEG_number, output_directory, file_prefix){

  # Check for dependencies
  CheckDependencies()

  # Start TOTAL timer
  tictoc::tic("TOTAL")

  # Set default disease
  if(missing(disease)){
    disease <- "Breast cancer"
  }

  # Set output directory
  if(missing(output_directory)) {
    output_directory <- "../output"
  }
  if(!file.exists(output_directory)){
    dir.create(output_directory)
  }
  if(!file.exists(file.path(output_directory, disease))){
    dir.create(file.path(output_directory, disease))
  }

  # Set file name
  if(missing(file_prefix)){
    file_prefix <- paste("omicsimulatorResults", disease, sample_number, top_DEG_number, sep = "_")
  }

  # Create cache directory
  if(!file.exists("cache")){
    dir.create("cache")
  }
  if(!file.exists(file.path("cache", disease))){
    dir.create(file.path("cache", disease))
  }

  # Set sample number [in case none is given as parameter]
  if(missing(sample_number)) {
    sample_number <- 10
  }

  # Set top differential expressed genes number [in case none is given as parameter]
  if(missing(top_DEG_number)) {
    top_DEG_number <- 100
  }


  ##########################################################################################################


  cat("----------------------------- \nOMICSIMULATOR START. \n----------------------------- \n")

  # Create input file from pathway

  # Start GET INPUT timer
  tictoc::tic("GET INPUT")

  cat("Create input file from pathway...\n")

  pathway_id <- GetPathwayID(disease)
  GenesFromPathway(pathway_id = pathway_id)

  # Get input data

  cat("Read input file...\n")

  input_data <- GetInputData()
  genes_variation <- input_data$genes_variation
  genes_dictionary <- input_data$genes_dictionary

  # Stop GET INPUT timer
  tictoc::toc()


  ##########################################################################################################



  # Start LOAD MATRICES timer
  tictoc::tic("LOAD MATRICES")

  # Load TCGA Matrix of specified TUMOR data
  sample_type <- "Primary solid Tumor"
  tcga_matrix_tumor <- LoadTCGAMatrix(disease, sample_type, sample_number)

  # Load TCGA Matrix of NORMAL tissue data
  sample_type <- "Solid Tissue Normal"
  tcga_matrix_normal <- LoadTCGAMatrix(disease, sample_type, sample_number)

  # Stop LOAD MATRICES timer
  tictoc::toc()

  # Start SIMULATION + DEA timer
  tictoc::tic("SIMULATION + DEA")

  # Do differential expression analysis (NORMAL + TUMOR)
  DEA_normal_tumor <- DEA(disease, sample_number, output_directory, tcga_matrix_normal, tcga_matrix_tumor, "NT_PT", file_prefix)
  top_DEG_real <- TopDEG(DEA_normal_tumor, top_DEG_number)



  ##########################################################################################################


  # Generate MAF File & get influenced genes
  influenced_genes <- GenerateMAF()

  # Simulate tumor

  #efo_code <- "EFO_0000305"
  #top_genes <- GetTopGenesFromOpentargets(efo_code, top_DEG_number)
  #genes_dictionary_from_opentargets <- hash::hash(top_genes, rep(1, top_DEG_number))

  genes_dictionary_from_eQTL <- hash::hash(influenced_genes, rep(1, length(influenced_genes)))

  cat("Simulate Counts: ")
  simulated_matrix <- SimulateCounts(tcga_matrix_normal, genes_dictionary_from_eQTL)

  # Do differential expression analysis (NORMAL + SIMULATED)
  DEA_normal_simulated <- DEA(disease, sample_number, output_directory, tcga_matrix_normal, simulated_matrix, "NT_Simulated", file_prefix)
  top_DEG_simulated <- TopDEG(DEA_normal_simulated, top_DEG_number)

  # Stop SIMULATION + DEA timer
  tictoc::toc()


  ##########################################################################################################


  # Start CORRELATION MATRIX timer
  tictoc::tic("CORRELATION MATRIX")

  cat("Load Correlation Matrix:\n")

  # Build correlation matrix to find coexpressed genes

  if(!file.exists(file.path("cache", disease, paste("cor_normal", disease, "_sample=", sample_number, ".rds", sep="")))){

    cat("--> Create Correlation Matrix...\n")

    cor_normal <- propagate::bigcor(t(tcga_matrix_normal[1:40000, 1:sample_number]), size = 10000, fun = "cor")

    # Save as RDS file in cache
    saveRDS(cor_normal, file.path("cache", disease, paste("cor_normal", disease, "_sample=", sample_number, ".rds", sep="")))

    cat("--> DONE. \n")
  }
  else{

    cor_normal <- readRDS(file.path("cache", disease, paste("cor_normal", disease, "_sample=", sample_number, ".rds", sep="")))

    cat("Loaded from save file.\n")
  }

  cat("Check for coexpressed genes... \n")

  # Find row numbers of top_DEG_real
  top_DEG_real_row_numbers <- NULL

  for(row_number in 1:nrow(tcga_matrix_normal[1:40000, 1:sample_number])){

    if(is.element(rownames(tcga_matrix_normal)[row_number], top_DEG_real)){

      top_DEG_real_row_numbers <- c(top_DEG_real_row_numbers, row_number)

    }
  }

  cat("DONE. \n")

  # Check for coexpressed genes
  cat("Get coexpressed genes numbers... \n")
  coexpressed_genes_numbers <- NULL
  threshold <- 0.8


  # Get only relevant rows
  cor_coexpressed <- cor_normal[top_DEG_real_row_numbers, ]

  # Create progress bar
  progress_bar <- txtProgressBar(min = 0, max = nrow(cor_coexpressed), style = 3)

  for(row_number in 1:nrow(cor_coexpressed)){

    setTxtProgressBar(progress_bar, row_number)

    for(col_number in (row_number + 1):ncol(cor_coexpressed)){

      cor_value <- cor_coexpressed[row_number, col_number]

      if((!is.na(cor_value)) && (cor_value >= 0.8)){

        coexpressed_genes_numbers <- c(coexpressed_genes_numbers, col_number)

      }

    }

  }

  close(progress_bar)

  coexpressed_genes_numbers <- unique(coexpressed_genes_numbers)


  # Get coexpressed genes symbol
  cat("Get coexpressed genes symbol... \n")

  coexpressed_genes <- NULL
  rownames_tcga_matrix_normal <- rownames(tcga_matrix_normal)

  for(row_number in 1:length(rownames_tcga_matrix_normal)){

    if(is.element(row_number, coexpressed_genes_numbers)){

      coexpressed_genes <- c(coexpressed_genes, rownames_tcga_matrix_normal[row_number])

    }
  }

  cat("Total coexpressed: ", length(coexpressed_genes), " \n")

  cat("DONE. \n")

  # Combine coexpressed genes with top_DEG_real
  top_DEG_real <- unique(c(coexpressed_genes, top_DEG_real))

  # Stop CORRELATION MATRIX timer
  tictoc::toc()


  ##########################################################################################################



  # Start CORRESPONDENCE timer
  tictoc::tic("CORRESPONDENCE")

  # Combine top differentially expressed genes

  overlapping_genes <- intersect(top_DEG_simulated, top_DEG_real)
  cat("Overlapping Genes: ", overlapping_genes, "\n")

  top_DEG <- unique(c(top_DEG_simulated, top_DEG_real))
  correspondence <- (length(overlapping_genes) * 100)/length(top_DEG)
  cat(correspondence, " % correspondence in top expressed genes. (", length(overlapping_genes), " overlapping genes)\n")

  # Save correspondence to output file
  correspondence_list <- list(correspondence, length(overlapping_genes), overlapping_genes)
  #write(correspondence_list, file = file.path(output_directory, disease, paste(file_prefix, "_Correspondence.txt", sep="")), sep = "")
  capture.output(correspondence_list, file = file.path(output_directory, disease, paste(file_prefix, "_Correspondence.txt", sep="")))

  # Print influenced genes in top genes
  influenced_genes <- ls(genes_dictionary_from_opentargets)
  cat("Influenced genes in top genes: ", intersect(influenced_genes, top_DEG), "\n")
  cat("Influenced genes in Overlapping genes: ", intersect(influenced_genes, overlapping_genes), "\n")

  # Stop CORRESPONDENCE timer
  tictoc::toc()


  ##########################################################################################################


  # Start BARPLOT timer
  tictoc::tic("BARPLOT")

  cat("Create barplot...")

  # Create a barplot of gene expression rates
  CreateBarplot(disease, output_directory, file_prefix, simulated_matrix, tcga_matrix_tumor, top_DEG, sample_number)

  cat("DONE. \n")

  # Stop BARPLOT timer
  tictoc::toc()


  ##########################################################################################################


  # Start VCF timer
  tictoc::tic("VCF")

  # Generate VCF file
  if(!is.null(genes_variation)){
    SampleGeneVariations(disease, output_directory, genes_variation, file_prefix)
  }
  else{
    cat("Note: There are no gene variations.\n")
  }

  # Stop VCF timer
  tictoc::toc()


  ##########################################################################################################


  cat("----------------------------- \nOMICSIMULATOR DONE. \n----------------------------- \n")

  # Stop TOTAL timer
  tictoc::toc()

  # Return values
  return(list(DEA_normal_tumor, DEA_normal_simulated))

}

