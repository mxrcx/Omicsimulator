#' Main function which calls all relevant to simulate the influence of variations on gene expression
#'
#' Main function which calls all relevant to simulate the influence of variations on gene expression
#'
#' @return List; List containing the DEA results of both normal-tumor and normal-simulated
#' @export
Omicsimulator <-
  function(disease, sample_number, top_DEG_number, output_directory, file_prefix, threshold_eQTls)
  { ... }
#'
#' @param disease String; name of a disease to get a KEGG pathway for
#' @param output_directory String; directory of the output files
#' @param top_DEG_number Integer; number of top differential expressed genes
#' @param sample_number Integer; the number of samples which are simulated
#' @param file_prefix String; name of the results output file
#' @param threshold_eQTls Integer; Percentage of used eQTLs
#'
#' @examples
#' Omicsimulator(disease = "Breast cancer", sample_number = 10, top_DEG_number = 100, threshold_eQTls = 0.7)
Omicsimulator <- function(disease, sample_number, top_DEG_number, output_directory, file_prefix, threshold_eQTls){

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

  # Set threshold for eQTLs [in case none is given as parameter]
  if(missing(threshold_eQTls)) {
    threshold_eQTls <- 0.7
  }


  ##########################################################################################################


  cat("----------------------------- \nOMICSIMULATOR START. \n----------------------------- \n")


  # Load TCGA Matrix of specified TUMOR data
  sample_type <- "Primary solid Tumor"
  tcga_matrix_tumor <- LoadTCGAMatrix(disease, sample_type, sample_number)

  # Load TCGA Matrix of NORMAL tissue data
  sample_type <- "Solid Tissue Normal"
  tcga_matrix_normal <- LoadTCGAMatrix(disease, sample_type, sample_number)

  # Do differential expression analysis (NORMAL + TUMOR)
  DEA_normal_tumor <- DEA(disease, sample_number, output_directory, tcga_matrix_normal, tcga_matrix_tumor, "NT_PT", file_prefix)
  top_DEG_real <- TopDEG(DEA_normal_tumor, top_DEG_number)

  # Generate MAF File
  tumor_sample_barcodes <- colnames(tcga_matrix_normal)
  simulated_matrix <- GenerateMAF(threshold_eQTls = threshold_eQTls, tumor_sample_barcodes = tumor_sample_barcodes, output_directory = output_directory, file_prefix = file_prefix, disease = disease, tcga_matrix_normal = tcga_matrix_normal, tcga_matrix_tumor = tcga_matrix_tumor)

  # Do differential expression analysis (NORMAL + SIMULATED)
  DEA_normal_simulated <- DEA(disease, sample_number, output_directory, tcga_matrix_normal, simulated_matrix, "NT_Simulated", file_prefix)
  top_DEG_simulated <- TopDEG(DEA_normal_simulated, top_DEG_number)


  ##########################################################################################################


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
  #influenced_genes <- ls(genes_dictionary_from_opentargets)
  #cat("Influenced genes in top genes: ", intersect(influenced_genes, top_DEG), "\n")
  #cat("Influenced genes in Overlapping genes: ", intersect(influenced_genes, overlapping_genes), "\n")


  ##########################################################################################################


  cat("Create barplot...")

  # Create a barplot of gene expression rates
  CreateBarplot(disease, output_directory, file_prefix, simulated_matrix, tcga_matrix_tumor, top_DEG, sample_number)

  cat("DONE. \n")


  ##########################################################################################################

  genes_variation <- NULL

  # Generate VCF file
  if(!is.null(genes_variation)){
    SampleGeneVariations(disease, output_directory, genes_variation, file_prefix)
  }
  else{
    cat("Note: There are no gene variations.\n")
  }


  ##########################################################################################################


  cat("----------------------------- \nOMICSIMULATOR DONE. \n----------------------------- \n")

  # Stop TOTAL timer
  tictoc::toc()

  # Return values
  return(list(DEA_normal_tumor, DEA_normal_simulated))

}

