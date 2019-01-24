#' Main function which calls all relevant to simulate the influence of variations on gene expression
#'
#' Main function which calls all relevant to simulate the influence of variations on gene expression
#'
#' @return None
#' @export
Omicsimulator <-
function(input_file, disease, output_directory, sample_number, top_DEG_number)
{ ... }
#'
#' @param input_file String; txt file which contains gene symbols with variations and gene symbols whose read count is influenced
#' @param disease String; name of a disease to get a KEGG pathway for
#' @param output_directory String; directory of the output files
#' @param top_DEG_number Integer; number of top differential expressed genes
#' @param sample_number Integer; the number of samples which are simulated
#'
#' @examples
#' Omicsimulator(disease = "Breast cancer", sample_number = 10, top_DEG_number = 20)
#' @import fdrtool SimSeq hash stringi vcfR GenomicRanges biomaRt data.table
Omicsimulator <- function(disease, sample_number, top_DEG_number, output_directory, input_file){

  # Check for dependencies
  CheckDependencies()

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

  # Set default disease
  if(missing(disease)){
    disease <- "Breast cancer"
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

  # Create input file [in case none is given as parameter]
  if(missing(input_file)) {

    cat("Create input file from pathway...")

    # Create input file from pathway [in case none is given as parameter]
    pathway_id <- GetPathwayID(disease)

    GenesFromPathway(pID = pathway_id)

    cat("DONE. \n")
  }

  # Get input data

  cat("Read input file...")
  input_data <- GetInputData()
  genes_variation <- input_data$genes_variation
  genes_dictionary <- input_data$genes_dictionary
  print(ls(genes_dictionary))
  cat("DONE. \n")


  ##########################################################################################################


  # Load TCGA Matrix of specified TUMOR data
  sample_type <- "Primary solid Tumor"
  tcga_matrix_tumor <- LoadTCGAMatrix(disease, sample_type, sample_number)

  # Load TCGA Matrix of NORMAL tissue data
  sample_type <- "Solid Tissue Normal"
  tcga_matrix_normal <- LoadTCGAMatrix(disease, sample_type, sample_number)

  # Do differential expression analysis (NORMAL + TUMOR)
  DEA_normal_tumor <- DEA(disease, sample_number, output_directory, tcga_matrix_normal, tcga_matrix_tumor, "NT_PT")
  top_DEG_real <- TopDEG(DEA_normal_tumor, top_DEG_number)


  ##########################################################################################################


  # Simulate tumor
  genes_dictionary_from_dea <- hash::hash(top_DEG_real, rep(1, top_DEG_number))
  simulated_matrix <- SimulateCounts(tcga_matrix_normal, genes_dictionary_from_dea)
  #simulated_matrix <- SimulateCounts(tcga_matrix_normal, genes_dictionary)

  # Do differential expression analysis (NORMAL + SIMULATED)
  DEA_normal_simulated <- DEA(disease, sample_number, output_directory, tcga_matrix_normal, simulated_matrix, "NT_Simulated")
  top_DEG_simulated <- TopDEG(DEA_normal_simulated, top_DEG_number)


  ##########################################################################################################


  cat("Create Correlation Matrices:\n")

  # Correlation Matrix
  cor_normal <- cor(tcga_matrix_normal, method = "pearson", use = "complete.obs")
  cor_tumor <- cor(tcga_matrix_tumor, method = "pearson", use = "complete.obs")
  cor_simulated <- cor(simulated_matrix, method = "pearson", use = "complete.obs")

  print(round(cor_normal, 3))
  print(round(cor_tumor, 3))
  print(round(cor_simulated, 3))

  if(!require(Hmisc)){install.packages("Hmisc")}
  rcorr_simulated <- Hmisc::rcorr(simulated_matrix)

  print(rcorr_simulated)

  cat("DONE. \n")


  ##########################################################################################################


  # Combine top differentially expressed genes
  top_DEG <- unique(c(top_DEG_simulated, top_DEG_real))
  cat((((top_DEG_number * 2) - length(top_DEG)) * 100) / top_DEG_number, " % correspondence in top expresses genes. (", (top_DEG_number * 2) - length(top_DEG), " overlapping genes)\n")

  overlapping_genes <- intersect(top_DEG_simulated, top_DEG_real)
  cat("Overlapping Genes: ", overlapping_genes, "\n")

  influenced_genes <- ls(genes_dictionary)
  cat("Influenced genes in top genes: ", intersect(influenced_genes, top_DEG), "\n")
  cat("Influenced genes in Overlapping genes: ", intersect(influenced_genes, overlapping_genes), "\n")


  ##########################################################################################################


  # Filter both matrices by top differentially expressed genes
  simulated_matrix <- simulated_matrix[row.names(simulated_matrix) %in% top_DEG, ]
  tcga_matrix_tumor <- tcga_matrix_tumor[row.names(tcga_matrix_tumor) %in% top_DEG, ]

  # Merge matrices
  merged_matrix <- merge(simulated_matrix, tcga_matrix_tumor, by = "row.names", all = TRUE)
  row.names(merged_matrix) <- merged_matrix[, 1]
  merged_matrix <- merged_matrix[, 2:ncol(merged_matrix)]

  # Create barplot
  average_count_simulated <- as.integer.Array(rowMeans(merged_matrix[, 1:sample_number]))
  average_count_tcga <- as.integer.Array(rowMeans(merged_matrix[, (sample_number + 1): (2 * sample_number)]))
  average_count <- c(average_count_simulated, average_count_tcga)
  genes <- rownames(merged_matrix)
  type <-c(rep("Simulated PT", length(average_count_simulated)), rep("Real PT", length(average_count_tcga)))

  data <-data.frame(genes, average_count)
  p <- ggplot2::ggplot(data, aes(genes, average_count)) +
    geom_bar(stat = "identity", aes(fill = type), position = "dodge") +
    xlab("Genes") +
    ylab("Average expression count") +
    ggtitle("Gene expression of top differentially expressed genes") +
    theme_bw()

  print(p)
  ggsave(file.path(output_directory, disease, paste("Top_DGE_barplot_", disease, "_sample=", sample_number, ".png", sep="")), height = 50, width = 100, units = "cm", limitsize = FALSE)

  # Merge both Matrices
  # merged_matrix <- MergeMatrices(disease, sample_number, output_directory, simulated_matrix, tcga_matrix)

  # Do DGE analysis
  # DGE(disease, sample_number, output_directory, merged_matrix)

  # Generate VCF file
  if(!is.null(genes_variation)){
    SampleGeneVariations(disease, output_directory, genes_variation)
  }
  else{
    cat("Note: There are no gene variations.\n")
  }

  cat("----------------------------- \nOMICSIMULATOR DONE. \n----------------------------- \n")
}
