#' Create a barplot of gene expression rates
#'
#' @param disease String; cancer name
#' @param output_directory String; name of the output directory
#' @param file_prefix String; name of the results output file
#' @param simulated_matrix
#' @param tcga_matrix_tumor
#' @param top_DEG
#' @param sample_number Integer; the number of samples which are simulated
#'
#' @return
#' @export
#'
#' @examples
CreateBarplot <- function(disease, output_directory, file_prefix, simulated_matrix, tcga_matrix_tumor, top_DEG, sample_number){

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
  ggplot2::ggsave(file.path(output_directory, disease, paste(file_prefix, "_TopDGE_Barplot", sep="")), device = "png", height = 50, width = 100, units = "cm", limitsize = FALSE)

}
