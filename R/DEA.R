#' DEA script
#'
#' @param disease
#' @param sample_number
#' @param output_directory
#' @param count_matrix
#' @param description
#'
#' @return
#' @export
#'
#' @examples
DEA <- function(disease, sample_number, output_directory, count_matrix_control, count_matrix_test, description){

  cat("Do Differential Expression Analysis: \n")

  if(!file.exists(file.path("cache", disease, paste("DEA_results", disease, "_description=", description, "_sample=", sample_number, ".rds", sep="")))){

    # Merge Matrices
    count_matrix <- merge(count_matrix_control[, 1:sample_number], count_matrix_test[, 1:sample_number], by = "row.names", all = TRUE)
    count_matrix <- count_matrix[complete.cases(count_matrix), ]

    # Set correct rownames for merged Matrix
    row.names(count_matrix) <- count_matrix[, 1]
    count_matrix <- count_matrix[, 2:ncol(count_matrix)]
    count_matrix <- as.matrix(count_matrix)

    # Extract Sample-Names
    names <- colnames(count_matrix)

    # Extract Gene-Symbols
    genes <- rownames(count_matrix)

    # Choose Groups
    group <- as.factor(rep(c("Control", "Test"), times = c(sample_number, sample_number)))
    #group <- as.factor(rep(1, times = sample_number))

    # Create col_data
    col_data <- data.frame(samples = names, group = group)

    # IMPORTANT! Change mode of matrix --> only Integer as counts allowed
    mode(count_matrix) <- "integer"

    # Create ddsMat
    ddsMat <- DESeq2::DESeqDataSetFromMatrix(countData = count_matrix, colData = col_data, design = ~ group)

    # Filter samples with low expression
    dds_filtered <- ddsMat[ , colSums(counts(ddsMat)) > 5]

    # Filter genes with low expression
    dds_filtered <- dds_filtered[rowSums(counts(dds_filtered)) > 1, ]

    # Normalization
    dds_filtered = DESeq2::estimateSizeFactors(dds_filtered)

    # Get p-Values
    prep <- DESeq2::DESeq(dds_filtered, fitType = 'local', quiet = TRUE, parallel = TRUE)
    res <- DESeq2::results(prep)
    print(res)

    # Write DGE analysis data to csv files
    if(!file.exists(file.path(output_directory, disease, description))){
      dir.create(file.path(output_directory, disease, description))
    }
    count_matrix <- cbind(count_matrix, res$pvalue)
    write.table(count_matrix, file = file.path(output_directory, disease, description, paste(description, "_count_matrix_", disease, "_sample=", sample_number, ".csv", sep="")), quote=FALSE, sep =";", row.names = TRUE, col.names = NA)
    write.table(res, file = file.path(output_directory, disease, description, paste(description, "_DGE_results_", disease, "_sample=", sample_number, ".csv", sep="")), quote=FALSE, sep =";", row.names = TRUE, col.names = NA)

    # log-Transforamtion (vst instead of rlog, because it is faster (vst needs >= 1000 rows))
    #dds_transformed = DESeq2::varianceStabilizingTransformation(dds_filtered, fitType = 'local')

    # Build heatmap
    #mat  <- SummarizedExperiment::assay(dds_transformed)
    #mat  <- mat - rowMeans(mat)
    #annotation_c <- as.data.frame(SummarizedExperiment::colData(dds_transformed)[,c("group")])
    #colnames(annotation_c) <- c("Group")
    #rownames(annotation_c) <- colnames(dds_transformed)
    #pheatmap::pheatmap(mat, annotation_col = annotation_c, filename = file.path(output_directory, paste(description, "_pheatmap_exprVal_", disease, "_sample=", sample_number, ".pdf", sep="")))

    # Save as RDS file in cache
    saveRDS(res, file.path("cache", disease, paste("DEA_results", disease, "_description=", description, "_sample=", sample_number, ".rds", sep="")))

    cat("--> DONE. \n")
  }
  else{

    res <- readRDS(file.path("cache", disease, paste("DEA_results", disease, "_description=", description, "_sample=", sample_number, ".rds", sep="")))

    cat("Loaded from save file.\n")
  }

  return(res)
}
