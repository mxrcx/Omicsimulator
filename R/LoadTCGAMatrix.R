#' Load the TCGA Matrix of specified disease
#'
#' @param sample_number Integer; the number of samples which are simulated
#' @param disease String; cancer name
#' @param sample_type String; type of samples ('Primary solid Tumor' or 'Solid Tissue Normal')
#'
#' @return tcga_matrix
#' @export
#'
#' @examples
LoadTCGAMatrix <- function(disease, sample_type, sample_number){

  cat("Get TCGA Matrix...")

  files <- list.files(path = file.path("cache", disease), pattern = paste("tcga_", disease, "_", sample_type, "_sample=.*\\.rds", sep=""))
  print(files)

  if(length(files) > 0){

    # Load stored data from file
    tcga_matrix <- readRDS(file.path("cache", disease, files))

    if(ncol(tcga_matrix) >= sample_number){

      tcga_matrix <- tcga_matrix[, 1:sample_number]

      return(tcga_matrix)

    }
    else{

      file.remove(file.path("cache", disease, files))

    }
  }

  # New TCGA Query
  cancerProject <- GetCancerProject(disease)

  query <- TCGAbiolinks::GDCquery(project = cancerProject,
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "HTSeq - Counts",
                    sample.type = sample_type)

  samples <- getResults(query,cols=c("cases"))

  query <- TCGAbiolinks::GDCquery(project = cancerProject,
                                  data.category = "Transcriptome Profiling",
                                  data.type = "Gene Expression Quantification",
                                  workflow.type = "HTSeq - Counts",
                                  #barcode = c(dataSmTP_short, dataSmNT_short)
                                  barcode = samples[1:sample_number])

  TCGAbiolinks::GDCdownload(query)

  tcga_rna_seq_data <- TCGAbiolinks::GDCprepare(query)

  tcga_matrix <- SummarizedExperiment::assay(tcga_rna_seq_data) # , "raw_count")

  # Transform Ensemble IDs to Gene Symbols
  mart <- biomaRt::useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  genes_ensemble <- rownames(tcga_matrix)
  bm_result <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id","hgnc_symbol"), values = genes_ensemble, mart = mart)

  for(pos in 1:length(genes_ensemble)){
    if(complete.cases(bm_result[pos:pos, ])){
      genes_ensemble[pos] <- bm_result[pos, 'hgnc_symbol']
    }
  }
  rownames(tcga_matrix) <- genes_ensemble

  # Save as RDS file in cache
  saveRDS(tcga_matrix, file.path("cache", disease, paste("tcga_", disease, "_", sample_type, "_sample=", sample_number, ".rds", sep="")))

  return(tcga_matrix)

}
