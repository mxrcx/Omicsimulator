#' Generates a vcf file based on a given list of gene symbols
#'
#' Gets gene range of every gene with variations, Searches for possible variations in this range, Writes vcf entry for first variation. Prints complete vcf file.
#'
#' @param disease
#' @param output_directory
#' @param genes_variation
#'
#' @return None
#' @export
#'
#' @examples
#' SampleGeneVariations(genes.variation)
SampleGeneVariations <- function(disease, output_directory, genes_variation){

  cat("Calculate gene variations...")

  if(file.exists(file.path(output_directory, disease, paste("gene_variations_", disease, ".vcf.gz", sep="")))){

    cat("SAVE FILE ALREADY EXISTING. \n")

  }
  else {

    # Read a given vcf file
    vcf_file <- system.file("extdata", "vcf_input.vcf", package = "omicsimulator")
    vcf <- vcfR::read.vcfR(vcf_file, verbose = FALSE)

    # Choose the size of variants and samples
    vcf <- vcf[1:length(genes_variation), 1:11]

    # Get gene ranges
    mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "may2009.archive.ensembl.org", path = "/biomart/martservice", archive = FALSE)
    gene_ranges <- biomaRt::getBM(attributes = c("hgnc_symbol", "chromosome_name",
                                        "start_position", "end_position"),
                         filters = c("hgnc_symbol"),
                         values=list(genes_variation),
                         mart=mart)

    snpmart <- biomaRt::useMart(biomart = "ENSEMBL_MART_SNP", dataset="hsapiens_snp")

    # Create data table
    data <- data.table::data.table(v1=numeric(), v2=numeric(), v3=numeric(), v4=numeric(), v5=numeric(), v6=numeric(), v7=numeric(), v8=numeric())

    chromosome_names <- NULL
    start_positions <- NULL
    end_positions <- NULL
    entries <- NULL
    for(row in 1:nrow(gene_ranges)){
      # Get variations
      chromosome_name <- gene_ranges$chromosome_name[row]
      start_position <- gene_ranges$start_position[row]
      end_position <- start_position + 100 #gene.ranges$end_position[row] - 0.9 * (gene.ranges$end_position[row] - gene.ranges$start_position[row])

      chromosome_names <- c(chromosome_name, chromosome_names)
      start_positions <- c(start_position, start_positions)
      end_positions <- c(end_position, end_positions)
    }

    pos_values <- list(chromosome_names, start_positions, end_positions)


    #variations <- biomaRt::getBM(attributes = c('refsnp_id','allele','chrom_start'),
     #                            filters = c('chr_name','start','end'),
      #                           values = pos_values,
       #                          mart = snpmart)
    #print(variations)

    # Manipulate the vcf and change the given values to own generated values
    for(row in 1:nrow(gene_ranges)){
      cat("\n \t Processing variation ",row ," of ", nrow(gene_ranges), "...")

      # Chromosome
      chrom <- gene_ranges$chromosome_name[row]

      # Get variations
      chromosome_name <- gene_ranges$chromosome_name[row]
      start_position <- gene_ranges$start_position[row]
      end_position <- start_position + 100 #gene.ranges$end_position[row] - 0.9 * (gene.ranges$end_position[row] - gene.ranges$start_position[row])

      value <- list(chromosome_name, start_position, end_position)
      values <- c(values, value)

      variations <- biomaRt::getBM(attributes = c('refsnp_id','allele','chrom_start'),
                                   filters = c('chr_name','start','end'),
                                   values = list(chromosome_name, start_position, end_position),
                                   mart = snpmart)

      # Position
      pos <- variations$chrom_start[row]

      # Variation id
      id <- variations$refsnp_id[row]

      ref <- stringi::stri_extract(variations$allele[row], regex='[^/]*') # reference
      alt <- stringi::stri_extract(variations$allele[row], regex='[^/]+$') # alternative

      # deafult columns
      qual <- 100
      filter <- 'PASS'
      info <- 'None'

      # Create data table entry
      entry <- cbind(chrom, POS=pos, ID=id, REF=ref, ALT=alt, QUAL=qual, FILTER=filter, INFO=info)
      if(nrow(data) == 0)
        data <- rbind(entry)
      else
        data <- rbind(data,entry)

      # Appearance in sample genome
      for(col in 2:6){
        vcf@gt[row, col] <- "0|0"
      }
      for(col in 7:11){
        vcf@gt[row, col] <- "1|1"
      }

      cat("DONE.")
    }

    vcf@fix <- data

    # Write manipulated vcf file
    vcfR::write.vcf(vcf, file = file.path(output_directory, disease, paste("gene_variations_", disease, ".vcf.gz", sep="")))

    cat("\nDONE. \n")
  }
}
