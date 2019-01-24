#' Read given gene symbols out of txt file
#'
#' Reads data out of txt file and puts it into a list, containing gene symbols of genes with variation, and a dictionary,
#' containing gene symbols of genes mapped to their influence on the expression value.
#' Every row of the txt file has to be formatted this way:
#' GeneSymbolOfGeneWithVariation[String] GeneSymbolOfGeneWithChangedExpressionValue[String] Influence[Characters --,-,++,+]
#'
#'
#' @return List of all gene symbols with variations and dictionary of gene symbols and the specific influence on their expression
#' @export
GetInputData <-
function()
{ ... }
#'
#' @examples
#' GetInputData()
GetInputData <- function(){

  # Read input out of rds file
  entries <- readRDS(file.path("cache", "inputFile_influence.rds"))

  # Create dictionary
  genes_variation <- NULL
  genes_influenced <- NULL
  influence <- NULL
  genes_dictionary <- NULL

  if(length(entries) == 0){

    cat("Note: There is no influence on this dataset.\n")

  }
  else{

    # Fill dictionary with input data
    input_data <- matrix(unlist(entries), ncol = 3, byrow = TRUE)

    for(row in 1:nrow(input_data)) {
      genes_variation <- c(genes_variation, toString(input_data[row, 1]))
      genes_influenced <- c(genes_influenced, toString(input_data[row, 2]))
      if(input_data[row, 3] == "+"){
        influence <- c(influence, 1)
      }
      else if(input_data[row, 3] == "++"){
        influence <- c(influence, 2)
      }
      else if(input_data[row, 3] == "-"){
        influence <- c(influence, -1)
      }
      else{
        influence <- c(influence, -2)
      }
    }

    genes_dictionary <- hash::hash(genes_influenced, influence)
  }

  return(list("genes_variation" = genes_variation, "genes_dictionary" = genes_dictionary))
}
