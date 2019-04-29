#' GetTopGenesFromOpentargets
#'
#' Access opentargets.org REST style API and download data
#'
#' @param efo_code String; coding of a specific cancer disease (e.g. Breast cancer: EFO_0000305)
#' @param top_DEG_number Integer; number of top differential expressed genes
#'
#' @return top genes
GetTopGenesFromOpentargets <- function(efo_code, top_DEG_number){

  # 1. GET request

  url_opentargets <- paste("https://api.opentargets.io/v3/platform/public/association/filter?disease=", efo_code, "&size=", top_DEG_number, sep = "")

  res_opentargets <- httr::GET(url = url_opentargets)


  # 2. Transform JSON into a list

  res_json <- httr::content(res_opentargets, as = "text")

  res_list <- jsonlite::fromJSON(res_json)


  # 3. Extract data

  res_data <- res_list$data

  res_selection <- cbind.data.frame(res_data$target$gene_info$symbol, res_data$association_score$overall)

  # 4. Get top genes

  treshold <- top_DEG_number

  top_genes <- res_selection[1:treshold, 1]

  return(top_genes)
}
