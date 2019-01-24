#' Return top differentially expressed genes
#'
#' @param results_DEG Matrix; differentaial analysis results
#' @param top_DEG_number Integer; number of top differential expressed genes
#'
#' @return
#' @export
#'
#' @examples
TopDEG <- function(results_DEG, top_DEG_number){

  # Order results by pvalue
  resOrdered <- results_DEG[order(results_DEG$pvalue),]
  print(resOrdered)

  # Extract top differential expressed genes
  top_DEG <- rownames(resOrdered)[1:top_DEG_number]

  return(top_DEG)
}
