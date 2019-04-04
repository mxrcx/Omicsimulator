#' CheckDependencies
#'
#' Check Dependecies needed for this project (install if neccessary + load dependencies)
#'
#' @return
#' @export
#'
#' @examples
CheckDependencies <- function(){

  if(!require(hash)){install.packages("hash")}
  if(!require(jsonlite)){install.packages("jsonlite")}
  if(!require(rlist)){install.packages("rlist")}
  if(!require(propagate)){install.packages("propagate")}
  if(!require(httr)){install.packages("httr")}
  if(!require(R.utils)){install.packages("R.utils")}
  if(!require(vcfR)){install.packages("vcfR")}
  if(!require(stringi)){install.packages("stringi")}
  if(!require(devtools)){install.packages("devtools")}
  if(!require(pheatmap)){install.packages("pheatmap")}
  if(!require(CePa)){install.packages("CePa")}
  if(!require(ggplot2)){install.packages("ggplot2")}
  if(!require(functional)){install.packages("functional")}

  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }

  if(!require(Rgraphviz)){
    BiocManager::install("Rgraphviz")
  }
  if(!require(biomaRt)){
    BiocManager::install("biomaRt")
  }
  if(!require(KEGGREST)){
    BiocManager::install("KEGGREST")
  }
  if(!require(KEGGgraph)){
    BiocManager::install("KEGGgraph")
  }
  if(!require(TCGAbiolinks)){
    BiocManager::install("TCGAbiolinks")
    #devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")
  }
  if(!require(DESeq2)){
    BiocManager::install("DESeq2")
  }
  if(!require(SummarizedExperiment)){
    BiocManager::install("SummarizedExperiment")
  }

  if(!require("hpar")){
    BiocManager::install("hpar")
  }
  if(!require("rols")){
    BiocManager::install("rols")
  }
}
