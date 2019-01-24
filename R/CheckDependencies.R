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
  if(!require(httr)){install.packages("httr")}
  if(!require(R.utils)){install.packages("R.utils")}
  if(!require(vcfR)){install.packages("vcfR")}
  if(!require(fdrtool)){install.packages("fdrtool")}
  if(!require(SimSeq)){install.packages("SimSeq")}
  if(!require(stringi)){install.packages("stringi")}
  if(!require(devtools)){install.packages("devtools")}
  if(!require(pheatmap)){install.packages("pheatmap")}
  if(!require(Rgraphviz)){
    source("https://bioconductor.org/biocLite.R")
    biocLite("Rgraphviz")
  }
  if(!require(CePa)){install.packages("CePa")}
  if(!require(biomaRt)){
    source("https://bioconductor.org/biocLite.R")
    biocLite("biomaRt")
  }
  if(!require(KEGGREST)){
    source("https://bioconductor.org/biocLite.R")
    biocLite("KEGGREST")
  }
  if(!require(KEGGgraph)){
    source("https://bioconductor.org/biocLite.R")
    biocLite("KEGGgraph")
  }
  if(!require(TCGAbiolinks)){
    #source("https://bioconductor.org/biocLite.R")
    #biocLite("TCGAbiolinks")
    devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")
  }
  if(!require(DESeq2)){
    source("https://bioconductor.org/biocLite.R")
    biocLite("DESeq2")
  }
  if(!require(SummarizedExperiment)){
    source("https://bioconductor.org/biocLite.R")
    biocLite("SummarizedExperiment")
  }
  if(!require(ggplot2)){
    install.packages("ggplot2")
  }
  if(!require("hpar")){
    if (!requireNamespace("BiocManager", quietly = TRUE)){
      install.packages("BiocManager")
    }
    BiocManager::install("hpar", version = "3.8")
  }
}
