#' Get Cancer project from TCGA
#'
#' @param disease
#'
#' @return String; code for the TCGA Cancer Project
#' @export
#'
#' @examples
GetCancerProject <- function(disease){

  # Map disease name to TCGA Cancer Project
  switch(disease,
         # Cancers: Overview
         'Pathways in cancer'={ cancerProject <- '05200' },
         'Central carbon metabolism in cancer'={ cancerProject <- '05230' },
         'Choline metabolism in cancer'={ cancerProject <- '05231' },
         'Transcriptional misregulation in cancer'={ cancerProject <- '05202' },
         'MicroRNAs in cancer'={ cancerProject <- '05206' },
         'Proteoglycans in cancer'={ cancerProject <- '05205' },
         'Chemical carcinogenesis'={ cancerProject <- '05204' },
         'Viral carcinogenesis'={ cancerProject <- '05203' },

         # Cancers: Specific types
         'Colorectal cancer'={ cancerProject <- '05210' },
         'Pancreatic cancer'={ cancerProject <- 'TCGA-PAAD' },
         'Hepatocellular carcinoma'={ cancerProject <- '05225' },
         'Gastric cancer'={ cancerProject <- '05226' },
         'Glioma'={ cancerProject <- '05214' },
         'Thyroid cancer'={ cancerProject <- '05216' },
         'Acute myeloid leukemia'={ cancerProject <- '...' },
         'Chronic myeloid leukemia'={ cancerProject <- '...' },
         'Basal cell carcinoma'={ cancerProject <- '05217' },
         'Melanoma'={ cancerProject <- 'TCGA-SKCM' },
         'Renal cell carcinoma'={ cancerProject <- '05211' },
         'Bladder cancer'={ cancerProject <- 'TCGA-BLCA' },
         'Prostate cancer'={ cancerProject <- '05215' },
         'Endometrial cancer'={ cancerProject <- '05213' },
         'Breast cancer'={ cancerProject <- 'TCGA-BRCA' },
         'Small cell lung cancer'={ cancerProject <- '05222' },
         'Non-small cell lung cancer'={ cancerProject <- '05223' }
  )

  return(cancerProject)
}
