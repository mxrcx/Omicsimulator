#' Mapping of disease name to pathway name
#'
#' Mapping of disease name to pathway name
#'
#' @param disease String; name of a disease to get a KEGG pathway for
#'
#' @return String; id of the KEGG pathway
#' @export
#'
#' @examples
GetPathwayID <- function(disease){

  # Map disease name to pathway name
  switch(disease,
         # Cancers: Overview
         'Pathways in cancer'={ pID <- '05200' },
         'Central carbon metabolism in cancer'={ pID <- '05230' },
         'Choline metabolism in cancer'={ pID <- '05231' },
         'Transcriptional misregulation in cancer'={ pID <- '05202' },
         'MicroRNAs in cancer'={ pID <- '05206' },
         'Proteoglycans in cancer'={ pID <- '05205' },
         'Chemical carcinogenesis'={ pID <- '05204' },
         'Viral carcinogenesis'={ pID <- '05203' },

         # Cancers: Specific types
         'Colorectal cancer'={ pID <- '05210' },
         'Pancreatic cancer'={ pID <- '05212' },
         'Hepatocellular carcinoma'={ pID <- '05225' },
         'Gastric cancer'={ pID <- '05226' },
         'Glioma'={ pID <- '05214' },
         'Thyroid cancer'={ pID <- '05216' },
         'Acute myeloid leukemia'={ pID <- '05221' },
         'Chronic myeloid leukemia'={ pID <- '05220' },
         'Basal cell carcinoma'={ pID <- '05217' },
         'Melanoma'={ pID <- '05218' },
         'Renal cell carcinoma'={ pID <- '05211' },
         'Bladder cancer'={ pID <- '05219' },
         'Prostate cancer'={ pID <- '05215' },
         'Endometrial cancer'={ pID <- '05213' },
         'Breast cancer'={ pID <- '05224' },
         'Small cell lung cancer'={ pID <- '05222' },
         'Non-small cell lung cancer'={ pID <- '05223' }
  )

  return(pID)
}
