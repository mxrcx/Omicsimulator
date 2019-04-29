#' Get gene expression relations between genes from a specific pathway
#'
#' Get gene expression relations between genes from a KEGG pathway and write them in input file, using pathway specific kgml file of KEGGREST database
#'
#' @param pID String; id of the KEGG pathway
#' @param species
#'
#' @export
#'
#' @examples
GenesFromPathway <- function(species, pathway_id){

  # Default species value
  if(missing(species)) {
    species <- "hsa"
  }

  # Get kgml file & create pathway
  kgml.file <- KEGGREST::keggGet(paste(species, pID, sep=""), c("kgml"))
  kgml.pathway <- KEGGgraph::parseKGML(kgml.file)

  # Create nodes and edges objects
  pNodes <- KEGGgraph::nodes(kgml.pathway)
  pEdges <- KEGGgraph::edges(kgml.pathway)

  # Filter edges for gene expression relations & create influnce file entries
  entry.number <- 1
  entries <- list()

  for(entry in 1:length(pEdges)){
    edge <- pEdges[[entry]]
    if(KEGGgraph::getType(edge) == "GErel"){
      if(KEGGgraph::getName(KEGGgraph::getSubtype(edge)[[1]]) == "expression"){

        gene1.symbol <- sub('\\,.*', '', KEGGgraph::getDisplayName(pNodes[[KEGGgraph::getEntryID(edge)[1]]]))
        gene2.symbol <- sub('\\,.*', '', KEGGgraph::getDisplayName(pNodes[[KEGGgraph::getEntryID(edge)[2]]]))

        entries[[entry.number]] <- list(gene1.symbol, gene2.symbol, "+") #paste(gene1.symbol, gene2.symbol, "+", sep=" ")
        entry.number <- entry.number + 1

      }
      if(KEGGgraph::getName(KEGGgraph::getSubtype(edge)[[1]]) == "repression"){

        gene1.symbol <- sub('\\,.*', '', KEGGgraph::getDisplayName(pNodes[[KEGGgraph::getEntryID(edge)[1]]]))
        gene2.symbol <- sub('\\,.*', '', KEGGgraph::getDisplayName(pNodes[[KEGGgraph::getEntryID(edge)[2]]]))

        entries[[entry.number]] <- list(gene1.symbol, gene2.symbol, "-") #paste(gene1.symbol, gene2.symbol, "-", sep=" ")
        entry.number <- entry.number + 1

      }
    }
  }


  # Write entries into input file
  #file.create('inst/extdata/influence_of_variations.txt')
  #file.create(system.file("extdata",'input_file.txt', package = "omicsimulator"))

  #for (entry in entries) {
    #write(entry,file=('inst/extdata/influence_of_variations.txt'),append=TRUE)
    #write(entry,file=system.file("extdata",'input_file.txt', package = "omicsimulator"),append=TRUE)
  #}

  saveRDS(entries, file.path("cache", "inputFile_influence.rds"))

  # Visualize Pathway
  # VisualizePathway(species, pID)
}
