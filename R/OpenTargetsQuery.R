# R-Script - Accessing associations from opentargets.org


###########################


# Get EFO id for a disease

# 1. GET request

search_term <- "breast"

url_efo<- paste("http://www.ebi.ac.uk/ols/api/search?q=", search_term, "&rows=10000", sep = "")

res_efo <- httr::GET(url = url_efo)

# 2. Transform JSON into a list

res_json_efo <- httr::content(res_efo, as = "text")

res_list_efo <- jsonlite::fromJSON(res_json_efo)

# 3. Extract data

res_terms_efo <- res_list_efo$response$docs

efo_list <- subset(res_terms_efo, ontology_prefix == 'EFO', select = c("ontology_name", "id", "label"))

print(efo_list)


###########################


# Access opentargets.org REST style API and download data

# 1. GET request

efo_code <- "EFO_0000305"

url_opentargets <- paste("https://api.opentargets.io/v3/platform/public/association/filter?disease=", efo_code, "&size=1000", sep = "")

res_opentargets <- httr::GET(url = url_opentargets)


# 2. Transform JSON into a list

res_json <- httr::content(res_opentargets, as = "text")

res_list <- jsonlite::fromJSON(res_json)


# 3. Extract data

res_data <- res_list$data

res_selection <- cbind.data.frame(res_data$target$gene_info$symbol, res_data$association_score$overall)

# 4. Get top genes

treshold <- 100

top_genes <- res_selection[1:treshold, 1]

print(top_genes)
