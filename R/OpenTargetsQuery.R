# Access opentargets.org REST style API and download data

# 1. GET request

url_opentargets <- "https://api.opentargets.io/v3/platform/public/association/filter?disease=EFO_0000305&size=1000"

res_opentargets <- httr::GET(url = url_opentargets)


# 2. Transform JSON into a list

res_json <- httr::content(res_opentargets, as = "text")

res_list <- jsonlite::fromJSON(res_json)


# 3. Extract data

res_data <- res_list$data

print(res_data)




