# This is the example, adapt as needed.

library(httr2)
library(purrr) # For better list handling

# Configuration
api_config <- list(
  base_url = "https://api.anecdata.org",
  # Replace with a real key from the Developer page
  api_key = readLines("api.env") 
)

fetch_observations <- function() {
  cat("--- Starting API Request ---\n")
  
  # Construct the request
  req <- request(api_config$base_url) %>%
    req_url_path("observations") %>%
    req_headers(
      "X-Api-Key" = api_config$api_key,
      "Content-Type" = "application/json",
      "Accept" = "application/json"
    ) %>%
    req_url_query( "project_id" = "ewa-pheno-lite",
                   "limit" = 8029) |> 
    req_retry(max_tries = 3)
  
  # Perform the request
  tryCatch({
    resp <- req_perform(req)
    
    # Parse the full JSON response
    response_body <- resp %>% resp_body_json()
    
    # Extract data and pagination
    observations <- response_body$data
    pagination <- response_body$pagination
    
    cat("Success! Found", pagination$totalItems, "total items.\n")
    cat("Showing current page:", pagination$currentPage, "\n\n")
    
    # Convert the 'data' part to a Data Frame
    # We use map_dfr from purrr to row-bind list elements into a data frame
    df <- map_dfr(observations, as.data.frame)
    
    print(head(df))
    return(df)
    
  }, http_401 = function(cnd) {
    cat("Status: 401 Unauthorized\n")
    cat("Message: Invalid API Key or inactive application.\n")
  }, error = function(cnd) {
    cat("An error occurred:\n")
    print(cnd)
  })
}

# Execute the function
observation_df <- fetch_observations()
