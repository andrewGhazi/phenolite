
    # This is the example, adapt as needed.

suppressPackageStartupMessages({
  library(httr2)
  library(purrr) # For better list handling
  # library(jsonlite)
  library(collapse)
  
})

# Configuration
api_config <- list(
  base_url = "https://api.anecdata.org",
  # Replace with a real key from the Developer page
  api_key = readLines("api.env") 
)

req <- request(api_config$base_url) %>%
  req_url_path("observations") %>%
  req_headers(
    "X-Api-Key" = api_config$api_key,
    "Content-Type" = "application/json",
    "Accept" = "application/json"
  ) %>%
  req_url_query( "project_id" = "ewa-pheno-lite",
                 "limit" = 1e6) |> # Doesn't seem to be a way to make this Inf?
  req_retry(max_tries = 3)

get_resp = function() {
  
  cli::cli_progress_step("Running API request...")
  
  tryCatch({
    resp <- req_perform(req)
    cli::cli_alert_success("Request succesful.")
  }, http_401 = function(cnd) {
    cat("Status: 401 Unauthorized\n")
    cat("Message: Invalid API Key or inactive application.\n")
  }, error = function(cnd) {
    cli::cli_alert_danger("An error occurred:\n")
    print(cnd)
  })
  
}
resp_file = path(data_dir, paste0("api_resp_", Sys.Date()), ext = "RData")

save(resp,
     file = resp_file,
     compress = "zstd",
     compression_level = 9)

cli::cli_alert_success("API response saved to {resp_file}")
