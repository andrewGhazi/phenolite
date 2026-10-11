# This script will update the dataset.  The last step deploy_app.R is not
# commited b/c it has creds, but it only runs two commands 1) connects to the
# EwA Connect Cloud account with rsconnect::connectCloudClientCredentials() with
# an appropriate credential and 2) runs rsconnect::deployApp() to republish the
# app with the new dataset

cli::start_app(output = "stdout")

library(fs)

proj_dir = "~/projects/phenolite"

setwd(proj_dir)

r_dir = path(proj_dir, "R")

data_dir = path(proj_dir, "data")

source(path(r_dir, "run_api_req.R"))

source(path(r_dir, "process_api_result.R"))

source(path(r_dir, "add_api_columns.R"))

source(path(r_dir, "deploy_app.R"))
