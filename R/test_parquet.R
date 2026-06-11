library(arrow)
library(stringr)
library(fastverse)
library(lubridate)
library(ggplot2)
library(fs)

d = fread("data/anecdata_export_EwA_Pheno_Lite_2026-04-27T02-24-48-852Z.csv") |> 
  janitor::clean_names()

d |> 
  qDF() |> 
  write_parquet('data/anecdata_export_EwA_Pheno_Lite_2026-04-27T02-24-48-852Z.parquet', 
                compression = 'gzip')
