
suppressPackageStartupMessages({
  
  library(collapse)
  library(purrr)
  library(data.table)
})
# resp_file = paste0("data/api_resp_", Sys.Date(), ".RData")
# Change this to latest data file to account for runs at midnight...

cli::cli_alert("Beginning to parse API response")

load(resp_file)

# load("data/full_api_req_2026-09-26.RData")

# RcppSimdJson instead?
# l = resp$body |> rawToChar() |> fromJSON()
l = resp$body |> RcppSimdJson::fparse()

nrows = l$data$rows |> sapply(\(x) length(x$fields)) 

un_df = l$data$user |> rbindlist(fill = TRUE) |> qDT()

api_d = data.table(rw = l$data$rows,
                   nr = nrows,
                   # un = l$data$user$username,
                   un = un_df$username,
                   # loc = l$data$location,
                   lat = l$data$location |> purrr::map_dbl("lat"), 
                   lng = l$data$location |> purrr::map_dbl("lng"), 
                   date = as.IDate(l$data$created))

api_d[, i:= .I]
# ^ Some rare entries have multiple entries ಠ_ಠ . i = entry index

messed_up = api_d |> sbt(nr %!=% 1)

# messed_up$rw[[1]]

to_unnest = api_d |>
  get_elem('rw') |>
  unlist2d(idcols = "i", DT = TRUE) |>
  slt(i, fields)
# 
# # rbindlist(to_unnest$fields)
# 
# to_unnest |> 
#   mtt(nrs = sapply(fields, fnrow),
#       ncs = sapply(fields, ncol)) |> 
#   roworder(nrs) |> 
#   fcount(nrs) # variable number of rows...

needed_value = c('species', 'leaf-phenophase', 'flower-phenophase') 

needed_cols = c("slug", "display_value")

has_needed = to_unnest$fields |> 
  map_lgl(\(x) all(needed_value 
                   %in%
                     x$slug)) 

# has_needed |> 
#   table()

stopifnot(all(has_needed))

slt_sbt = \(x) slt(x, needed_cols) |> sbt(slug %iin% needed_value)

tu2 = to_unnest[, fields := 
            map(fields, 
                slt_sbt)][,obs_i := .I]

unnested = tu2$fields |> 
  unlist2d() |> 
  pivot(names = "slug", values = "display_value",
        how = "wider") |> 
  qDT() |> 
  cbind(tu2 |> slt(i, obs_i))

api_out = unnested |>
  join(api_d |> slt(i, date, un, lat, long = lng), 
       verbose = FALSE,
       on = "i")

cli::cli_alert_info("API response contains {fnrow(api_out)} observations")

api_out_file = gsub("api_resp", "api_out", resp_file)

save(api_out, file = api_out_file,
     compress = "zstd", compression_level = 9)

cli::cli_alert_success("Parsed API output written to {api_out_file}")
