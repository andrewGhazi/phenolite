library(collapse)
library(purrr)
library(jsonlite)

load("data/full_api_req_2026-09-26.RData")

l = resp$body |> rawToChar() |> fromJSON()

nrows = l$data$rows |> sapply(\(x) length(x$fields)) 

api_d = data.table(rw = l$data$rows,
                   nr = nrows,
                   un = l$data$user$username,
                   date = as.IDate(l$data$created))

api_d[, i:= .I]
# ^ Some rare entries have multiple entries ಠ_ಠ . i = entry index

messed_up = api_d |> sbt(nr %!=% 1)

messed_up$rw[[1]]

to_unnest = api_d |> 
  get_elem('rw') |> 
  unlist2d(idcols = "i", DT = TRUE) |> 
  slt(i, fields)

# rbindlist(to_unnest$fields)

to_unnest |> 
  mtt(nrs = sapply(fields, fnrow),
      ncs = sapply(fields, ncol)) |> 
  roworder(nrs) |> 
  fcount(nrs) # variable number of rows...

needed_value = c('species', 'leaf-phenophase', 'flower-phenophase') 
needed_cols = c("slug", "display_value")

to_unnest$fields |> 
  map_lgl(\(x) all(needed_values 
                   %in%
                     x$slug)) |> 
  table()

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
  join(api_d |> slt(i, date, un), 
       verbose = FALSE,
       on = "i")
