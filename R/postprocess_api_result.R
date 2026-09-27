# Given an api output df, add the necessary columns to work with the app.
# Adapted from make_clean_names.R

suppressPackageStartupMessages({
  library(data.table)
  library(collapse)
  library(janitor)
  library(stringr)
  library(purrr)
  library(fs)
})

shiny_data_dir = "~/projects/phenolite/shiny/phenoshiny/data"

proj_dir = "~/projects/phenolite"

out_dir = path(proj_dir, 'output/api_out')

api_out = api_out |> 
  frename(\(x) gsub("-", "_", x))

print(head(api_out))

# add Genus ---------------------------------------------------------------

d_spec = api_out |> 
  slt(species) |> 
  fcount() |> 
  mtt(sp_pcs = strsplit(species, " "),
      n_pc = lengths(sp_pcs))

not_binom = d_spec |> 
  sbt(n_pc %!=% 2) |> 
  slt(-sp_pcs)

not_binom_out = path(out_dir, 'not_binom.tsv')

fwrite(not_binom,
       file = not_binom_out,
       sep = "\t")

cli::cli_alert_danger('Entries with non-binomial names (n = {sum(not_binom$N)}) written to: {not_binom_out}')
cli::cli_alert_warning('Please check those. First few:')

head(not_binom) |> print()

gs_df = d_spec |> 
  sbt(n_pc %==% 2) |> 
  mtt(g = map_chr(sp_pcs, 1) |> str_to_title(),
      s = map_chr(sp_pcs, 2) |> tolower()) |> 
  slt(species, g,s, N) |> 
  funique() |> 
  roworder(g)

# gs_df |> print(nrow = 352)  

gs_df$not_sci = FALSE

# cmn_i = c(28:29, 55, 65, 67, 154, 155, 167, 209, 293, 309, 348:349)
# 
# gs_df$species[cmn_i] |> dput()

# Some non-scientific species names that made it into the species column:
not_sci_nms = c("American witchazel", "American hornbeam", "beaked hazelnut", 
                "Black birch", "canada mayflower", "Golden Alexander", "grey birch", 
                "Highbush blueberry", "mountain laurel", "slippery elm", "Sweet birch", 
                "white oak", "Yellow birch")

setv(gs_df$not_sci,
     gs_df$species %iin% not_sci_nms,
     TRUE)

to_fix = gs_df[(not_sci)]

to_fix_out = path(out_dir,
                  'non_sci.tsv')
fwrite(to_fix, 
       file = to_fix_out,
       sep = "\t")

cli::cli_alert_danger("Non-scientific species names written to: {to_fix_out}")
cli::cli_alert_warning("Please fix those. First few:")

head(to_fix) |> print()

# These are clear:
clr = gs_df[!(not_sci)] |> 
  collapse::frename(species = binom, g = genus, s = species) |> 
  slt(-not_sci) 

gs_out = path(out_dir, 
              'Genus_species.tsv')

fwrite(clr,
       file = gs_out,
       sep = "\t")

# TODO: reactivate
# file_copy(gs_out,
#           shiny_data_dir,
#           overwrite = TRUE)

cli::cli_alert_success("Genus species table written to {gs_out}")

if (grepl("\\/", api_out$date[1])) {
  api_out = api_out |> 
    mtt(date = as.Date(date, format = "%m/%d/%Y"))
}

api_out = api_out |> 
  frename("binom" = "species") |> 
  join(clr, 
       on = "binom",
       verbose = FALSE) |> 
  mtt(wk = lubridate::week(date),
      yr = lubridate::year(date)) 

# add flower phenophase indicators ----------------------------------------

uniq_flw = api_out$flower_phenophase |> 
  unlist2d() |> 
  stringr::str_split(", ") |> 
  unlist() |> 
  funique() |> 
  qDT() |> 
  setColnames("flw_pheno") |> 
  sbt(flw_pheno %!=% "") |> 
  sbt(whichNA(flw_pheno, TRUE)) |> 
  mtt(cln_flw = janitor::make_clean_names(flw_pheno))

uniq_flw_out = path(out_dir, 'uniq_flw.tsv')

fwrite(uniq_flw,
       file = uniq_flw_out,
       sep = "\t")

# TODO: reactivate
# file_copy(uniq_flw_out,
#           shiny_data_dir,
#           overwrite = TRUE)

cli::cli_alert_success("Unique flower phenophases written to: {uniq_flw_out}")

# spaces and parentheses make the column names not syntactically valid

get_flw_j = function(x) {
  uniq_flw$flw_pheno %iin% x
}

idx_df = api_out |> 
  mtt(flw_split = stringr::str_split(flower_phenophase, ", "),
      flw_i = map(flw_split, get_flw_j), 
      i = 1:fnrow(api_out))

nz_idx = slt(idx_df, i, flw_i)[,unlist(flw_i), by = i] |> 
  qM()
# ^ This is probably the most compact way to store this, expand to unique
# columns just for convenience...

flw_mat = matrix(FALSE, nrow = fnrow(api_out),
                 ncol = fnrow(uniq_flw))

flw_mat[nz_idx] = TRUE

colnames(flw_mat) = uniq_flw$cln_flw

api_out = cbind(api_out, qDT(flw_mat))

# add leaf phenophase columns ---------------------------------------------

uniq_lf = api_out$leaf_phenophase |> 
  unlist2d() |> 
  stringr::str_split(", ") |> 
  unlist() |> 
  funique() |> 
  qDT() |> 
  setColnames("lf_pheno") |> 
  sbt(lf_pheno %!=% "") |> 
  sbt(whichNA(lf_pheno, TRUE)) |> 
  mtt(cln_lf = janitor::make_clean_names(lf_pheno))

uniq_lf_out = path(out_dir, 'uniq_lf.tsv')

fwrite(uniq_lf,
       file = uniq_lf_out,
       sep = "\t")

#TODO: reactivate
# file_copy(uniq_lf_out,
#           shiny_data_dir,
#           overwrite = TRUE)

cli::cli_alert_success("Unique flower phenophases written to: {uniq_lf_out}")

# spaces and parentheses make the column names not syntactically valid

get_lf_j = function(x) {
  uniq_lf$lf_pheno %iin% x
}

idx_df = api_out |> 
  mtt(lf_split = stringr::str_split(leaf_phenophase, ", "),
      lf_i = map(lf_split, get_lf_j), 
      i = 1:fnrow(api_out))

nz_idx = slt(idx_df, i, lf_i)[,unlist(lf_i), by = i] |> 
  qM()
# ^ This is probably the most compact way to store this, expand to unique
# columns just for convenience...

lf_mat = matrix(FALSE, nrow = fnrow(api_out),
                ncol = fnrow(uniq_lf))

lf_mat[nz_idx] = TRUE

colnames(lf_mat) = uniq_lf$cln_lf

api_out = cbind(api_out, qDT(lf_mat))

# write output ------------------------------------------------------------

out_f = path(out_dir, 'api_output_clean.csv')

api_out |> 
  frename("species_2" = "binom") |> # Just to be consistent with what app.R expects
  fwrite(file = out_f)

cli::cli_alert_success("Shiny input formatted file written to: {out_f}")

# TODO: reactivate
# dir_create('shiny/phenoshiny/data')
# 
# file_copy(out_f,
#           shiny_data_dir,
#           overwrite = TRUE)

cli::cli_alert_success("Shiny input formatted file copied to: {shiny_data_dir}")

old_f = dir_ls(shiny_data_dir, regexp = 'adj')

if (length(old_f) > 1) {
  
  file_info(old_f) |> 
    roworder(birth_time) |>
    head(1) |> 
    get_elem('path') |> 
    file_delete()
  
  cli::cli_alert_success("Older data file in shiny_data_dir deleted")
}
