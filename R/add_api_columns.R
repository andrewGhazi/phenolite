# Go through the api output file and add the columns the app needs like make_clean_names.R

cli::cli_alert_info("Beginning to add app-required columns to parsed API output...")

library(fs)
library(stringr)

api_out_file = gsub("api_resp", "api_out", resp_file)

load(api_out_file)

proj_dir = "~/projects/phenolite"

shiny_data_dir = path(proj_dir, "shiny/phenoshiny/data")

# out_f = f |> 
#   path_real() |> 
#   gsub("anecdata_export_EwA_Pheno_Lite",
#        "adj",
#        x = _)
# 
# d = fread(f) |> 
#   janitor::clean_names()

d = api_out |> janitor::clean_names()

# print(head(d))

# add Genus ---------------------------------------------------------------

d_spec = d |> 
  fcount() |> 
  mtt(sp_pcs = strsplit(species, " "),
      n_pc = lengths(sp_pcs))

not_binom = d_spec |> 
  sbt(n_pc %!=% 2) |> 
  slt(-sp_pcs)

not_binom_out = path(proj_dir, 'output/not_binom.tsv')

fwrite(not_binom,
       file = not_binom_out,
       sep = "\t")

cli::cli_alert_danger('Entries with non-binomial names written to: {not_binom_out}')
cli::cli_alert_warning('Please check those. First few:')

head(not_binom |> slt(species)) |> print()

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

to_fix_out = path(proj_dir,
                  'output/non_sci.tsv')
fwrite(to_fix, 
       file = to_fix_out,
       sep = "\t")

cli::cli_alert_danger("Non-scientific species names written to: {to_fix_out}")
cli::cli_alert_warning("Please fix those. First few:")

head(to_fix) |> print()

clr = gs_df[!(not_sci)] |> 
  collapse::frename(species = binom, g = genus, s = species) |> 
  slt(-not_sci) 

gs_out = path(proj_dir, 
              'output/Genus_species.tsv')

fwrite(clr,
       file = gs_out,
       sep = "\t")

file_copy(gs_out,
          shiny_data_dir,
          overwrite = TRUE)

cli::cli_alert_success("Genus species table written to {gs_out}")


if (grepl("\\/", d$date[1])) {
  d = d |> 
    mtt(date = as.Date(date, format = "%m/%d/%Y"))
}

d = d |> 
  join(clr, 
       on = c("species" = "binom"),
       verbose = FALSE) |> 
    mtt(wk = lubridate::week(date),
        yr = lubridate::year(date)) 

uniq_flw = d$flower_phenophase |> 
  unlist2d() |> 
  stringr::str_split(", ") |> 
  unlist() |> 
  funique() |> 
  qDT() |> 
  setColnames("flw_pheno") |> 
  sbt(flw_pheno %!=% "") |> 
  mtt(cln_flw = janitor::make_clean_names(flw_pheno))

uniq_flw_out = path(proj_dir, 'output/uniq_flw.tsv')

fwrite(uniq_flw,
       file = uniq_flw_out,
       sep = "\t")

file_copy(uniq_flw_out,
          shiny_data_dir,
          overwrite = TRUE)

cli::cli_alert_success("Unique flower phenophases written to: {uniq_flw_out}")

get_flw_j = function(x) {
  uniq_flw$flw_pheno %iin% x
}

idx_df = d |> 
  mtt(flw_split = stringr::str_split(flower_phenophase, ", "),
      flw_i = map(flw_split, get_flw_j), 
      i = 1:fnrow(d))

nz_idx = slt(idx_df, i, flw_i)[,unlist(flw_i), by = i] |> 
  qM()
# ^ This is probably the most compact way to store this, expand to unique
# columns just for convenience...

flw_mat = matrix(FALSE, nrow = fnrow(d),
                 ncol = fnrow(uniq_flw))

flw_mat[nz_idx] = TRUE

colnames(flw_mat) = uniq_flw$cln_flw

d = cbind(d, qDT(flw_mat))

# add leaf phenophase columns ---------------------------------------------

uniq_lf = d$leaf_phenophase |> 
  unlist2d() |> 
  stringr::str_split(", ") |> 
  unlist() |> 
  funique() |> 
  qDT() |> 
  setColnames("lf_pheno") |> 
  sbt(lf_pheno %!=% "") |> 
  mtt(cln_lf = janitor::make_clean_names(lf_pheno))

uniq_lf_out = path(proj_dir, 'output/uniq_lf.tsv')

fwrite(uniq_lf,
       file = uniq_lf_out,
       sep = "\t")

file_copy(uniq_lf_out,
          shiny_data_dir,
          overwrite = TRUE)

cli::cli_alert_success("Unique flower phenophases written to: {uniq_lf_out}")

get_lf_j = function(x) {
  uniq_lf$lf_pheno %iin% x
}

idx_df = d |> 
  mtt(lf_split = stringr::str_split(leaf_phenophase, ", "),
      lf_i = map(lf_split, get_lf_j), 
      i = 1:fnrow(d))

nz_idx = slt(idx_df, i, lf_i)[,unlist(lf_i), by = i] |> 
  qM()
# ^ This is probably the most compact way to store this, expand to unique
# columns just for convenience...

lf_mat = matrix(FALSE, nrow = fnrow(d),
                 ncol = fnrow(uniq_lf))

lf_mat[nz_idx] = TRUE

colnames(lf_mat) = uniq_lf$cln_lf

d = cbind(d, qDT(lf_mat))

d = d |> slt(colnames(d) %!=% "na")
# write output ------------------------------------------------------------

out_f = api_out_file |> 
  path_real() |> 
  stringr::str_replace_all(c("api_out" = "adj"))

# d |> 
#   fwrite(file = out_f, sep = "\t")
save(d, file = out_f, compress = "gzip")
# gzipped RData = fast enough (~10ms slower) AND space efficient compared to tsv

cli::cli_alert_success("Shiny input formatted file written to: {out_f}")

file_copy(out_f,
          shiny_data_dir,
          overwrite = TRUE)

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

