options(digits = 3)
library(fastverse)
library(cmdstanr)
library(gptoolsStan)
library(purrr)
library(stringr)

d = fread("data/anecdata_export_EwA_Pheno_Lite_2026-05-28T20-24-36-075Z.csv") |> 
  janitor::clean_names() |> 
  slt(species = species_2) |> 
  funique() |> 
  mtt(sp_pcs = strsplit(species, " "),
      n_pc = lengths(sp_pcs))

d |> 
  sbt(n_pc != 2)

gs_df = d |> 
  sbt(n_pc %==% 2) |> 
  mtt(g = map_chr(sp_pcs, 1) |> str_to_title(),
      s = map_chr(sp_pcs, 2) |> tolower()) |> 
  slt(species, g,s) |> 
  funique() |> 
  roworder(g)

gs_df |> print(nrow = 352)  

gs_df$not_sci = FALSE

cmn_i = c(28:29, 55, 65, 67, 154, 155, 167, 209, 293, 309, 348:349)

gs_df$species[cmn_i] |> dput()

gs_df$not_sci[cmn_i] = TRUE

to_fix = gs_df[(not_sci)]

clr = gs_df[!(not_sci)] |> 
  collapse::frename(species = binom, g = genus, s = species) |> 
  slt(-not_sci) 

clr |> 
  fwrite(file = 'output/Genus_species.tsv')

clr |> fcount(genus) |> roworder(N)

fread("data/anecdata_export_EwA_Pheno_Lite_2026-05-28T20-24-36-075Z.csv") |> 
  janitor::clean_names() |> 
  join(clr, on = c("species_2" = "binom")) |> 
  fcount(genus) |> 
  roworder(N) |> sbt(N > 52)

to_fix

d$n_pc |> table()

d |> sbt(n_pc == 1)

d |> sbt(n_pc == 3)

d |> sbt(n_pc == 4)
