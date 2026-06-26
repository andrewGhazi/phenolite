# This version uses rstan to update the lines

library(shiny)
library(ggplot2)
library(data.table)
library(collapse)
library(ggtext)
library(forcats)

source("load_rstan.R")

to_fct = \(x, genus) {
  x |> 
    fctr() |> 
    fct_relabel(\(x) paste0("*", x, "*")) |> 
    fct_lump_min(5, other_level = paste0("other *", genus, "*")) |> 
    fct_infreq() |> 
    fct_rev() 
}

print("reading...")
# d = fread("data/anecdata_export_EwA_Pheno_Lite_2026-05-28T20-24-36-075Z.csv") |> 
#   janitor::clean_names()

d = fread("data/adj_2026-05-28T20-24-36-075Z.csv")

uniq_lf = fread("data/uniq_lf.tsv") 

uniq_flw = fread('data/uniq_flw.tsv')

latest_wk = d |> 
  slt(wk, yr) |> 
  funique() |> 
  sbt(yr %==% fmax(yr)) |> 
  sbt(wk %==% fmax(wk))

genera = d |> 
  slt(genus) |> 
  fcount() |> 
  roworder(-N) |> 
  na_omit() |> 
  sbt(N >= 30) # can select from genera with >= 30 observations

pd = d |> 
  sbt(species_2 %==% "Quercus rubra") |> 
  slt(sp = species_2, date, open_flowers, fruits, leaf_buds_swelling) |> 
  pivot(ids = 1:2) |> 
  mtt(yr = lubridate::year(date)) 

lubridate::year(pd$date) <- 2026

pd |> 
  ggplot(aes(date, variable)) + 
  geom_jitter(aes(alpha = value),
              width = 0,
              height = .2,
              size = .7,
              pch = 16) + 
  facet_wrap(ncol = 1, 
             vars(yr)) + 
  labs(y = NULL,
       title = "*Quercus rubra*",
       x = NULL,
       alpha = "observed") + 
  theme_bw() + 
  theme(plot.title = element_markdown(),
        panel.grid.major.y = element_blank(),
        panel.spacing = unit(1, "pt"),
        strip.text = element_text(margin = margin(0,0,0,0))) + 
  scale_x_date(labels = scales::label_date("%b"),
               breaks = as.Date(paste0("2026-", 1:12, "-15"))) 

