options(digits = 3)
library(fastverse)
library(cmdstanr)
library(gptoolsStan)
library(ggplot2)
library(forcats)

d = fread("data/anecdata_export_EwA_Pheno_Lite_2026-05-28T20-24-36-075Z.csv") |> 
  janitor::clean_names()

oaks = d |> 
  mtt(wk = lubridate::week(date), 
      yr = lubridate::year(date)) |> 
  sbt(species_2 %like% "[Qq]uercus") 

z_df = expand.grid(yr = 2023:2026, 
                   wk = 1:52) |> 
  mtt(prop = 0, n = 0, k = 0) |> 
  sbt(!(yr == 2026 & wk > 19)) |> 
  mtt(from_zf = TRUE)

to_fct = \(x) {
  x |> 
    fctr() |> 
    fct_lump_min(5) |> 
    fct_infreq() |> 
    fct_rev() 
}

dodge = 0

oak_obs = oaks |> 
  mtt(yr = year(date),
      wk = lubridate::week(date),
      of_or_pl = (flower_phenophase %like% "Pollen release") | 
                 (flower_phenophase %like% "Open flowers")) |> 
  slt(date, species = species_2, lat, lng, leaf_pheno = leaf_phenophase,
      flw_pheno = flower_phenophase,
      wk:of_or_pl) |> 
  mtt(species = to_fct(species),
      pres = c("absent", "present")[of_or_pl + 1],
      d_hide  = as.IDate(paste0("2026-01-01")) + 7*wk - 3.5 + 
           dodge * (yr - 2023) - (yr_diff*dodge/2),
      yd = date) |> 
  roworder(of_or_pl) 

lubridate::year(oak_obs$yd) <- 2026

oak_obs |> 
  ggplot(aes(d_hide, species)) + 
  geom_point(pch = 15, 
             aes(color = pres),
             position = position_jitter(width=.5, height = .2),
             size = .6) + 
  scale_color_manual(values = c("grey", "black")) + 
  labs(color = "Pollen release | open flowers",
       x = NULL,
       y = NULL) + 
  scale_x_date(labels = scales::label_date("%b"),
               breaks = as.Date(paste0("2026-", 1:12, "-01"))) + 
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank()) 
  
oak_obs |> 
  ggplot(aes(yd, species)) + 
  geom_point(pch = 15, 
             aes(color = pres),
             position = position_jitter(width=0, height = .2),
             size = .6) + 
  scale_color_manual(values = c("grey", "black")) + 
  labs(color = "Pollen release | open flowers",
       x = NULL,
       y = NULL) + 
  scale_x_date(labels = scales::label_date("%b"),
               breaks = as.Date(paste0("2026-", 1:12, "-01"))) + 
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank()) 
  
