library(stringr)
library(fastverse)
library(lubridate)
library(ggplot2)
library(fs)

d = fread("data/anecdata_export_EwA_Pheno_Lite_2026-05-28T20-24-36-075Z.csv") |> 
  janitor::clean_names()

d |> 
  mtt(wk = lubridate::week(date), 
      yr = lubridate::year(date)) |> 
  sbt(leaf_phenophase %like% "Breaking" & 
        species_2 %like% "[Qq]uercus") |> 
  fcount(wk, yr) |>
  ggplot(aes(wk, N)) + 
  geom_point(aes(color = factor(yr))) + 
  geom_smooth(aes(color = factor(yr)),
              se = FALSE) + 
  labs(title = "Breaking leaf buds, weekly") + 
  xlim(c(0, 52))

d |> 
  sbt(species_2 %like% "[Qq]uercus") |> 
  mtt(dy = lubridate::yday(date)) |> 
  ggplot(aes(date, species_2)) + 
  geom_point()


d |> 
  sbt(species_2 %like% "[Qq]uercus" & 
        leaf_phenophase %like% "Breaking") |> 
  mtt(dy = lubridate::yday(date)) |> 
  ggplot(aes(dy, species_2)) + 
  geom_jitter(width = 0, height = .1, alpha = .1, pch = 16) + 
  xlim(c(0,365)) + 
  labs(x = "day of the year")
