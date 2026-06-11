library(stringr)
library(fastverse)
library(lubridate)
library(ggplot2)
library(fs)

d = fread("data/anecdata_export_EwA_Pheno_Lite_2026-04-27T02-24-48-852Z.csv") |> 
  janitor::clean_names()

# I don't have family-level annotations. Let's just do oaks.
oaks = d |> 
  sbt(species_2 %like% "[Qq]uercus") |> 
  slt(date, species = species_2, lat, long = lng,
      leaf_phenophase, flower_phenophase)

# df of 0's for months with no observations
z_df = expand.grid(yr = 2023:2026, 
                   mo = 1:12) |> 
  mtt(prop = 0, n = 0, k = 0) |> 
  sbt(!(yr == 2026 & mo > 4))

oak_obs = oaks |> 
  mtt(yr = year(date),
      mo = month(date),
      of_or_pl = (flower_phenophase %like% "Pollen release") | 
                 (flower_phenophase %like% "Open flowers")) |> 
  gby(yr, mo) |> 
  smr(prop = fmean(of_or_pl),
      k = fsum(of_or_pl),
      n = fnobs(of_or_pl)) 

# Complete observations with rows for 0 entries
oak_cpl = join(z_df, oak_obs, 
               on = c("yr", "mo"), 
               how = "anti") |> 
  rbind(oak_obs) |> 
  roworder(yr, mo) |> 
  mtt(q = 1 - prop,
      se = fifelse(n > 0, sqrt(prop*q / sqrt(n)), 0),
      a = fifelse(n > 0,   k+1, 0),
      b = fifelse(n > 0, n-k+1, 0),
      post_lo = fifelse(n > 3, qbeta(.1, a, b), 0),
      post_hi = fifelse(n > 3, qbeta(.9, a, b), 0))

dodge = 1.5

plot_input = oak_cpl |> 
  mtt(d = as.IDate(paste0("2026-", mo, "-15")) + # change everything to 2026, then hide the year
           dodge * (yr - 2023) - (3 * dodge/2)) |> 
  mtt(yr = factor(yr)) 
  
plot_input |> 
  ggplot(aes(d, prop)) + 
  geom_line(aes(color = yr, 
                group = yr),
            alpha = .5) + 
  geom_point(data = plot_input |> sbt(n %==% 0),
             aes(color = yr)) + 
  geom_segment(
    data = plot_input |> sbt(n>0),
    aes(color = yr,
        y = post_lo, yend = post_hi,
        xend = d)
    ) + 
  geom_point(data = plot_input |> sbt(n>0),
             aes(color = yr)) + 
  ylim(c(NA,1)) + 
  scale_x_date(labels = scales::label_date("%b"),
               breaks = as.Date(paste0("2026-", 1:12, "-15"))) + 
  theme_bw() + 
  labs(title = "Quercus species with open flowers or pollen release",
       subtitle = "80% posterior interval with Beta(1,1) prior",
       caption = "Data: EwA Pheno Lite",
       x = NULL,
       y = "proportion",
       color = NULL) + 
  scale_color_manual(values = pals::parula(8)[c(1,3,5,7)]) + 
  theme(text = element_text(family = "Arial"),
        panel.grid.minor = element_blank())

ggsave(path("output", "quercus_openfl_pollen.png"),
       w = 6, h = 4)
