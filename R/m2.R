options(digits = 3)
library(fastverse)
library(cmdstanr)
library(gptoolsStan)
library(ggplot2)

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

oak_obs = oaks |> 
  mtt(yr = year(date),
      wk = lubridate::week(date),
      of_or_pl = (flower_phenophase %like% "Pollen release") | 
        (flower_phenophase %like% "Open flowers")) |> 
  gby(yr, wk) |> 
  smr(prop = fmean(of_or_pl),
      k = fsum(of_or_pl),
      n = fnobs(of_or_pl)) |> 
  mtt(from_zf = FALSE)

zf_cnts = join(z_df, oak_obs, 
               on = c("yr", "wk"), 
               how = "anti") |> 
  rbind(oak_obs) |> 
  roworder(yr, wk) |> 
  qDT()

dodge = 1

yr_diff = fmax(zf_cnts$yr) - fmin(zf_cnts$yr)

plot_input = zf_cnts |> 
  mtt(d = as.IDate(paste0("2026-01-01")) + 7*wk - 3.5 + 
        dodge * (yr - 2023) - (yr_diff*dodge/2)) |> 
  mtt(yr = factor(yr)) 

p = plot_input |> 
  ggplot(aes(d, prop)) + 
  geom_point(aes(color = yr, group = yr)) + 
  scale_color_manual(values = pals::parula(8)[c(1,3,5,7)]) + 
  scale_x_date(labels = scales::label_date("%b"),
               breaks = as.Date(paste0("2026-", 1:12, "-01"))) + 
  theme_bw()

m1 = cmdstan_model("stan/m1.stan",
                   include_paths = gptools_include_path())

agg_cnts = zf_cnts |> 
  gby(wk) |> 
  smr(n_i = fsum(n),
      y = fsum(k)) |> 
  mtt(og = TRUE)

infl = 4

infill = data.table(wk = seq(1,52, by = 1/infl), 
                    n_i = 0, y = 0,
                    og = FALSE) |> 
  sbt(!(wk %% 1 == 0))

to_model = rbind(agg_cnts, infill) |> 
  roworder(wk)

to_model[, i := .I]

m2 = cmdstan_model("stan/m2.stan",
                   include_paths = gptools_include_path())

dl = list(n = 52,
          n_fill = fnrow(to_model),
          n_i = to_model$n_i,
          obs_idx = to_model[(og)] |> get_elem('i'),
          y = to_model$y,
          epsilon = 1e-8)

f2 = m2$sample(data = dl, parallel_chains = 4,
               refresh = 500)

p_dt = f2$summary('proba') |> 
  qDT() |> 
  mtt(d = as.Date("2026-01-01") + 1:(fnrow(to_model)) * (7/4) - 3.5/4 ) |> 
  cbind(to_model |> slt(wk))

zf_cnts[k>0] |> roworder(wk) |> sbt(prop > .5)

zf_cnts |> sbt(wk %iin% 20:21)

obs_prop = to_model[(og)] |> 
  mtt(p = y/n_i) |> 
  mtt(d = as.Date("2026-01-01") + 1:(52) * (7) - 3.5 ) 

plot_input |> 
  ggplot(aes(d)) + 
  geom_point(aes(y = prop, color = yr, group = yr)) + 
  geom_line(data = obs_prop,
             aes(y = p),
            col = "blue") +
  scale_color_manual(values = pals::parula(8)[c(1,3,5,7)]) + 
  scale_x_date(labels = scales::label_date("%b"),
               breaks = as.Date(paste0("2026-", 1:12, "-01"))) + 
  geom_vline(xintercept = as.Date("2026-05-17")) + # at wk=20
  theme_bw() + 
  geom_ribbon(data = p_dt,
              aes(ymin = q5,
                  ymax = q95),
              alpha = .2) + 
  geom_line(data = p_dt,
            aes(y = mean),
            color = "red") 
