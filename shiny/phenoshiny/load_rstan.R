library(rstan)

m = stan_model(file = "binom_gp.stan", save_dso = FALSE, auto_write = FALSE)

# Temporarily edited makeconf_path() to match shinyapps: https://discourse.mc-stan.org/t/compile-stan-model-in-shiny-app/10022/18
# Run this in a vanilla R session (in the shiny app's directory), not RStudio:
library(fs)

setwd("~/projects/phenolite/shiny/phenoshiny")

m = stan_model(file = "binom_gp.stan", auto_write = TRUE)

lp = get("dso_last_path", m@dso@.CXXDSOMISC)

dlp = lp |> dirname()

shiny_dir = "~/projects/phenolite/shiny/phenoshiny/stan/binom_gp/"

shiny_pth = "stan/binom_gp/"

file_copy(dir_ls(dlp), shiny_dir,
          overwrite = TRUE)

assign("dso_last_path",
       path(shiny_pth, basename(lp)),
       envir = m@dso@.CXXDSOMISC)

mod_envir = get("module", m@dso@.CXXDSOMISC)

pN = get("packageName", mod_envir)

pN[['path']] <- path(shiny_pth, basename(pN[['path']]))

assign('packageName', pN,
       envir = mod_envir)

assign("module", mod_envir,
       envir = m@dso@.CXXDSOMISC)

print(m@dso)

saveRDS(m, 'binom_gp.rds')
# Let's hope this precompiled model runs on shinyapps.io servers

m = readRDS('binom_gp.rds')
print(get("dso_last_path", m@dso@.CXXDSOMISC))
