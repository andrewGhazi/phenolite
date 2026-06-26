library(rstan)

setwd("~/projects/phenolite/shiny/phenoshiny/")

m = stan_model(file = "binom_gp.stan")

saveRDS(m, 'binom_gp.rds')


