library(rstan)

# m = stan_model(file = "binom_gp.stan", save_dso = FALSE, auto_write = FALSE)

m = readRDS('binom_gp.rds')
