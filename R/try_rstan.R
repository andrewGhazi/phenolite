library(rstan)

# The model is so fast that maybe cmdstanr is comparatively slow just because it
# has to interact with disk. Try rstan.

fit1 = stan(file = "~/projects/phenolite/stan/binom_gp.stan",
            data = dl)

m = stan_model(file = "~/projects/phenolite/stan/binom_gp.stan")

init_fun = \() {
  pvec = c(1,
           fsd(agg_cnts$qinit),
           fmean(agg_cnts$qinit), 
           W(agg_cnts$qinit),
           alloc(0, fnrow(zf_cnts)))
  
  list(ell = pvec[1],
       wk_sigma = pvec[2],
       intercept = pvec[3],
       wkv = pvec[1:52 + 3],
       yrwkv = alloc(0, fnrow(zf_cnts)))
}

optimizing(m, data = dl, init = init_fun)

bench::mark({
  optimizing(m, data = dl, init = init_fun)
})

# Yep, about 6x faster
