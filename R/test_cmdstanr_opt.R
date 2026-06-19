library(fastverse)
library(cmdstanr)

# Well writing the function myself is pretty hard, let's try it with stan for
# the sake of comparison.

dl = list(N = fnrow(zf_cnts),
          k = zf_cnts$k,
          n = zf_cnts$n,
          wk_i = zf_cnts$wk)

m = cmdstan_model("~/projects/phenolite/stan/binom_gp.stan")

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

f = m$optimize(dl, 
               show_messages = FALSE, 
               init = init_fun)


br = bench::mark(e1 = {
  f = m$optimize(dl, 
                 show_messages = FALSE, 
                 init = init_fun)
}, e2 = {
  f = m$optimize(dl, 
                 show_messages = FALSE)
}, e3 = {
  f = m$optimize(dl, 
                 show_messages = FALSE, 
                 init = init_fun, 
                 algorithm = 'lbfgs',
                 tol_obj = 1e-2)
}, #e4 = {
#   
#   f = m$optimize(dl, 
#                  show_messages = FALSE, 
#                  init = init_fun, 
#                  algorithm = 'newton')
# }, 
check = FALSE)

# about 20x faster

print(br)
plot(br)

# ests = f$summary(variables = c("ell", "wk_sigma", "intercept", "wkv", "yrwkv"))

f$draws(format = 'data.frame', 
        variables =  c("ell", "wk_sigma", "intercept", "wkv", "yrwkv")) |> unlist() |> head(-3)

# 136ms...

plot(ests$estimate)
