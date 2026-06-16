library(fastverse)
library(Rcpp)
library(RcppArmadillo)

get_L = \(l, wk_sigma) {
  cov_mat = matrix(0, nr = 52, nc = 52)
  # cov_mat = diag(52)
  
  j = 1
  
  for (i in 2:27) {
    cov_mat[i,j] = -(i-1)^2 / (2*l^2)
  }
  
  cov_mat[27:52,j] = cov_mat[27:2,1]
  
  cv = cov_mat[,j]
  
  for (j in 2:52) {
    cv = cv[c(52, 1:51)]
    cov_mat[,j] = cv
  }
  
  cov_mat = wk_sigma^2 * exp(cov_mat)
  
  L = chol(cov_mat) |> t()
}

sourceCpp(file = "~/projects/phenolite/src/temporal_L.cpp",
          cacheDir = "~/projects/phenolite/src/cch")

ll = function(pvec, N, k, n, wk_i) {
  l = pvec[1]
  # print(head(pvec))
  
  wk_sigma = pvec[2]
  
  int = pvec[3]
  
  n_hp = 3
  # print(l)
  # print(wk_sigma)
  
  wkv = pvec[1:52 + n_hp]
  
  yrwkv = tail(pvec, N)
  
  L = getL(l, wk_sigma)
  
  wk_ll = dnorm(solve(L, matrix(wkv, nc = 1)), log = TRUE) |> sum()
  # ^ check this...
  
  lp = dexp(wk_sigma, rate = 10, log = TRUE) +
    dexp(l, rate = .5, log = TRUE) + 
    sum(dnorm(yrwkv, sd = .5, log = TRUE)) + 
    wk_ll + 
    dnorm(int, log = TRUE)
  
  # wk_sigma - prior mean of 0.1. Don't let the wk effects eat everything
  # l prior mean: 2 week length scale
  # std normal prior on yrwkv effects. Don't understand why these don't eat everything though...
  
  ll = dbinom(k, 
              size = n, 
              prob = plogis(wkv[wk_i] + yrwkv + int), 
              log = TRUE) |> 
    sum()
  
  lp + ll
  
}
# 
# pvec = c(1,1, 0, alloc(0, 52+178))
# 
# ll(pvec, zfc)
# 
# pvec = c(2,1, 0, alloc(0, 52+178))
# 
# ll(pvec, zfc)
# 
# (opt_res = optim(pvec, ll, zfc = zf_cnts,
#       control = list(fnscale = -1,
#                      maxit = 1000),
#       lower = c(.3,.1, alloc(-Inf, 52+178)),
#       method = "L-BFGS-B"))
# # TODO, make this optim faster with Stan or something...
# 
# opt_res$par |> plot()
# 
# eff = opt_res$par |> tail(-3)
# wkv_fit = eff |> head(52)
# yrwkv_fit = eff |> tail(-52)
# 
# zf_cnts = zf_cnts |> 
#   mtt(fitted = plogis(wkv_fit[zf_cnts$wk] + yrwkv_fit + opt_res$par[3]))


