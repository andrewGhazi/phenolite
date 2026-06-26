library(rstan)

# Try to adapt rstan::stan_model
# https://github.com/stan-dev/rstan/blob/3ff4a321a98f80c8365ec01e5524116a6c183668/rstan/rstan/R/rstan.R#L19

stanc_ret = stanc(file = "binom_gp.stan")


model_cppname <- stanc_ret$model_cppname
model_name <- stanc_ret$model_name
model_code <- stanc_ret$model_code
model_cppcode <- stanc_ret$cppcode
model_cppcode <- paste("#ifndef MODELS_HPP",
                       "#define MODELS_HPP",
                       "#define STAN__SERVICES__COMMAND_HPP",
                       "#include <rstan/rstaninc.hpp>",
                       model_cppcode,
                       "#endif",
                       sep = '\n')


inc <- paste("#include <Rcpp.h>\n",
             "using namespace Rcpp;\n",
             model_cppcode ,
             rstan:::get_Rcpp_module_def_code(model_cppname),
             sep = '')

fx = cxxfunction(sig = sig, body = body, plugin = plugin, includes = includes,
                 settings = settings, ..., verbose = verbose)