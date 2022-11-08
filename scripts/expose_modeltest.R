## clearing the environment
rm(list = ls())  
gc()    

library(rstan)
library(tidyverse)
####################################################################################
## testing models from stan
stanmodel_file <- file.path("stan_models", "rtem_test.stan")
expose_stan_functions(stanmodel_file)

ts_seq <- round(sample(seq(70, 100, length.out = 10), 5, replace = FALSE)) %>% sort()
agebmtseq <- rep(49, length(ts_seq))
init_cond <- c(0.2 * 9e4, 0.8*9e4, 0.2 * 5e4, 0.8*5e4, 0.2 * 8e5, 0.8*8e5, 0.2 * 4e5, 0.8*4e5,
               0,0,0,0,0,0,0,0)
params <- c(0.1, 0.08, 0.03, 0.05, 0.12, 0.005, 0.02, 0.01)
ode_sol <- solve_ode_chi(ts_seq, agebmtseq, init_cond, params)

stan_pred_df <- data.frame("time_seq" = ts_seq,
                           "y_pred" = matrix(unlist(ode_sol), nrow = length(ts_seq), byrow = TRUE)) 
