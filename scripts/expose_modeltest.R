## clearing the environment
rm(list = ls())  
gc()    

library(rstan)
library(tidyverse)
####################################################################################
## testing models from stan
stanmodel_file <- file.path("stan_models", "MAP_rtem_naiTreg.stan")
expose_stan_functions(stanmodel_file)

ts_seq <- round(sample(seq(70, 100, length.out = 10), 3, replace = FALSE)) %>% sort()
agebmtseq <- rep(49, length(ts_seq))
init_cond <- c(0.2 * 9e4, 0.8*9e4, 0.2 * 5e4, 0.8*5e4, 0.2 * 8e5, 0.8*8e5, 0.2 * 4e5, 0.8*4e5,
               0,0,0,0,0,0,0,0)
#init_cond <- c(exp(8.5), exp(10), exp(8.75), exp(12.4), exp(10.5), exp(12), exp(9.9), exp(11.14),
#               0,0,0,0)
params <- c(0.002, 0.0003, 0.15, 0.003, 0.03, 0.004, 0.03, 0.01)

global_parms <- c(params,init_cond)
ode_sol <- solve_ode_chi(ts_seq, agebmtseq, init_cond, params)
math_reduce(global_parms, local_params = c(0), x_r = 90, x_i = 80)

data_pred <- math_reduce(global_parms, local_params = c(0), x_r=solve_time, x_i = unique_times_counts$age.at.BMT)

# time sequence for predictions specific to age bins within the data
ts_pred1 <- 10^seq(log10(66), log10(200), length.out = 300)
ts_pred2 <- 10^seq(log10(91), log10(330), length.out = 300)
ts_pred3 <- 10^seq(log10(90), log10(350), length.out = 300)
ts_pred4 <- 10^seq(log10(174), log10(450), length.out = 300)
tb_pred1 <- rep(45, 300)
tb_pred2 <- rep(66, 300)
tb_pred3 <- rep(76, 300)
tb_pred4 <- rep(118, 300)

ode_pred1 <- solve_ode_chi(ts_pred1, tb_pred1, init_cond, params)
ode_pred2 <- solve_ode_chi(ts_pred2, tb_pred2, init_cond, params)
ode_pred3 <- solve_ode_chi(ts_pred3, tb_pred3, init_cond, params)
ode_pred4 <- solve_ode_chi(ts_pred4, tb_pred4, init_cond, params)

stan_pred_df1 <- data.frame("time_seq" = ts_pred1,
                           "y_pred" = matrix(unlist(ode_pred1), nrow = length(ts_pred1), byrow = TRUE)) %>%
  mutate(counts_thy = y_pred.1 + y_pred.2 + y_pred.7 + y_pred.8 + y_pred.9 + y_pred.10 + y_pred.15 + y_pred.16,
         counts_per = y_pred.3 + y_pred.4 + y_pred.5 + y_pred.6 + y_pred.11 + y_pred.12 + y_pred.13 + y_pred.14,
         total_counts = counts_thy + counts_per,
         ageBMT_bin = 'agebin1') %>%
  select(time_seq, ageBMT_bin, contains('counts'))


stan_pred_df2 <- data.frame("time_seq" = ts_pred2,
                            "y_pred" = matrix(unlist(ode_pred2), nrow = length(ts_pred2), byrow = TRUE)) %>%
  mutate(counts_thy = y_pred.1 + y_pred.2 + y_pred.7 + y_pred.8 + y_pred.9 + y_pred.10 + y_pred.15 + y_pred.16,
         counts_per = y_pred.3 + y_pred.4 + y_pred.5 + y_pred.6 + y_pred.11 + y_pred.12 + y_pred.13 + y_pred.14,
         total_counts = counts_thy + counts_per,
         ageBMT_bin = 'agebin2') %>%
  select(time_seq, ageBMT_bin, contains('counts'))


stan_pred_df3 <- data.frame("time_seq" = ts_pred3,
                            "y_pred" = matrix(unlist(ode_pred3), nrow = length(ts_pred3), byrow = TRUE)) %>%
  mutate(counts_thy = y_pred.1 + y_pred.2 + y_pred.7 + y_pred.8 + y_pred.9 + y_pred.10 + y_pred.15 + y_pred.16,
         counts_per = y_pred.3 + y_pred.4 + y_pred.5 + y_pred.6 + y_pred.11 + y_pred.12 + y_pred.13 + y_pred.14,
         total_counts = counts_thy + counts_per,
         ageBMT_bin = 'agebin3') %>%
  select(time_seq, ageBMT_bin, contains('counts'))


stan_pred_df4 <- data.frame("time_seq" = ts_pred4,
                            "y_pred" = matrix(unlist(ode_pred4), nrow = length(ts_pred4), byrow = TRUE)) %>%
  mutate(counts_thy = y_pred.1 + y_pred.2 + y_pred.7 + y_pred.8 + y_pred.9 + y_pred.10 + y_pred.15 + y_pred.16,
         counts_per = y_pred.3 + y_pred.4 + y_pred.5 + y_pred.6 + y_pred.11 + y_pred.12 + y_pred.13 + y_pred.14,
         total_counts = counts_thy + counts_per,
         ageBMT_bin = 'agebin4') %>%
  select(time_seq, ageBMT_bin, contains('counts'))


Counts_pred <- rbind(stan_pred_df1, stan_pred_df2, stan_pred_df3, stan_pred_df4)

init_line <- sapply(ts_pred1, solve_init, init_cond=init_cond, parms=params)
init_df <- data.frame("time_seq" = ts_pred1,
                      "y_pred" = matrix(unlist(init_line), nrow = length(ts_pred4), byrow = TRUE)) %>%
  mutate(counts_thy = y_pred.1 + y_pred.2 + y_pred.7 + y_pred.8 + y_pred.9 + y_pred.10 + y_pred.15 + y_pred.16,
         counts_per = y_pred.3 + y_pred.4 + y_pred.5 + y_pred.6 + y_pred.11 + y_pred.12 + y_pred.13 + y_pred.14,
         total_counts = counts_thy + counts_per,
         ageBMT_bin = 'agebin4') %>%
  select(time_seq, ageBMT_bin, contains('counts'))

ggplot() +
  geom_line(data = init_df, aes(x = time_seq, y = total_counts), col='darkred', size=3) +
  geom_line(data = Counts_pred, aes(x = time_seq, y = total_counts, color = ageBMT_bin), size=1.2) +
  labs(title=paste('Total counts of thymic naive Tregs'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(60, 450) , trans="log10", breaks=c(10, 30, 100, 300))+
  scale_y_continuous() 











  

