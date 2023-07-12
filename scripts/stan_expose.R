rm(list = ls()); gc()

library(rstan)
library(tidyverse)

## model
rstan::expose_stan_functions("stan_models/MAP_Incumbent.stan")
params <- c(1.363198e+05, 0.01, 0.04, 0.02, 0.002)
par_inc <- c(0.3, 0.01, 0.05, 0.02, 10, 11, 9, 10)
theta <- c(0)
x_i <- c(60)
x_r <- c(100)
math_reduce(par_inc, theta, x_r, x_i)

## ts
ts_pred1 <- 10^seq(log10(66), log10(450), length.out = 300)
ts_pred2 <- 10^seq(log10(91), log10(450), length.out = 300)
ts_pred3 <- 10^seq(log10(90), log10(450), length.out = 300)
ts_pred4 <- 10^seq(log10(174), log10(450), length.out = 300)
tb_pred1 <- rep(45, 300)
tb_pred2 <- rep(66, 300)
tb_pred3 <- rep(76, 300)
tb_pred4 <- rep(118, 300)


theta_spline(40, params)
g_age(0, params)
ageseq <- seq(0, 100, 0.1)
gvec <- sapply(ageseq, g_age, parms=params)
ggplot()+
  geom_line(aes(x=ageseq, y=gvec))

logit_transf <- function(x){log(x/(1-x))}
asinsq_transf <- function(x){asin(sqrt(x))}


chi_vec <- sapply(solve_time - ageAtBMT, Chi_spline)
total_counts <- N_pooled_time(solve_time, ageAtBMT, params)
donor_counts <- N_donor_time(solve_time, ageAtBMT, params)
Nfd_pred1 <- donor_counts/(total_counts * chi_vec)

asinsq_transf <- function(x){asin(sqrt(x))}

ggplot()+
  geom_point(aes(x = solve_time, y=Nfd_pred1)) + ylim(0,1)


host_counts_mean = N_host_time(solve_time, ageAtBMT, params);
host_ki_counts = U_host_time(solve_time, ageAtBMT, params);
host_ki_mean = host_ki_counts/host_counts_mean;
donor_ki_counts = U_donor_time(solve_time, ageAtBMT, params);
donor_ki_mean = donor_ki_counts/donor_counts;


ggplot()+
  geom_point(aes(x = solve_time, y=host_ki_mean)) +
  geom_point(aes(x = solve_time, y=donor_ki_mean), col=2) +
  ylim(0,1)


asinsq_transf(donor_ki_mean)
asinsq_transf(host_ki_mean)



