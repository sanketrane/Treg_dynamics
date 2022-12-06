## clearing the environment
rm(list = ls())  
gc()    

library(rstan)
library(tidyverse)
####################################################################################
## testing models from stan
stanmodel_file <- file.path("stan_models", "MAP_Incumbent_simple.stan")
expose_stan_functions(stanmodel_file)

ts_seq <- round(sample(seq(70, 100, length.out = 10), 3, replace = FALSE)) %>% sort()
agebmtseq <- rep(49, length(ts_seq))
#init_cond <- c(0.2 * 9e4, 0.8*9e4, 0.2 * 5e4, 0.8*5e4, 0.2 * 8e5, 0.8*8e5, 0.2 * 4e5, 0.8*4e5,
#               0,0,0,0,0,0,0,0)
#init_cond <- c(exp(9.1), exp(10.6), exp(8.8), exp(10.4), exp(10.3), exp(10.7), exp(10.4), exp(11.4),
            #   0,0,0,0)
init_cond <- c(exp(11.1), exp(14.2),  exp(7.8), exp(11.8), 0,0)

params <- c(psi=0.011, alpha=0.95, delta_D=0.029, beta=0.018)

global_parms <- c(params,init_cond)
math_reduce(global_parms, local_params = c(0), x_r = ts_seq, x_i = agebmtseq)

#data_pred <- math_reduce(global_parms, local_params = c(0), x_r=solve_time, x_i = unique_times_counts$age.at.BMT)

# time sequence for predictions specific to age bins within the data
ts_pred1 <- 10^seq(log10(66), log10(450), length.out = 300)
ts_pred2 <- 10^seq(log10(91), log10(450), length.out = 300)
ts_pred3 <- 10^seq(log10(90), log10(450), length.out = 300)
ts_pred4 <- 10^seq(log10(174), log10(450), length.out = 300)
tb_pred1 <- rep(45, 300)
tb_pred2 <- rep(66, 300)
tb_pred3 <- rep(76, 300)
tb_pred4 <- rep(118, 300)

ode_pred1 <- math_reduce(global_parms, local_params = c(0), x_r = ts_pred1, x_i = tb_pred1)
ode_pred2 <- math_reduce(global_parms, local_params = c(0), x_r = ts_pred2, x_i = tb_pred2)
ode_pred3 <- math_reduce(global_parms, local_params = c(0), x_r = ts_pred3, x_i = tb_pred3)
ode_pred4 <- math_reduce(global_parms, local_params = c(0), x_r = ts_pred4, x_i = tb_pred4)

#data_pred_df <- data.frame("time_seq" = solve_time,
#                            "y_pred" = matrix(unlist(data_pred), nrow = length(solve_time), byrow = TRUE)) %>%
#  rename(counts_thy = y_pred.1,
#         counts_per = y_pred.2,
#         Nfd_thy = y_pred.3,
#         Nfd_per = y_pred.4,
#         donor_ki_thy = y_pred.7,
#         donor_ki_per = y_pred.9,
#         host_ki_thy = y_pred.5,
#         host_ki_per = y_pred.6)

stan_pred_df1 <- data.frame("time_seq" = ts_pred1,
                           "y_pred" = matrix(unlist(ode_pred1), nrow = length(ts_pred1), byrow = TRUE)) %>%
  rename(counts_thy = y_pred.1,
         counts_per = y_pred.2,
         Nfd_thy = y_pred.3,
         Nfd_per = y_pred.4) %>%
  bind_cols(ageBMT_bin = 'agebin1')

stan_pred_df2 <- data.frame("time_seq" = ts_pred2,
                            "y_pred" = matrix(unlist(ode_pred2), nrow = length(ts_pred2), byrow = TRUE)) %>%
  rename(counts_thy = y_pred.1,
         counts_per = y_pred.2,
         Nfd_thy = y_pred.3,
         Nfd_per = y_pred.4) %>%
  bind_cols(ageBMT_bin = 'agebin2')

stan_pred_df3 <- data.frame("time_seq" = ts_pred3,
                            "y_pred" = matrix(unlist(ode_pred3), nrow = length(ts_pred3), byrow = TRUE)) %>%
  rename(counts_thy = y_pred.1,
         counts_per = y_pred.2,
         Nfd_thy = y_pred.3,
         Nfd_per = y_pred.4) %>%
  bind_cols(ageBMT_bin = 'agebin3')

stan_pred_df4 <- data.frame("time_seq" = ts_pred4,
                            "y_pred" = matrix(unlist(ode_pred4), nrow = length(ts_pred4), byrow = TRUE)) %>%
  rename(counts_thy = y_pred.1,
         counts_per = y_pred.2,
         Nfd_thy = y_pred.3,
         Nfd_per = y_pred.4) %>%
  bind_cols(ageBMT_bin = 'agebin4')

Counts_pred <- rbind(stan_pred_df1, stan_pred_df2, stan_pred_df3, stan_pred_df4)

ggplot() +
  #geom_point(data = counts_data, aes(x = age.at.S1K, y = Thymus, color = ageBMT_bin)) +
  geom_line(data = Counts_pred, aes(x = time_seq, y = counts_thy, color = ageBMT_bin), linetype=2) +
  #geom_point(data = counts_data, aes(x = age.at.S1K, y = Periphery, color = ageBMT_bin)) +
  geom_line(data = Counts_pred, aes(x = time_seq, y = counts_per, color = ageBMT_bin)) +
  labs(title=paste('Total counts of thymic naive Tregs'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(60, 450) , trans="log10", breaks=c(10, 30, 100, 300))+
  scale_y_log10() 

ggplot() +
  geom_point(data = Nfd_data, aes(x = age.at.S1K, y = Thymus, col=ageBMT_bin)) +
  geom_line(data = Counts_pred, aes(x = time_seq, y = Nfd_thy, color = ageBMT_bin)) +
  labs(title=paste('Nfd thymic naive Tregs'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(60, 450), breaks=c(10, 30, 100, 300))

ggplot() +
  geom_point(data = Nfd_data, aes(x = age.at.S1K, y = Periphery, col=ageBMT_bin)) +
  geom_line(data = Counts_pred, aes(x = time_seq, y = Nfd_per, color = ageBMT_bin)) +
  labs(title=paste('Nfd peripheral naive Tregs'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(60, 450), breaks=c(10, 30, 100, 300))

ggplot() +
  geom_point(data = donorki_data, aes(x = age.at.S1K, y = Thymus), col=4) +
  geom_line(data = Counts_pred, aes(x = time_seq, y = donor_ki_thy), col=4) +
  geom_point(data = hostki_data, aes(x = age.at.S1K, y = Thymus), col=2) +
  geom_line(data = Counts_pred, aes(x = time_seq, y = host_ki_thy), col=2) +
  labs(title=paste('Nfd peripheral naive Tregs'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(60, 450) , trans="log10", breaks=c(10, 30, 100, 300)) +
  facet_wrap(.~ageBMT_bin)





library(deSolve)
library(tidyverse)
theta_spline <- function(Time, psi){psi * 10^6.4 * exp(-0.024 * (Time - 49))}
chi_spline <- function(Time){0.82 * (1 - exp(-0.063 * (Time - 10)))}
donor_eps_spline <- function(Time){exp(- 0.068 * (Time - 49)) + 0.38}

shm_chi <- function (t, state, parms) {
  eps_host = 0.32
  kloss = 1/3.5
  with(as.list(c(state, parms)),{
    dy1 <- theta_spline(t, psi) * (1- chi_spline(t - 40)) * eps_host + rho_D * (2 * y2 + y1) + beta * y3 - (kloss + alpha + delta_D) * y1
    dy2 <- theta_spline(t, psi) * (1- chi_spline(t - 40)) * (1 - eps_host) + kloss * y1 + beta * y4  - (rho_D + alpha + delta_D) * y2
    dy3 <- alpha * y1 + rho_D * (2 * y4 + y3) - (kloss + beta + delta_D) * y3
    dy4 <- alpha * y2 + kloss * y3 - (rho_D + beta + delta_D) * y4
    dy5 <- alpha * y7 + rho_I * (2 * y6 + y5) - (kloss + beta + rho_I) * y5
    dy6 <- alpha * y8 + kloss * y5 - (rho_I + beta + rho_I) * y6
    dy7 <- beta * y5 + rho_I * (2 * y8 + y7) - (kloss + alpha + rho_I) * y7
    dy8 <- beta * y6+ kloss * y7 - (rho_I + alpha + rho_I) * y8
    
    dy9 <- theta_spline(t, psi) * chi_spline(t - 40) * donor_eps_spline(t) + rho_D * (2 * y10 + y9) + beta * y11  - (kloss + alpha + delta_D) * y9
    dy10 <- theta_spline(t, psi) * Chi_spline(t - 40) * (1 - donor_eps_spline(t)) + kloss * y9 + beta * y12  - (rho_D + alpha + delta_D) * y10
    dy11 <- alpha * y9 + rho_D * (2 * y12 + y11) - (kloss + beta + delta_D) * y11
    dy12 <- alpha * y10 + kloss * y11 - (rho_D + beta + delta_D) * y12
    
    list(c(dy1, dy2, dy3, dy4, dy5, dy6, dy7, dy8, dy9, dy10, dy11, dy12))
    })
}

ts_seq <- round(sample(seq(70, 100, length.out = 10), 3, replace = FALSE)) %>% sort()
agebmtseq <- rep(49, length(ts_seq))
state <- c("y1"=exp(9.1), "y2"= exp(10.6), "y3"=exp(8.8), "y4"=exp(10.4), "y5"=exp(10.3), "y6"=exp(10.7),
               "y7"=exp(10.4), "y8"=exp(11.4), "y9"=0, "y10"=0, "y11"=0, "y12"=0)
params <- c(psi=0.0027, rho_D=0.00039, alpha=0.14, delta_D=0.0054, rho_I=0.037, beta=0.0023)

ode_solkve <- ode(state, ts_seq, shm_chi, params)



  

