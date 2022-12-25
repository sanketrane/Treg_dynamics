## clearing the environment
rm(list = ls())  
gc()    

library(rstan)
library(deSolve)
library(tidyverse)
####################################################################################
## testing models from stan
stanmodel_file <- file.path("stan_models", "MAP_Incumbent.stan")
expose_stan_functions(stanmodel_file)

ts_seq <- round(sample(seq(70, 100, length.out = 10), 3, replace = FALSE)) %>% sort()
agebmtseq <- rep(49, length(ts_seq))
#init_cond <- c(0.2 * 9e4, 0.8*9e4, 0.2 * 5e4, 0.8*5e4, 0.2 * 8e5, 0.8*8e5, 0.2 * 4e5, 0.8*4e5,
#               0,0,0,0,0,0,0,0)
init_cond <- c("y1"=exp(9.70302424), "y2"= exp(12.60218721), "y3"=exp(9.64950000), "y4"=exp(12.45474339), "y5"=exp(10.43676166),
               "y6"=exp(10.86395695), "y7"=exp(10.03221739), "y8"=exp(10.52207092), "y9"=0, "y10"=0, "y11"=0, "y12"=0)
#init_cond <- c("y1"=exp(9.6), "y2"= exp(14.3), "y3"=exp(9.1), "y4"=exp(11.9), "y5"=0, "y6"=0)

params <- c(psi=0.03021270, rho_D=0.01920767, alpha=0.07980150, delta_D=0.05992757, rho_I=0.02313391, beta=0.03173747)
#params <- c(psi=0.011, alpha=0.83, delta_D=0.027, delta_I=0.001, beta=0.015)

#data_pred <- math_reduce(global_parms, local_params = c(0), x_r=solve_time, x_i = unique_times_counts$age.at.BMT)

thy0 <- init_cond[1] + init_cond[2] + init_cond[7] + init_cond[8] + init_cond[9] + init_cond[10]
per0 <- init_cond[3] + init_cond[4] + init_cond[5] + init_cond[6] + init_cond[11] + init_cond[12]


# time sequence for predictions specific to age bins within the data
ts_pred1 <- 10^seq(log10(66), log10(450), length.out = 300)
ts_pred2 <- 10^seq(log10(91), log10(450), length.out = 300)
ts_pred3 <- 10^seq(log10(90), log10(450), length.out = 300)
ts_pred4 <- 10^seq(log10(174), log10(450), length.out = 300)
tb_pred1 <- rep(45, 300)
tb_pred2 <- rep(66, 300)
tb_pred3 <- rep(76, 300)
tb_pred4 <- rep(118, 300)

ode_pred1 <- solve_ode_chi(ts_pred1, tb_pred1, init_cond, params)
ode_pred2 <- solve_ode_chi(ts_pred2, tb_pred2, init_cond, params)
ode_pred3 <- solve_ode_chi(ts_pred3, tb_pred3, init_cond, params)
ode_pred4 <- solve_ode_chi(ts_pred4, tb_pred4, init_cond, params)

chivec1 <- sapply(ts_pred1-tb_pred1, Chi_spline)
chivec2 <- sapply(ts_pred2-tb_pred2, Chi_spline)
chivec3 <- sapply(ts_pred3-tb_pred3, Chi_spline)
chivec4 <- sapply(ts_pred4-tb_pred4, Chi_spline)

stan_pred_df1 <- data.frame("time_seq" = ts_pred1,
                           "y_pred" = matrix(unlist(ode_pred1), nrow = length(ts_pred1), byrow = TRUE)) %>%
  mutate(counts_thy = y_pred.1 + y_pred.2 + y_pred.7 + y_pred.8 + y_pred.9 + y_pred.10,
         counts_per = y_pred.3 + y_pred.4 + y_pred.5 + y_pred.6 + y_pred.11 + y_pred.12,
         total_counts = counts_thy + counts_per,
         Nfd_thy = (y_pred.9 + y_pred.10)/(counts_thy * chivec1),
         Nfd_per = (y_pred.11 + y_pred.12)/(counts_per * chivec1),
         donor_ki_thy = (y_pred.9)/(y_pred.9 + y_pred.10),
         donor_ki_per = (y_pred.11)/(y_pred.11 + y_pred.12),
         host_ki_thy = (y_pred.1 + y_pred.7)/(y_pred.1 + y_pred.2 + y_pred.7 + y_pred.8),
         host_ki_per = (y_pred.3 + y_pred.5)/(y_pred.3 + y_pred.4 + y_pred.5 + y_pred.6),
         ageBMT_bin = 'agebin1') %>%
  select(time_seq, ageBMT_bin, contains('thy'), contains('per'))


stan_pred_df2 <- data.frame("time_seq" = ts_pred2,
                            "y_pred" = matrix(unlist(ode_pred2), nrow = length(ts_pred2), byrow = TRUE)) %>%
  mutate(counts_thy = y_pred.1 + y_pred.2 + y_pred.7 + y_pred.8 + y_pred.9 + y_pred.10,
         counts_per = y_pred.3 + y_pred.4 + y_pred.5 + y_pred.6 + y_pred.11 + y_pred.12,
         total_counts = counts_thy + counts_per,
         Nfd_thy = (y_pred.9 + y_pred.10)/(counts_thy * chivec2),
         Nfd_per = (y_pred.11 + y_pred.12)/(counts_per * chivec2),
         donor_ki_thy = (y_pred.9)/(y_pred.9 + y_pred.10),
         donor_ki_per = (y_pred.11)/(y_pred.11 + y_pred.12),
         host_ki_thy = (y_pred.1 + y_pred.7)/(y_pred.1 + y_pred.2 + y_pred.7 + y_pred.8),
         host_ki_per = (y_pred.3 + y_pred.5)/(y_pred.3 + y_pred.4 + y_pred.5 + y_pred.6),
         ageBMT_bin = 'agebin2') %>%
  select(time_seq, ageBMT_bin, contains('thy'), contains('per'))

stan_pred_df3 <- data.frame("time_seq" = ts_pred3,
                            "y_pred" = matrix(unlist(ode_pred3), nrow = length(ts_pred3), byrow = TRUE)) %>%
  mutate(counts_thy = y_pred.1 + y_pred.2 + y_pred.7 + y_pred.8 + y_pred.9 + y_pred.10,
         counts_per = y_pred.3 + y_pred.4 + y_pred.5 + y_pred.6 + y_pred.11 + y_pred.12,
         total_counts = counts_thy + counts_per,
         Nfd_thy = (y_pred.9 + y_pred.10)/(counts_thy * chivec3),
         Nfd_per = (y_pred.11 + y_pred.12)/(counts_per * chivec3),
         donor_ki_thy = (y_pred.9)/(y_pred.9 + y_pred.10),
         donor_ki_per = (y_pred.11)/(y_pred.11 + y_pred.12),
         host_ki_thy = (y_pred.1 + y_pred.7)/(y_pred.1 + y_pred.2 + y_pred.7 + y_pred.8),
         host_ki_per = (y_pred.3 + y_pred.5)/(y_pred.3 + y_pred.4 + y_pred.5 + y_pred.6),
         ageBMT_bin = 'agebin3') %>%
  select(time_seq, ageBMT_bin, contains('thy'), contains('per'))


stan_pred_df4 <- data.frame("time_seq" = ts_pred4,
                            "y_pred" = matrix(unlist(ode_pred4), nrow = length(ts_pred4), byrow = TRUE)) %>%
  mutate(counts_thy = y_pred.1 + y_pred.2 + y_pred.7 + y_pred.8 + y_pred.9 + y_pred.10,
         counts_per = y_pred.3 + y_pred.4 + y_pred.5 + y_pred.6 + y_pred.11 + y_pred.12,
         total_counts = counts_thy + counts_per,
         Nfd_thy = (y_pred.9 + y_pred.10)/(counts_thy * chivec4),
         Nfd_per = (y_pred.11 + y_pred.12)/(counts_per * chivec4),
         donor_ki_thy = (y_pred.9)/(y_pred.9 + y_pred.10),
         donor_ki_per = (y_pred.11)/(y_pred.11 + y_pred.12),
         host_ki_thy = (y_pred.1 + y_pred.7)/(y_pred.1 + y_pred.2 + y_pred.7 + y_pred.8),
         host_ki_per = (y_pred.3 + y_pred.5)/(y_pred.3 + y_pred.4 + y_pred.5 + y_pred.6),
         ageBMT_bin = 'agebin4') %>%
  select(time_seq, ageBMT_bin, contains('thy'), contains('per'))



Stan_pred <- rbind(stan_pred_df1, stan_pred_df2, stan_pred_df3, stan_pred_df4)

ggplot() +
  geom_hline(yintercept = thy0, col='darkred', size=2)+
  geom_hline(yintercept = per0, col='navy', size=2)+
  geom_line(data = Stan_pred, aes(x = time_seq, y = counts_thy, color = ageBMT_bin), linetype=2) +
  geom_line(data = Stan_pred, aes(x = time_seq, y = counts_per, color = ageBMT_bin), linetype=1) +
  geom_point(data = counts_data, aes(x = age.at.S1K, y = total_counts, color = ageBMT_bin), size=2) +
  labs(title=paste('Total counts of thymic naive Tregs'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(60, 450) , trans="log10", breaks=c(10, 30, 100, 300))+
  scale_y_log10() 

ggplot() +
  geom_point(data = Nfd_data, aes(x = age.at.S1K, y = Nfd, col=ageBMT_bin)) +
  geom_line(data = Stan_pred, aes(x = time_seq, y = Nfd_thy, color = ageBMT_bin), size=0.7) +
  labs(title=paste('Nfd thymic naive Tregs'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(60, 450), breaks=c(10, 30, 100, 300)) +
  facet_wrap(.~ location)

ggplot() +
  geom_point(data = Nfd_data, aes(x = age.at.S1K, y = Nfd, col=ageBMT_bin)) +
  geom_line(data = Stan_pred, aes(x = time_seq, y = Nfd_per, color = ageBMT_bin), size=0.7) +
  labs(title=paste('Nfd thymic naive Tregs'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(60, 450), breaks=c(10, 30, 100, 300)) +
  facet_wrap(.~ location)



ggplot() +
  geom_point(data=donorki_data, aes(x = age.at.S1K, y = Ki67_naiveTregs_thy), col=2) +
  geom_line(data = Stan_pred, aes(x = time_seq, y = donor_ki_thy), col=2) +
  geom_point(data=hostki_data, aes(x = age.at.S1K, y = Ki67_naiveTregs_thy), col=4) +
  geom_line(data = Stan_pred, aes(x = time_seq, y = host_ki_thy), col=4) +
  labs(title=paste('Nfd peripheral naive Tregs'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(60, 450) , trans="log10", breaks=c(10, 30, 100, 300)) +
  facet_wrap(.~ ageBMT_bin) + ylim(0, 0.5)


ggplot() +
  geom_point(data=donorki_data, aes(x = age.at.S1K, y = Ki67_naiveTregs_periph), col=2) +
  geom_line(data = Stan_pred, aes(x = time_seq, y = donor_ki_per), col=2) +
  geom_point(data=hostki_data, aes(x = age.at.S1K, y = Ki67_naiveTregs_periph), col=4) +
  geom_line(data = Stan_pred, aes(x = time_seq, y = host_ki_per), col=4) +
  labs(title=paste('Nfd peripheral naive Tregs'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(60, 450) , trans="log10", breaks=c(10, 30, 100, 300)) +
  facet_wrap(.~ ageBMT_bin) + ylim(0, 0.5)






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
#


### R ode solver predictions
Theta_spline <- function(Time, psi){psi * (10^6.407133491 * exp(-0.002387866 * (Time - 49)))}
chi_spline <- function(Time){ifelse(Time-10 >0, 
                                    0.81548689 * (1 - exp(-0.06286984 * (Time - 10))), 0)}
Donor_eps_spline <- function(Time){0}#exp(- 0.06799028 * (Time - 49)) + 0.37848471}


theta_stanvec <- sapply(ts_pred1, theta_spline, psi=0.011)
theta_Rvec <- sapply(ts_pred1, Theta_spline, psi=0.011)

ggplot()+
  geom_line(aes(x=ts_pred1, y=theta_stanvec), col=2, size=3)+
  geom_line(aes(x=ts_pred1, y=theta_Rvec), col=4)+
  scale_y_log10()

Chi_stanvec <- sapply(ts_pred1,Chi_spline)
Chi_Rvec <- sapply(ts_pred1, chi_spline)

ggplot()+
  geom_line(aes(x=ts_pred1, y=Chi_stanvec), col=2, size=3)+
  geom_line(aes(x=ts_pred1, y=Chi_Rvec), col=4)



shm_chi <- function (Time, ageatBMT, y, parms) {
  eps_host = 0.326611
  kloss = 1/3.5
  with(as.list(c(y, parms)),{
    dy1 <- Theta_spline(Time, psi) * (1- chi_spline(Time - ageatBMT)) * eps_host + rho_D * (2 * y2 + y1) + beta * y3 - (kloss + alpha + delta_D) * y1
    dy2 <- Theta_spline(Time, psi) * (1- chi_spline(Time - ageatBMT)) * (1 - eps_host) + kloss * y1 + beta * y4  - (rho_D + alpha + delta_D) * y2
    dy3 <- alpha * y1 + rho_D * (2 * y4 + y3) - (kloss + beta + delta_D) * y3
    dy4 <- alpha * y2 + kloss * y3 - (rho_D + beta + delta_D) * y4
    dy5 <- alpha * y7 + rho_I * (2 * y6 + y5) - (kloss + beta + rho_I) * y5
    dy6 <- alpha * y8 + kloss * y5 - (rho_I + beta + rho_I) * y6
    dy7 <- beta * y5 + rho_I * (2 * y8 + y7) - (kloss + alpha + rho_I) * y7
    dy8 <- beta * y6 + kloss * y7 - (rho_I + alpha + rho_I) * y8
    
    dy9 <- Theta_spline(Time, psi) * chi_spline(Time - ageatBMT) * donor_eps_spline(Time)+ rho_D * (2 * y10 + y9) + beta * y11  - (kloss + alpha + delta_D) * y9
    dy10 <- Theta_spline(Time, psi) * chi_spline(Time - ageatBMT) * (1 - donor_eps_spline(Time)) + kloss * y9 + beta * y12  - (rho_D + alpha + delta_D) * y10
    dy11 <- alpha * y9 + rho_D * (2 * y12 + y11) - (kloss + beta + delta_D) * y11
    dy12 <- alpha * y10 + kloss * y11 - (rho_D + beta + delta_D) * y12
    
    list(c(dy1, dy2, dy3, dy4, dy5, dy6, dy7, dy8,
      dy9, dy10, dy11, dy12))
  })
}



ode(y=init_cond, times=c(40, 45), func=shm_chi, parms=params, ageatBMT=40)[2,2:13]
solve_init(45, init_cond, params)

solve_ode_chi(66, 45, init_cond, params)
ode(y=init_cond1, times=c(45, 66), func=shm_chi, parms=params, ageatBMT=45)[2,2:13]

init_predstan1 <- solve_init(45, init_cond, params)
init_stan1 <- c(init_predstan1[1] + init_predstan1[9], init_predstan1[2] + init_predstan1[10], init_predstan1[3] + init_predstan1[11],
                init_predstan1[4] + init_predstan1[12], init_predstan1[5], init_predstan1[6], init_predstan1[7], init_predstan1[8],
                y9=0,y10=0,y11=0,y12=0)

init_pred1 <- ode(y=init_cond, times=c(40, 45), func=shm_chi, parms=params, ageatBMT=40)[2,2:13]
init_cond1 <- c(init_pred1[1] + init_pred1[9], init_pred1[2] + init_pred1[10], init_pred1[3] + init_pred1[11],
                init_pred1[4] + init_pred1[12], init_pred1[5], init_pred1[6], init_pred1[7], init_pred1[8],
                y9=0,y10=0,y11=0,y12=0)
R_ode_pred1 <- data.frame(ode(y=init_cond1, times=c(45, ts_pred1), func=shm_chi, parms=params, ageatBMT=45)) %>%
  mutate(time_seq = time,
         counts_thy = y1 + y2 + y7 + y8 + y9 + y10,
         counts_per = y3 + y4 + y5 + y6 + y11 + y12,
         total_counts = counts_thy + counts_per,
         Nfd_thy = (y9 + y10)/(counts_thy * chivec4),
         Nfd_per = (y11 + y12)/(counts_per * chivec4),
         donor_ki_thy = (y9)/(y9 + y10),
         donor_ki_per = (y11)/(y11 + y12),
         host_ki_thy = (y1 + y7)/(y1 + y2 + y7 + y8),
         host_ki_per = (y3 + y5)/(y3 + y4 + y5 + y6),
         ageBMT_bin = 'agebin1') %>%
  select(time_seq, ageBMT_bin, contains('thy'), contains('per'))

init_pred2 <- ode(y=init_cond, times=c(40, 66), func=shm_chi, parms=params, ageatBMT=40)[2,2:13]
init_cond2 <- c(init_pred2[1] + init_pred2[9], init_pred2[2] + init_pred2[10], init_pred2[3] + init_pred2[11],
                init_pred2[4] + init_pred2[12], init_pred2[5], init_pred2[6], init_pred2[7], init_pred2[8],
                y9=0,y10=0,y11=0,y12=0)

R_ode_pred2 <- data.frame(ode(y=init_cond2,  times=c(66, ts_pred2), func=shm_chi, parms=params, ageatBMT=66)) %>%
  mutate(time_seq = time,
         counts_thy = y1 + y2 + y7 + y8 + y9 + y10,
         counts_per = y3 + y4 + y5 + y6 + y11 + y12,
         total_counts = counts_thy + counts_per,
         Nfd_thy = (y9 + y10)/(counts_thy * chivec4),
         Nfd_per = (y11 + y12)/(counts_per * chivec4),
         donor_ki_thy = (y9)/(y9 + y10),
         donor_ki_per = (y11)/(y11 + y12),
         host_ki_thy = (y1 + y7)/(y1 + y2 + y7 + y8),
         host_ki_per = (y3 + y5)/(y3 + y4 + y5 + y6),
         ageBMT_bin = 'agebin2') %>%
  select(time_seq, ageBMT_bin, contains('thy'), contains('per'))


init_pred3 <-  ode(y=init_cond, times=c(40, 76), func=shm_chi, parms=params, ageatBMT=40)[2,2:13]
init_cond3 <- c(init_pred3[1] + init_pred3[9], init_pred3[2] + init_pred3[10], init_pred3[3] + init_pred3[11],
                init_pred3[4] + init_pred3[12], init_pred3[5], init_pred3[6], init_pred3[7], init_pred3[8],
                y9=0,y10=0,y11=0,y12=0)

R_ode_pred3 <-data.frame(ode(y=init_cond3,  times=c(76, ts_pred3), func=shm_chi, parms=params, ageatBMT=76)) %>%
  mutate(time_seq = time,
         counts_thy = y1 + y2 + y7 + y8 + y9 + y10,
         counts_per = y3 + y4 + y5 + y6 + y11 + y12,
         total_counts = counts_thy + counts_per,
         Nfd_thy = (y9 + y10)/(counts_thy * chivec4),
         Nfd_per = (y11 + y12)/(counts_per * chivec4),
         donor_ki_thy = (y9)/(y9 + y10),
         donor_ki_per = (y11)/(y11 + y12),
         host_ki_thy = (y1 + y7)/(y1 + y2 + y7 + y8),
         host_ki_per = (y3 + y5)/(y3 + y4 + y5 + y6),
         ageBMT_bin = 'agebin3') %>%
  select(time_seq, ageBMT_bin, contains('thy'), contains('per'))


init_pred4 <-  ode(y=init_cond, times=c(40, 118), func=shm_chi, parms=params, ageatBMT=40)[2,2:13]
init_cond4 <- c(init_pred4[1] + init_pred4[9], init_pred4[2] + init_pred4[10], init_pred4[3] + init_pred4[11],
                init_pred4[4] + init_pred4[12], init_pred4[5], init_pred4[6], init_pred4[7], init_pred4[8],
                y9=0,y10=0,y11=0,y12=0)

R_ode_pred4 <- data.frame(ode(y=init_cond4,  times=c(118, ts_pred4), func=shm_chi, parms=params, ageatBMT=118)) %>%
  mutate(time_seq = time,
         counts_thy = y1 + y2 + y7 + y8 + y9 + y10,
         counts_per = y3 + y4 + y5 + y6 + y11 + y12,
         total_counts = counts_thy + counts_per,
         Nfd_thy = (y9 + y10)/(counts_thy * chivec4),
         Nfd_per = (y11 + y12)/(counts_per * chivec4),
         donor_ki_thy = (y9)/(y9 + y10),
         donor_ki_per = (y11)/(y11 + y12),
         host_ki_thy = (y1 + y7)/(y1 + y2 + y7 + y8),
         host_ki_per = (y3 + y5)/(y3 + y4 + y5 + y6),
         ageBMT_bin = 'agebin4') %>%
  select(time_seq, ageBMT_bin, contains('thy'), contains('per'))


R_pred <- rbind(R_ode_pred1, R_ode_pred2, R_ode_pred3, R_ode_pred4)


ggplot() +
  geom_line(data = Stan_pred, aes(x = time_seq, y = counts_thy, color = ageBMT_bin), linetype=2, size=2) +
  geom_line(data = R_pred, aes(x = time_seq, y = counts_thy, color = ageBMT_bin), linetype=2,) +
  geom_line(data = Stan_pred, aes(x = time_seq, y = counts_per, color = ageBMT_bin), linetype=1, size=2) +
  geom_line(data = R_pred, aes(x = time_seq, y = counts_per, color = ageBMT_bin), linetype=1) +
  labs(title=paste('Total counts of thymic naive Tregs'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(60, 450) , trans="log10", breaks=c(10, 30, 100, 300))+
  scale_y_log10() 

ggplot() +
  #geom_point(data = Nfd_data, aes(x = age.at.S1K, y = Thymus, col=ageBMT_bin)) +
  geom_point(data = R_pred, aes(x = time_seq, y = Nfd_thy, color = ageBMT_bin), size=0.7) +
  labs(title=paste('Nfd thymic naive Tregs'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(60, 450), breaks=c(10, 30, 100, 300))

ggplot() +
  geom_point(data = Nfd_data, aes(x = age.at.S1K, y = Periphery, col=ageBMT_bin)) +
  geom_point(data = Counts_pred, aes(x = time_seq, y = Nfd_per, color = ageBMT_bin)) +
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





