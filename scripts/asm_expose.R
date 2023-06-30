rm(list = ls()); gc()

library(rstan)
library(tidyverse)

## model
rstan::expose_stan_functions("stan_models/MAP_asm_deltavar.stan")
params <- c(1.363198e+05, 6.041325e-03, 1.708735e-02)
theta <- c(0)
x_i <- c(60)
x_r <- c(100)
math_reduce(params, theta, x_r, x_i)

## ts
ts_pred_chi1 <- 10^seq(log10(58), log10(450), length.out = 300)
ts_pred_chi2 <- 10^seq(log10(75), log10(450), length.out = 300)
ts_pred_chi3 <- 10^seq(log10(101), log10(450), length.out = 300)
tb_pred1 <- rep(54, 300)
tb_pred2 <- rep(71, 300)
tb_pred3 <- rep(97, 300)

chi_vec <- sapply(ts_pred_chi1 - tb_pred1, Chi_spline)
total_counts <- N_pooled_time(ts_pred_chi1, tb_pred1, params)
donor_counts <- N_donor_time(ts_pred_chi1, tb_pred1, params)
Nfd_pred1 <- donor_counts/(total_counts * chi_vec)

ggplot()+
  geom_point(aes(x = ts_pred_chi1, y=Nfd_pred1), size=2) + ylim(0,1)

## delta varying ASM in R
#peripheral turnover rate lambda
lambda_exp <- function(a,l,r) {
  l*exp(-r*a)
}

sp_counts <- function(Time){
  ta=1;
  theta0  =  4.3E5;  theta_f = 1.8E3;  n = 2.1;   X = 30.0;   q = 3.7;
  ifelse(Time < ta, theta0,
         theta0 + (theta_f * (Time - ta)^n) * (1 - (((Time - ta)^q)/((X^q) + ((Time - ta)^q))))
  )
}

chi_func <- function(Time){
  chiEst = 0.83155297;
  qEst = 0.047555868;
  chiEst * (1 - exp(-qEst * Time))
 }

#Initial age distribution g(a)

g <- function(a, parms) { 
  N0 = parms[1]
  ifelse(a>=0 && a<=1, 10^N0, 0)
}

#total thymic output
theta <- function(Time, parms) {
  N0 = parms[1]
  psi = 10^N0/sp_counts(1)
  #print(psi)
  return(psi * sp_counts(Time))
}

G <- function(a,tb,parms) { 
  l = exp(parms[2])
  r = exp(parms[3])
  ta=1
  
  ifelse(a < tb-ta, theta(tb-ta-a,parms)  * exp(- l/r* (1-exp(-r*a))), 0)
        # ifelse(a >= tb - ta, g(a-tb+ta,parms) * exp(- l/r * (exp(-r*(a-tb+ta))-exp(-r*a))), 0))
}  

G_t <- function(a,t,tb,parms){ 
  l = exp(parms[2])
  r = exp(parms[3])
  ta=1
  
  ifelse(a<= t-tb,theta(t-ta-a,parms) * exp(- l/r* (1-exp(-r*a))),
         ifelse(a<=t, G(a-t+tb,tb,parms) * exp(- l/r* (exp(-r*(a-t+tb))-exp(-r*a))), 0))
} 


N_Gt <- function(t,tb,parms){
  integrate(G_t, lower = 0, upper = t, t=t,tb=tb, parms=parms)$value
}

N_counts <- function(Time, t0, parms){
  mapply(N_Gt, t=Time, tb=t0, MoreArgs = list(parms=parms))
}

pars <- c(log10(1.363198e+05), log(6.041325e-03), log(1.708735e-02))

#normalised donor fraction
theta_d <- function(Time, parms) {
  tC = 42
  m = 2
  ifelse (Time <= tC, theta(Time,parms)* (Time/tC)^m,
          theta(Time,parms))
}

theta_td_exp <- function(a,t,tb, parms) {
  l = exp(parms[2])
  r = exp(parms[3])
  ta=1
  ifelse(a <= t-tb, theta_d(t-ta-a,parms) * exp(- l/r* (1-exp(-r*a))), 0)
}

Ndonor_exp <- function(t,tb, parms) {
  integrate(theta_td_exp, lower= 0, upper=t-tb, t=t, tb=tb, parms=parms)$value/
    N_Gt(t, tb, parms)
  
}

Ncounts_exp <- function(Time, t0, params){
  
  func_value <- function(Time, t0, params){
    integrate(theta_td_exp, lower= 0, upper=t-tb, t=t, tb=tb, parms=parms)$value/
      N_Gt(t, tb, parms)
  }
  return(mapply(func_value, Time, t0, MoreArgs = list(params)))
}

library(parallel)
ncounts <- mcmapply(N_Gt, t=chimera_data$age.at.S1K, tb=chimera_data$age.at.BMT, MoreArgs = list(pars), mc.cores=4 )
Ncounts_stan <- N_pooled_time(chimera_data$age.at.S1K, chimera_data$age.at.BMT, params)
ggplot()+
  geom_point(aes(x = chimera_data$age.at.S1K, y=ncounts), col=2, size=2) +
  geom_point(aes(x = chimera_data$age.at.S1K, y=Ncounts_stan), col=4, size=2)
  
ncounts1 <- sapply(ts_pred_chi1, N_Gt, t0=tb_pred1[1], parms=pars)
ncounts2 <- sapply(ts_pred_chi2, N_Gt, tb=tb_pred2[1], parms=pars)
ncounts3 <- sapply(ts_pred_chi3, N_Gt, tb=tb_pred3[1], parms=pars)

donor_counts1 <- sapply(ts_pred_chi1, Ndonor_exp, tb=tb_pred1[1], parms=pars)
donor_counts2 <- sapply(ts_pred_chi2, Ndonor_exp, tb=tb_pred2[1], parms=pars)
donor_counts3 <- sapply(ts_pred_chi3, Ndonor_exp, tb=tb_pred3[1], parms=pars)

chi_vec1 <- chi_func(ts_pred_chi1 - tb_pred1)
chi_vec2 <- chi_func(ts_pred_chi2 - tb_pred2)
chi_vec3 <- chi_func(ts_pred_chi3 - tb_pred3)

Nfd1 = donor_counts1/(chi_vec1)
Nfd2 = donor_counts2/(chi_vec2)
Nfd3 = donor_counts3/(chi_vec3)

ggplot()+
  geom_point(aes(x = ts_pred_chi1, y=donor_counts1), col=2, size=2) +
  geom_point(aes(x = ts_pred_chi2, y=donor_counts1), col=4, size=2) +
  geom_point(aes(x = ts_pred_chi3, y=donor_counts1), col=6, size=2) +
  geom_point(aes(x = ts_pred_chi1, y=Nfd_pred1), col=3, size=2) +
  ylim(0, 1.02)


ggplot()+
  geom_point(aes(x = ts_pred_chi1, y=Nfd1), col=2, size=2) +
  geom_point(aes(x = ts_pred_chi2, y=Nfd2), col=4, size=2) +
  geom_point(aes(x = ts_pred_chi3, y=Nfd3), col=6, size=2) +
  geom_point(aes(x = ts_pred_chi1, y=Nfd_pred1), col=3, size=2) 

asinsqrt <- function(x){
  a=1.2
  return(asin(sqrt(x)/sqrt(a)))
}


logitTransformation <- function(x){
  b = 1.2
  answer = log(x/(b-x));
}
#fitting log N_total & logit N_fd to the NIMR data
LL_NCD4_NIMR <- function(param, boot_data) { 
  k  <- length(param)        #numner of unknown parameters 
  n1 <- nrow(boot_data)           #number of observations in dataset1
  n2 <- nrow(boot_data)           #number of observations in dataset2
  n  <- n1+n2
  
  ncounts_pred <- mcmapply(N_Gt, t=chimera_data$age.at.S1K, tb=chimera_data$age.at.BMT, MoreArgs = list(param), mc.cores=1)
  nfd_pred <- mcmapply(Ndonor_exp, t=chimera_data$age.at.S1K, tb=chimera_data$age.at.BMT, MoreArgs = list(param), mc.cores=1 )
  
  R1 <- sum((log(boot_data$total_counts) - log(ncounts_pred))^2)   #SSR for dataset1
  R2 <- sum((logitTransformation(boot_data$Nfd) - logitTransformation(nfd_pred))^2)  #SSR for dataset2
  
  #log-likelihood ignoring all the terms dependent only on the number of observations n
  #matrix multipltication of residual and transpose of residuals
  logl <- -(n1/2)*log(R1) - (n2/2)*log(R2)
  
  aiccd4NIMR <<- -2*logl + 2*k
  
  return(-logl)     #since optim minimizes the function by default, ML
} 

fit_NCD4_NIMR <- optim(par=c(5.4, -3, -4.5), fn=LL_NCD4_NIMR, boot_data=chimera_data)
fit_NCD4_NIMR
aiccd4NIMR

ncounts <- mcmapply(N_Gt, t=ts_pred_chi1, tb=tb_pred1, 
                    MoreArgs = list(fit_NCD4_NIMR$par), mc.cores=4 )
nfd1 <- mcmapply(Ndonor_exp, t=ts_pred_chi1, tb=tb_pred1, 
                    MoreArgs = list(fit_NCD4_NIMR$par), mc.cores=4 )
nfd2 <- mcmapply(Ndonor_exp, t=ts_pred_chi2, tb=tb_pred2, 
                MoreArgs = list(fit_NCD4_NIMR$par), mc.cores=4 )
nfd3 <- mcmapply(Ndonor_exp, t=ts_pred_chi3, tb=tb_pred3, 
                MoreArgs = list(fit_NCD4_NIMR$par), mc.cores=4 )


ggplot()+
  geom_line(aes(x = ts_pred_chi1, y=ncounts), col=2, size=1.5) +
  geom_point(aes(x = chimera_data$age.at.S1K, y=chimera_data$total_counts), col=2, size=2)+
  scale_y_log10()


##### main plots
ageBMT_names <- c(`ageBMT_group1` = "7-9wks",
                  `ageBMT_group2` = "9-11wks", 
                  `ageBMT_group3` = "11-25wks")

ggplot()+
  geom_line(aes(x = ts_pred_chi1, y=nfd1), col=2, size=1.5) +
  geom_line(aes(x = ts_pred_chi2, y=nfd2), col=7, size=1.5) +
  geom_line(aes(x = ts_pred_chi3, y=nfd3), col=4, size=1.5) +
  geom_point(aes(x = chimera_data$age.at.S1K, y=chimera_data$Nfd, col=chimera_data$ageBMT_bin), size=2)+
  scale_color_manual(name= "Age at BMT", values = c(2, 7, 4), labels = ageBMT_names)+
  labs(x = 'Host age', y = NULL, title = paste0("Normalised donor fraction in naive ", toupper(Population), " T cells")) +
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) + 
  scale_x_continuous(limits = c(50, 450), breaks = c(10, 30, 100, 300, 450)) +
  myTheme + theme(legend.position = c(0.85, 0.15))










