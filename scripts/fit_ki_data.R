## clearing the environment
rm(list = ls())  
gc()    

library(rstan)
library(deSolve)
library(tidyverse)
library(parallel)
####################################################################################
## importing data to be fitted 
counts_file <- file.path("data", "Counts_naiTreg.csv")
counts_data <- read.csv(counts_file) %>% 
  arrange(age.at.S1K) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 56, 'agebin1',
                             ifelse(age.at.BMT <= 70, 'agebin2',
                                    ifelse(age.at.BMT <= 84, 'agebin3', 'agebin4')))) %>%
  gather(c(Thymus, Periphery), key='location', value = 'total_counts')

Nfd_file <- file.path("data", "Nfd_naiTreg.csv")
Nfd_data <- read.csv(Nfd_file) %>% 
  arrange(age.at.S1K)%>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 56, 'agebin1',
                             ifelse(age.at.BMT <= 70, 'agebin2',
                                    ifelse(age.at.BMT <= 84, 'agebin3', 'agebin4')))) %>%
  gather(c(Thymus, Periphery), key='location', value = 'Nfd')

hostki_file <- file.path("data", "hostKi67_naiTreg.csv")
hostki_data <- read.csv(hostki_file) %>% 
  arrange(age.at.S1K) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 56, 'agebin1',
                             ifelse(age.at.BMT <= 70, 'agebin2',
                                    ifelse(age.at.BMT <= 84, 'agebin3', 'agebin4'))),
         subcomp='Host') %>%
  gather(c(Thymus, Periphery), key='location', value = 'prop_ki')

donorki_file <- file.path("data", "donorKi67_naiTreg.csv")
donorki_data <- read.csv(donorki_file) %>% 
  arrange(age.at.S1K) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 56, 'agebin1',
                             ifelse(age.at.BMT <= 70, 'agebin2',
                                    ifelse(age.at.BMT <= 84, 'agebin3', 'agebin4'))),
         subcomp='Donor')%>%
  gather(c(Thymus, Periphery), key='location', value = 'prop_ki')

ki_data <- rbind(donorki_data, hostki_data)

## Unique time points with indices to map
data_time_donorki <- donorki_data$age.at.S1K 
data_time_hostki <- hostki_data$age.at.S1K 
tb_time_donorki <- donorki_data$age.at.BMT 
tb_time_hostki <- hostki_data$age.at.BMT

### R ode solver predictions
Theta_spline <- function(Time, psi){psi * (10^6.407133491 * exp(-0.002387866 * (Time - 49)))}
chi_spline <- function(Time){ifelse(Time-10 >0, 
                                    0.81548689 * (1 - exp(-0.06286984 * (Time - 10))), 0)}
Donor_eps_spline <- function(Time){exp(- 0.06799028 * (Time - 49)) + 0.37848471}

chivec1 <- sapply(ts_pred1-tb_pred1, chi_spline)
chivec2 <- sapply(ts_pred2-tb_pred2, chi_spline)
chivec3 <- sapply(ts_pred3-tb_pred3, chi_spline)
chivec4 <- sapply(ts_pred4-tb_pred4, chi_spline)


## ODE System
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
    
    dy9 <- Theta_spline(Time, psi) * chi_spline(Time - ageatBMT) * Donor_eps_spline(Time)+ rho_D * (2 * y10 + y9) + beta * y11  - (kloss + alpha + delta_D) * y9
    dy10 <- Theta_spline(Time, psi) * chi_spline(Time - ageatBMT) * (1 - Donor_eps_spline(Time)) + kloss * y9 + beta * y12  - (rho_D + alpha + delta_D) * y10
    dy11 <- alpha * y9 + rho_D * (2 * y12 + y11) - (kloss + beta + delta_D) * y11
    dy12 <- alpha * y10 + kloss * y11 - (rho_D + beta + delta_D) * y12
    
    list(c(dy1, dy2, dy3, dy4, dy5, dy6, dy7, dy8,
      dy9, dy10, dy11, dy12))
  })
}

## initial conditions at ageatBMTs
sol_init <- function(parms, tb){
  rho_D = parms[1]; rho_I = parms[2];
  f1 = parms[3]; f2 = parms[4]; f3 = parms[5]; f4 = parms[6];
  ta = 40.0;
  
  ## parameter -- estimates and new definitions
  init_cond <- c("y1"= f1 * exp(9.7), "y2"= (1-f1) *  exp(9.7), 
                 "y3"= f2 *  exp(14.5), "y4"= (1-f2) *  exp(14.5), 
                 "y5"= f3 *  exp(9.25), "y6"= (1-f3) *  exp(9.25),
                 "y7"= f4 *  exp(11.45), "y8" = (1-f4) *  exp(11.45), 
                 "y9"=0, "y10"=0, "y11"=0, "y12"=0)
  
  params_est <- c(psi=0.008099747, alpha=0.638427240, rho_D = rho_D,
                  delta_D=0.020695076 + rho_D, rho_I = rho_I, beta=0.011862542)
  
  return(ode(y=init_cond, times=c(ta, tb), func=shm_chi, parms=params_est, ageatBMT=ta)[2,2:13])
}

params_new <- c(0.0004, 0.04, 0.3, 0.6, 0.4, 0.3)

sol_ode_par <- function(Time, tb, parms){
  rho_D = parms[1]; rho_I = parms[2];
  params_est <- c(psi=0.008099747, alpha=0.638427240, rho_D = rho_D,
                  delta_D=0.020695076 + rho_D, rho_I = rho_I, beta=0.011862542)
  
    init_tb <- sol_init(parms = parms, tb = tb)
    
    init_cond <- c(init_tb[1] + init_tb[9], init_tb[2] + init_tb[10], init_tb[3] + init_tb[11],
                   init_tb[4] + init_tb[12], init_tb[5], init_tb[6], init_tb[7], init_tb[8],
                   y9=0,y10=0,y11=0,y12=0)

  return(ode(y=init_cond, times=c(tb, Time), func=shm_chi, parms=params_est, ageatBMT=tb)[2, 2:13])
}


#### data transformation functions
logit_trans <- function(x){
  
  log(x/(1-x))
}

expit_trans <- function(x){
  
  exp(x)/(1 + exp(x))
}



##fitting log N_total & logit ki_prop to the transformed NCD4 data
LL_fit <- function(params, boot_data1, boot_data2) { 
   ## model predictions
  Donor_ode_sol <- mcmapply(sol_ode_par, Time=data_time_donorki, 
                            tb=tb_time_donorki, MoreArgs = list(params), mc.cores = 13)
  
  Donor_sol_df <- data.frame(t(Donor_ode_sol)) %>%
    mutate(donor_ki_thy = (y9)/(y9 + y10),
           donor_ki_per = (y11)/(y11 + y12)) %>%
    select(contains('donor'))
  
  Host_ode_sol<- mcmapply(sol_ode_par, Time=data_time_hostki, 
                          tb=tb_time_hostki, MoreArgs = list(params), mc.cores = 15)
  
  Host_sol_df <- data.frame(t(Host_ode_sol)) %>%
    mutate(host_ki_thy = (y9)/(y9 + y10),
           host_ki_per = (y11)/(y11 + y12)) %>%
    select(contains('host'))
  
  ## logit transformed ki67 proportions from the model
  donor_ki_thy_pred <- logit_trans(Donor_sol_df$donor_ki_thy)
  donor_ki_per_pred <- logit_trans(Donor_sol_df$donor_ki_per)
  host_ki_thy_pred  <- logit_trans(Host_sol_df$host_ki_thy)
  host_ki_per_pred  <- logit_trans(Host_sol_df$host_ki_per)
  
  ### data
  donor_ki_thy_obs <- logit_trans(boot_data1$Ki67_naiveTregs_thy)
  donor_ki_per_obs <- logit_trans(boot_data1$Ki67_naiveTregs_per)
  host_ki_thy_obs  <- logit_trans(boot_data2$Ki67_naiveTregs_thy)
  host_ki_per_obs  <- logit_trans(boot_data2$Ki67_naiveTregs_per)
  
  ## calculating the sun of squared residuals
  ssqres1 <- sum((donor_ki_thy_obs - donor_ki_thy_pred)^2)
  ssqres2 <- sum((donor_ki_per_obs - donor_ki_per_pred)^2)
  ssqres3 <- sum((host_ki_thy_obs - host_ki_thy_pred)^2)
  ssqres4 <- sum((host_ki_per_obs - host_ki_per_pred)^2)
  
  k  <- length(params)                #number of unknown parameters 
  n1 <- nrow(boot_data1)         #number of observations in dataset1
  n2 <- nrow(boot_data2)           #number of observations in dataset1
  
  #cost function
  #log-likelihood ignoring all the terms dependent only on the number of observations n
  #matrix multipltication of residual and transpose of residuals
  logl <-  - (n1/2)* log(ssqres1) - (n1/2)* log(ssqres2) - 
    (n2/2)* log(ssqres3)  - (n2/2)* log(ssqres4)  
  
  return(-logl)     #since optim minimizes the function by default, ML
} 

optim_fit <- optim(par = params_new, fn = LL_fit, boot_data1 = donorki_data, boot_data2=hostki_data,
                   method=c("Nelder-Mead"), control = list(trace = 1, maxit = 2000))

optim_fit$par
par_est <- optim_fit$par

AIC_est <- 2 * length(optim_fit$par) + 2 * optim_fit$value
print(paste0("Estimated parameters are: ", optim_fit$par))
print(paste0("AIC value is: ", AIC_est))


# time sequence for predictions specific to age bins within the data
ts_pred1 <- 10^seq(log10(66), log10(450), length.out = 300)
ts_pred2 <- 10^seq(log10(91), log10(450), length.out = 300)
ts_pred3 <- 10^seq(log10(90), log10(450), length.out = 300)
ts_pred4 <- 10^seq(log10(174), log10(450), length.out = 300)
tb_pred1 <- rep(45, 300)
tb_pred2 <- rep(66, 300)
tb_pred3 <- rep(76, 300)
tb_pred4 <- rep(118, 300)

Donor_ode_pred1 <- mcmapply(sol_ode_par, Time=ts_pred1, tb=tb_pred1, MoreArgs = list(params_new), mc.cores = 13)

Donor_pred1_df <- data.frame(t(Donor_ode_pred1)) %>%
  mutate(time_seq = ts_pred1,
         donor_ki_thy = (y9)/(y9 + y10),
         donor_ki_per = (y11)/(y11 + y12)) %>%
  select(time_seq, contains('donor'))

Host_ode_pred1 <- mcmapply(sol_ode_par, Time=ts_pred1, tb=tb_pred1, MoreArgs = list(params_new), mc.cores = 10)

Host_pred1_df <- data.frame(t(Host_ode_pred1)) %>%
  mutate(time_seq = data_time_hostki,
         donor_ki_thy = (y9)/(y9 + y10),
         donor_ki_per = (y11)/(y11 + y12)) %>%
  select(time_seq, contains('host'))

Donor_ode_pred1 <- mcmapply(sol_ode_par, Time=ts_pred1, tb=tb_pred1, MoreArgs = list(params_new), mc.cores = 13)

Donor_pred1_df <- data.frame(t(Donor_ode_pred1)) %>%
  mutate(time_seq = ts_pred1,
         donor_ki_thy = (y9)/(y9 + y10),
         donor_ki_per = (y11)/(y11 + y12)) %>%
  select(time_seq, contains('donor'))

Host_ode_pred1 <- mcmapply(sol_ode_par, Time=ts_pred1, tb=tb_pred1, MoreArgs = list(params_new), mc.cores = 10)

Host_pred1_df <- data.frame(t(Host_ode_pred1)) %>%
  mutate(time_seq = data_time_hostki,
         donor_ki_thy = (y9)/(y9 + y10),
         donor_ki_per = (y11)/(y11 + y12)) %>%
  select(time_seq, contains('host'))

Donor_ode_pred1 <- mcmapply(sol_ode_par, Time=ts_pred1, tb=tb_pred1, MoreArgs = list(params_new), mc.cores = 13)

Donor_pred1_df <- data.frame(t(Donor_ode_pred1)) %>%
  mutate(time_seq = ts_pred1,
         donor_ki_thy = (y9)/(y9 + y10),
         donor_ki_per = (y11)/(y11 + y12)) %>%
  select(time_seq, contains('donor'))

Host_ode_pred1 <- mcmapply(sol_ode_par, Time=ts_pred1, tb=tb_pred1, MoreArgs = list(params_new), mc.cores = 10)

Host_pred1_df <- data.frame(t(Host_ode_pred1)) %>%
  mutate(time_seq = data_time_hostki,
         donor_ki_thy = (y9)/(y9 + y10),
         donor_ki_per = (y11)/(y11 + y12)) %>%
  select(time_seq, contains('host'))


myTheme <- theme(text = element_text(size = 12), axis.text = element_text(size = 12), axis.title =  element_text(size = 12, face = "bold"),
                 plot.title = element_text(size=12, face = 'bold',  hjust = 0.5), legend.text = element_text(size=12),
                 legend.title = element_text(size = 12, face = 'bold'), legend.background = element_blank())

# setting ggplot theme for rest fo the plots
theme_set(theme_bw())

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove + after exponent, if exists. E.g.: (e^+2 -> e^2)
  l <- gsub("e\\+","e",l)  
  # turn the 'e' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # convert 1x10^ or 1.000x10^ -> 10^
  l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  # return this as an expression
  parse(text=l)
}

log10minorbreaks=as.numeric(1:10 %o% 10^(4:8))

legn_labels <- c('6-8', '8-10', '10-12', '12-25')
fac_labels <- c(`agebin1`= '6-8 weeks', `agebin2`= '8-10 weeks', `agebin3`= '10-12 weeks', `agebin4`= '12-25 weeks')

## Total counts
Counts_pred <- R_pred %>%
  select(time_seq, ageBMT_bin, contains('counts')) %>%
  rename(Thymus = counts_thy,
         Periphery = counts_per) %>%
  gather(c(Thymus, Periphery), key='location', value = 'total_counts')

p1 <- ggplot() +
  geom_line(data = Counts_pred, aes(x = time_seq, y = total_counts, color = ageBMT_bin)) +
  geom_point(data = counts_data, aes(x = age.at.S1K, y = total_counts, color = ageBMT_bin), size=2) +
  labs(title=paste('Total counts of naive Tregs'),  y=NULL, x= "Host age (days)") + 
  scale_color_discrete(name="Host age at \n BMT (Wks)", labels=legn_labels)+
  scale_x_continuous(limits = c(60, 450) , trans="log10", breaks=c(10, 30, 100, 300))+
  scale_y_continuous(limits = c(5e3, 5e6), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  facet_wrap(~ factor(location, levels = c('Thymus', "Periphery")))+
  guides(fill = 'none') + myTheme 


# Normalised donor fractions
Nfd_pred <- R_pred %>%
  select(time_seq, ageBMT_bin, contains('Nfd')) %>%
  rename(Thymus = Nfd_thy,
         Periphery = Nfd_per) %>%
  gather(c(Thymus, Periphery), key='location', value = 'Nfd')

p2 <- ggplot() +
  geom_line(data = Nfd_pred, aes(x = time_seq, y = Nfd, color = ageBMT_bin)) +
  geom_point(data = Nfd_data, aes(x = age.at.S1K, y = Nfd, color = ageBMT_bin), size=2) +
  labs(x = "Host age (days)", y = NULL, title = "Normalised Chimerism in naive Tregs") +
  scale_color_discrete(name="Host age at \n BMT (Wks)", labels=legn_labels)+
  scale_x_continuous(limits = c(60, 450), breaks = c(0,100,200,300, 400, 500))+
  scale_y_continuous(limits =c(0, 1.02), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) + 
  facet_wrap(~ factor(location, levels = c('Thymus', "Periphery")))+
  guides(fill='none')+ myTheme

# Thymic Ki67 fractions
ki_thy_pred <- R_pred %>%
  select(time_seq, ageBMT_bin, contains('ki_thy')) %>%
  rename(Donor = donor_ki_thy,
         Host = host_ki_thy) %>%
  gather(c(Donor, Host), key='subcomp', value = 'prop_ki')


p3 <- ggplot() +
  geom_line(data = ki_thy_pred, aes(x = time_seq, y = prop_ki*100, color = subcomp)) +
  geom_point(data = filter(ki_data, location == "Thymus"), aes(x = age.at.S1K, y = prop_ki*100, color = subcomp), size=1.5) +
  labs(x = "Host age (days)", y = NULL, title = "% Ki67hi in thymic naive Tregs") +
  scale_x_continuous(limits = c(60, 450), breaks = c(0,100,200,300, 400, 500))+
  scale_y_continuous(limits =c(0, 50), breaks = c(0, 10, 20, 30, 40, 50))+ 
  facet_wrap(~ ageBMT_bin, scales = 'free', labeller = as_labeller(fac_labels))+
  guides(fill='none') + myTheme + theme(legend.title = element_blank())

# Peripheral Ki67 fractions
ki_per_pred <- R_pred %>%
  select(time_seq, ageBMT_bin, contains('ki_per')) %>%
  rename(Donor = donor_ki_per,
         Host = host_ki_per) %>%
  gather(c(Donor, Host), key='subcomp', value = 'prop_ki')


p4 <- ggplot() +
  geom_line(data = ki_thy_pred, aes(x = time_seq, y = prop_ki*100, color = subcomp)) +
  geom_point(data = filter(ki_data, location == "Periphery"), aes(x = age.at.S1K, y = prop_ki*100, color = subcomp), size=1.5) +
  labs(x = "Host age (days)", y = NULL, title = "% Ki67hi in peripheral naive Tregs") +
  scale_x_continuous(limits = c(60, 450), breaks = c(0,100,200,300, 400, 500))+
  scale_y_continuous(limits =c(0, 50), breaks = c(0, 10, 20, 30, 40, 50))+ 
  facet_wrap(~ ageBMT_bin, scales = 'free', labeller = as_labeller(fac_labels))+
  guides(fill='none') + myTheme + theme(legend.title = element_blank())

cowplot::plot_grid(p1, p2, nrow = 2)


