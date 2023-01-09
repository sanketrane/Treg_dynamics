## clearing the environment
rm(list = ls())  
gc()    

library(rstan)
library(deSolve)
library(tidyverse)
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

init_cond_full <- c("y1"=6.59029693, "y2"= 7.83911652,
                    "y3"= 6.33645005, "y4"= 8.68433688, 
                    "y5"= 6.89305477, "y6"= 6.99702152,
                    "y7"= 6.38125476, "y8" = 6.6390128,
                    "y9"=0, "y10"=0, "y11"=0, "y12"=0)

init_cond_simple <- c("y1"=exp(9.7), "y2"= exp(14.5), "y3"=exp(9.25), "y4"=exp(11.45), "y5"=0, "y6"=0)

params_full <- c(psi=0.00702977401, rho_D=0.0285546524452, alpha=0.562172401,
                 delta_D=0.0464934056, rho_I=0.14219259882, beta=0.01012512249)
params_simple <- c(psi=0.011, alpha=0.83, delta_D=0.027, delta_I=0, beta=0.012)

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


### R ode solver predictions
Theta_spline <- function(Time, psi){psi * (10^6.407133491 * exp(-0.002387866 * (Time - 49)))}
chi_spline <- function(Time){ifelse(Time-10 >0, 
                                    0.81548689 * (1 - exp(-0.06286984 * (Time - 10))), 0)}
Donor_eps_spline <- function(Time){0}#exp(- 0.06799028 * (Time - 49)) + 0.37848471}

chivec1 <- sapply(ts_pred1-tb_pred1, chi_spline)
chivec2 <- sapply(ts_pred2-tb_pred2, chi_spline)
chivec3 <- sapply(ts_pred3-tb_pred3, chi_spline)
chivec4 <- sapply(ts_pred4-tb_pred4, chi_spline)

## ODE System
shm_simple <- function (Time, ageatBMT, y, parms) {
  with(as.list(c(y, parms)),{
    # Host compartment
    #Thymic displaceable
    dy1 <- Theta_spline(Time, psi) * (1- chi_spline(Time - ageatBMT)) + beta * y2 - (alpha + delta_D) * y1
    #Peripheral displaceable
    dy2 <- alpha * y1 - (beta + delta_D) * y2
    #Peripheral Incumbent
    dy3 <- alpha * y4 - (beta + delta_I) * y3
    #Thymics Incumbent
    dy4 <- beta * y3 - (alpha + delta_I) * y4
    
    # Donor compartment
    #Thymic displaceable
    dy5 <- Theta_spline(Time, psi) * chi_spline(Time - ageatBMT) + beta * y6  - (alpha + delta_D) * y5
    #Peripheral displaceable
    dy6 <- alpha * y5 - (beta + delta_D) * y6
    
    list(c(dy1, dy2, dy3, dy4, dy5, dy6))
  })
}


## ODE System
shm_full <- function (Time, ageatBMT, y, parms) {
  eps_host = 0
  kloss = 1/3.5
  with(as.list(c(y, parms)),{
    # Host compartment
    dy1 <- Theta_spline(Time, psi) * (1- chi_spline(Time - ageatBMT)) * eps_host + rho_D * (2 * y2 + y1) + beta * y3 - (kloss + alpha + delta_D) * y1
    dy2 <- Theta_spline(Time, psi) * (1- chi_spline(Time - ageatBMT)) * (1 - eps_host) + kloss * y1 + beta * y4  - (rho_D + alpha + delta_D) * y2
    dy3 <- alpha * y1 + rho_D * (2 * y4 + y3) - (kloss + beta + delta_D) * y3
    dy4 <- alpha * y2 + kloss * y3 - (rho_D + beta + delta_D) * y4
    dy5 <- alpha * y7 + rho_I * (2 * y6 + y5) - (kloss + beta + rho_I) * y5
    dy6 <- alpha * y8 + kloss * y5 - (rho_I + beta + rho_I) * y6
    dy7 <- beta * y5 + rho_I * (2 * y8 + y7) - (kloss + alpha + rho_I) * y7
    dy8 <- beta * y6 + kloss * y7 - (rho_I + alpha + rho_I) * y8
    
    # Donor compartment
    dy9 <- Theta_spline(Time, psi) * chi_spline(Time - ageatBMT) * Donor_eps_spline(Time)+ rho_D * (2 * y10 + y9) + beta * y11  - (kloss + alpha + delta_D) * y9
    dy10 <- Theta_spline(Time, psi) * chi_spline(Time - ageatBMT) * (1 - Donor_eps_spline(Time)) + kloss * y9 + beta * y12  - (rho_D + alpha + delta_D) * y10
    dy11 <- alpha * y9 + rho_D * (2 * y12 + y11) - (kloss + beta + delta_D) * y11
    dy12 <- alpha * y10 + kloss * y11 - (rho_D + beta + delta_D) * y12
    
    list(c(dy1, dy2, dy3, dy4, dy5, dy6, dy7, dy8,
      dy9, dy10, dy11, dy12))
  })
}


## initial conditions at ageatBMTs
init_pred_simple1 <- ode(y=init_cond_simple, times=c(40, 45), func=shm_simple, parms=params_simple, ageatBMT=40)[2,2:7]
init_cond_simple1 <- c(init_pred_simple1[1] + init_pred_simple1[5], init_pred_simple1[2] + init_pred_simple1[6], init_pred_simple1[3],
                       init_pred_simple1[4], y5=0,y6=0)
R_simple_pred1 <- data.frame(ode(y=init_cond_simple1, times=c(45, ts_pred1), func=shm_full, parms=params_simple, ageatBMT=45)) %>%
  filter(time != 45) %>%
  mutate(time_seq = time,
         counts_thy = y1 + y2 + y7 + y8 + y9 + y10,
         counts_per = y3 + y4 + y5 + y6 + y11 + y12,
         total_counts = counts_thy + counts_per,
         Nfd_thy = (y9 + y10)/(counts_thy * chivec4),
         Nfd_per = (y11 + y12)/(counts_per * chivec4),
         ageBMT_bin = 'agebin1') %>%
  select(time_seq, ageBMT_bin, contains('thy'), contains('per'))

init_pred2 <- ode(y=init_cond, times=c(40, 66), func=shm_full, parms=params, ageatBMT=40)[2,2:13]
init_cond2 <- c(init_pred2[1] + init_pred2[9], init_pred2[2] + init_pred2[10], init_pred2[3] + init_pred2[11],
                init_pred2[4] + init_pred2[12], init_pred2[5], init_pred2[6], init_pred2[7], init_pred2[8],
                y9=0,y10=0,y11=0,y12=0)

R_ode_pred2 <- data.frame(ode(y=init_cond2,  times=c(66, ts_pred2), func=shm_full, parms=params, ageatBMT=66)) %>%
  filter(time != 66) %>%
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


init_pred3 <-  ode(y=init_cond, times=c(40, 76), func=shm_full, parms=params, ageatBMT=40)[2,2:13]
init_cond3 <- c(init_pred3[1] + init_pred3[9], init_pred3[2] + init_pred3[10], init_pred3[3] + init_pred3[11],
                init_pred3[4] + init_pred3[12], init_pred3[5], init_pred3[6], init_pred3[7], init_pred3[8],
                y9=0,y10=0,y11=0,y12=0)

R_ode_pred3 <-data.frame(ode(y=init_cond3,  times=c(76, ts_pred3), func=shm_full, parms=params, ageatBMT=76)) %>%
  filter(time != 76) %>%
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


init_pred4 <-  ode(y=init_cond, times=c(40, 118), func=shm_full, parms=params, ageatBMT=40)[2,2:13]
init_cond4 <- c(init_pred4[1] + init_pred4[9], init_pred4[2] + init_pred4[10], init_pred4[3] + init_pred4[11],
                init_pred4[4] + init_pred4[12], init_pred4[5], init_pred4[6], init_pred4[7], init_pred4[8],
                y9=0,y10=0,y11=0,y12=0)

R_ode_pred4 <- data.frame(ode(y=init_cond4,  times=c(118, ts_pred4), func=shm_full, parms=params, ageatBMT=118)) %>%
  filter(time != 118) %>%
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

## initial conditions at ageatBMTs
init_pred1 <- ode(y=init_cond_full, times=c(40, 45), func=shm_full, parms=params_full, ageatBMT=40)[2,2:13]
init_cond1 <- c(init_pred1[1] + init_pred1[9], init_pred1[2] + init_pred1[10], init_pred1[3] + init_pred1[11],
                init_pred1[4] + init_pred1[12], init_pred1[5], init_pred1[6], init_pred1[7], init_pred1[8],
                y9=0,y10=0,y11=0,y12=0)
R_ode_pred1 <- data.frame(ode(y=init_cond1, times=c(45, ts_pred1), func=shm_full, parms=params_full, ageatBMT=45)) %>%
  filter(time != 45) %>%
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

init_pred2 <- ode(y=init_cond_full, times=c(40, 66), func=shm_full, parms=params_full, ageatBMT=40)[2,2:13]
init_cond2 <- c(init_pred2[1] + init_pred2[9], init_pred2[2] + init_pred2[10], init_pred2[3] + init_pred2[11],
                init_pred2[4] + init_pred2[12], init_pred2[5], init_pred2[6], init_pred2[7], init_pred2[8],
                y9=0,y10=0,y11=0,y12=0)

R_ode_pred2 <- data.frame(ode(y=init_cond2,  times=c(66, ts_pred2), func=shm_full, parms=params_full, ageatBMT=66)) %>%
  filter(time != 66) %>%
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


init_pred3 <-  ode(y=init_cond_full, times=c(40, 76), func=shm_full, parms=params_full, ageatBMT=40)[2,2:13]
init_cond3 <- c(init_pred3[1] + init_pred3[9], init_pred3[2] + init_pred3[10], init_pred3[3] + init_pred3[11],
                init_pred3[4] + init_pred3[12], init_pred3[5], init_pred3[6], init_pred3[7], init_pred3[8],
                y9=0,y10=0,y11=0,y12=0)

R_ode_pred3 <-data.frame(ode(y=init_cond3,  times=c(76, ts_pred3), func=shm_full, parms=params_full, ageatBMT=76)) %>%
  filter(time != 76) %>%
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


init_pred4 <-  ode(y=init_cond_full, times=c(40, 118), func=shm_full, parms=params_full, ageatBMT=40)[2,2:13]
init_cond4 <- c(init_pred4[1] + init_pred4[9], init_pred4[2] + init_pred4[10], init_pred4[3] + init_pred4[11],
                init_pred4[4] + init_pred4[12], init_pred4[5], init_pred4[6], init_pred4[7], init_pred4[8],
                y9=0,y10=0,y11=0,y12=0)

R_ode_pred4 <- data.frame(ode(y=init_cond4,  times=c(118, ts_pred4), func=shm_full, parms=params_full, ageatBMT=118)) %>%
  filter(time != 118) %>%
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


