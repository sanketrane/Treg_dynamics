## Fitting splines to SP4 data to estimate the means

## cleaning environment and cache
rm(list = ls()); gc()

#Loading required libraries
library(tidyverse)

# importing data ---------------------------------
### We consider thymic FoxP3 negative SP4 as the source for Tregs (naive Tregs) 
### for all purposes henceforth -- SP4 = FoxP3 negative SP4

## loading pre-processed data and munging
### different data frames because the number of viable observations per dataset are different
SP4_counts_df <- read_csv("data/Counts_thymicSource.csv") %>%
  rename(SP4_counts = FoxP3_Neg_SP4) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 63, 'agegroup1',
                             ifelse(age.at.BMT <= 88, 'agegroup2',
                                    'agegroup3')))

SP4_chimerism_df <- read_csv("data/Chimerism_thymicSource.csv") %>%
  rename(SP4_chi = FoxP3_Neg_SP4) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 63, 'agegroup1',
                             ifelse(age.at.BMT <= 88, 'agegroup2',
                                    'agegroup3'))) %>%
  group_by(ageBMT_bin) %>%
  mutate(MeanAgeBMT = mean(age.at.BMT)) %>%
  filter(mouse.ID != 213374,
         mouse.ID != 202893)

SP4_donorki67_df <- read_csv("data/donorKi67_thymicSource.csv") %>%
  rename(SP4_ki = FoxP3_Neg_SP4) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 63, 'agegroup1',
                             ifelse(age.at.BMT <= 88, 'agegroup2',
                                    'agegroup3')))
SP4_hostki67_df <- read_csv("data/hostKi67_thymicSource.csv") %>%
  rename(SP4_ki = FoxP3_Neg_SP4) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 63, 'agegroup1',
                             ifelse(age.at.BMT <= 88, 'agegroup2',
                                    'agegroup3')))


## vectors depicting timeseq for predictions 
timeseq <- seq(63, 450)  ## host age
dpt_seq <- seq(14, 300)   ## days post BMT
  

## fd_fit
## phenomenological function
chi_spline <- function(Time, chiEst, qEst){
  Chi = ifelse((Time - 14) < 0, 0,
               chiEst * (1-exp(-qEst * (Time - 14))))
  return(Chi)
}

## LL function
SP4_chi_nlm <- nls(SP4_chi ~ chi_spline(age.at.S1K - age.at.BMT, chiEst, qEst),
                  data =  SP4_chimerism_df,
                  start = list(chiEst=0.84, qEst=0.1))
SP4_chi_pars <- coef(SP4_chi_nlm)

# prediction
SP4chi_fit <- data.frame(dpt_seq, "y_sp" = chi_spline(dpt_seq, SP4_chi_pars[1], SP4_chi_pars[2]))
SP4chi_fit_m <- data.frame(dpt_seq, "y_sp" = chi_spline(dpt_seq, 0.8, SP4_chi_pars[2]))
ggplot() + 
  geom_point(data=SP4_chimerism_df, aes(x=age.at.S1K-age.at.BMT, y=SP4_chi, col=ageBMT_bin), size =2) +
  geom_line(data = SP4chi_fit, aes(x = dpt_seq , y = y_sp), col=4, size =1) + 
  #geom_line(data = SP4chi_fit_m, aes(x = dpt_seq , y = y_sp), col=4, size =1) + 
  #geom_text(data = SP4_chimerism_df, aes(x = age.at.S1K-age.at.BMT , y = SP4_chi, label=mouse.ID), col=4, size =4) + 
  ylim(0,1) + #scale_x_log10(limits=c(10, 750)) +
  labs(title = 'Donor fraction in FoxP3 negative SP4 cells',  y=NULL,  x = 'Days post BMT') 

## Counts fit
## phenomenological function
counts_spline <- function(age.at.S1K, basl, nu){
  Time = age.at.S1K   ###
  t0 = 49
  return(10^basl * exp(-nu * (Time - t0)))
}

SP4_counts_nlm <- nls(log(SP4_counts) ~ log(counts_spline(age.at.S1K, basl, nu)),
                   data =  SP4_counts_df,
                   start = list(basl=6, nu=0.03))
SP4_counts_pars <- coef(SP4_counts_nlm)

SP4counts_fit <- data.frame(timeseq, "y_sp" = counts_spline(timeseq, SP4_counts_pars[1], SP4_counts_pars[2]))

ggplot() + 
  geom_point(data=SP4_counts_df, aes(x=age.at.S1K, y=SP4_counts), col=4, size =2) +
  geom_line(data = SP4counts_fit, aes(x = timeseq , y = y_sp), col=4, size =1) + 
  scale_y_log10(limits=c(1e5, 2e7)) +
  labs(title = 'Total counts of FoxP3 negative SP4 cells',  y=NULL,  x = 'Host age (days)') 

### Donor Ki67 fit
## phenomenological function
ki_spline <- function(Time, eps_0, eps_f){
  t0 = 49
  Ki_val = exp(- eps_f * (Time-49)) + eps_0
  return(Ki_val)
}

## LL function
SP4_donorKi_nlm <- nls(SP4_ki ~ ki_spline(age.at.S1K, eps_0, eps_f),
                   data =  SP4_donorki67_df,
                   start = list(eps_0=0.35, eps_f=0.02))
SP4_donorKi_pars <- coef(SP4_donorKi_nlm)

# prediction
SP4_donorKi_fit <- data.frame(timeseq, "y_sp" = ki_spline(timeseq, SP4_donorKi_pars[1], SP4_donorKi_pars[2]))
ggplot() + 
  geom_point(data=SP4_donorki67_df, aes(x=age.at.S1K, y=SP4_ki), col=4, size =2) +
  geom_line(data = SP4_donorKi_fit, aes(x = timeseq , y = y_sp), col=4, size =1) + 
  ylim(0,1) + labs(title = 'Ki67Hi fraction in donor-BM-derived FoxP3 negative SP4 cells',  y=NULL,  x = 'Days post BMT') 


ggplot() + 
  geom_point(data=SP4_hostki67_df, aes(x=age.at.S1K, y=SP4_ki), col=4, size =2) +
  geom_hline(yintercept = 0.326611, col=4, size =1) + 
  ylim(0,1) + labs(title = 'Ki67Hi fraction in host-BM-derived FoxP3 negative SP4 cells',  y=NULL,  x = 'Days post BMT') 




# exporting plot as PDF
pdf(file = file.path("out_fit", paste(colnames(cnts)[4], "IndPlots%03d.pdf", sep = "")),
    width = 12, height = 4.5, onefile = FALSE, useDingbats = FALSE )

cowplot::plot_grid(p2, p1, nrow = 1, labels = c("A", "B"))

dev.off()

