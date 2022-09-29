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
SP4_counts <- read_csv("data/Counts_thymicSource.csv") %>%
  rename(SP4_counts = FoxP3_Neg_SP4) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 63, 'agegroup1',
                             ifelse(age.at.BMT <= 88, 'agegroup2',
                                    'agegroup3')))
SP4_chimerism <- read_csv("data/Chimerism_thymicSource.csv") %>%
  rename(SP4_chi = FoxP3_Neg_SP4) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 63, 'agegroup1',
                             ifelse(age.at.BMT <= 88, 'agegroup2',
                                    'agegroup3')))
SP4_donorki67 <- read_csv("data/donorKi67_thymicSource.csv") %>%
  rename(SP4_ki = FoxP3_Neg_SP4) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 63, 'agegroup1',
                             ifelse(age.at.BMT <= 88, 'agegroup2',
                                    'agegroup3')))
SP4_hostki67 <- read_csv("data/hostKi67_thymicSource.csv") %>%
  rename(SP4_ki = FoxP3_Neg_SP4) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 63, 'agegroup1',
                             ifelse(age.at.BMT <= 88, 'agegroup2',
                                    'agegroup3')))


qplot(x=age.at.S1K, y = SP4_ki, data = SP4_donorki67, col=ageBMT_bin) + ylim(0,1)
  scale_y_log10(limits=c(0.1, 1e0))

### spline functions
count_spline <- function(age.at.S1K, basl, nu){
  Time = age.at.S1K   ### t0 = 66 days
  return(10^basl * exp(-nu * Time))
}

ki_spline <- function(Time, eps_0, q, n){
  dpt=Time-66
  Ki_val = #exp(- eps_f * (Time + A)) + eps_0
    eps_0 * (1 - ((dpt^q)/(1 + (dpt^n))))
  return(Ki_val)
}

timeseq <- seq(60, 450)
qplot(x=timeseq, y=ki_spline(timeseq, 0.7, 0.5, 4))




## LL function for counts
count_obj <- function(params) {
  basl <- params[["basl"]]
  nu  <- params[["nu"]]
  
  ## SSR
  SP4_counts %>%
    mutate(count_fit = count_spline(age.at.S1K,
                                    basl = basl,
                                    nu = nu),
           log_sr = (log(SP4_counts) - log(count_fit))^2) %>%
    summarize(log_ssr = sum(log_sr)) %>%
    unlist()
}

# initial guess for LL calculation
params <- list(basl = 6, nu=0.03)

## model fit
counts_fit <- optim(params, count_obj, control = list(trace = TRUE, maxit = 1000))
counts_pars_est <- counts_fit$par


logit_transf <- function(x){log(x/(1-x))}

## LL function for ki proportions
ki_obj <- function(params, dataf) {
  eps_0 <- params[["eps_0"]]
  eps_f  <- params[["eps_f"]]
  A  <- params[["A"]]
  
  ## SSR
  dataf %>%
    mutate(ki_fit = ki_spline(age.at.S1K,
                              eps_0 = eps_0,
                              eps_f = eps_f, 
                              A),
           SqRes = (logit_transf(SP4_ki) - logit_transf(ki_fit))^2) %>%
    summarize(SSR = sum(SqRes)) %>%
    unlist()
}

# initial guess for LL calculation
params <- list(eps_0 = 0.3, eps_f = 0.01, A=10)

## model fit
ki_fit_donor <- optim(params, ki_obj, dataf = SP4_donorki67,
                      control = list(trace = TRUE, maxit = 1000))
ki_donor_est <- ki_fit_donor$par

## time course for prediction
ts_p <- seq(66, 450)

## prediction plots
ki_fit_donor <- data.frame(ts_p, "y_sp" = ki_spline(ts_p, ki_donor_est[1], ki_donor_est[2],  ki_donor_est[3]))

ggplot() + 
  geom_point(data = SP4_donorki67, aes(age.at.S1K, SP4_ki),  col=4, size =2.5) + 
  geom_line(data = ki_fit_donor, aes(x = ts_p, y = y_sp), col = "#feb236", size =1) + 
  scale_y_continuous(limits = c(0, 1)) +
  labs(y =NULL, x = 'Mouse age (days)', title = 'Proportions of Ki67+ cells within SP8 compartment') +
  theme_bw()

## prediction plots
#counts_fit_spline <- data.frame(ts_p, "y_sp" = count_spline(ts_p, counts_pars_est[1], counts_pars_est[2], counts_pars_est[3], counts_pars_est[4], counts_pars_est[5]))

counts_fit_spline <- data.frame(ts_p, "y_sp" = count_spline(ts_p, log(4.3e5), 1.8e3, 2.1, 30, 3.7))

p2 <- ggplot(data_df) + 
  geom_point(aes(days, counts), col=4, size =2.5) + 
  #geom_vline(xintercept = 49, col = 'darkred', linetype=2)+
  #geom_line(data = counts_fit_spline, aes(x = ts_p, y = y_sp), col = "#feb236", size =1) + 
  scale_x_log10(limits= c(3, 350), breaks = c(3,10,30,100, 300)) + 
  scale_y_continuous(limits = c(5e4, 2e6), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  labs(y =NULL, x = 'Mouse age (days)', title = 'Counts of SP8 T cells') +
  myTheme


# exporting plot as PDF
pdf(file = file.path("out_fit", paste(colnames(cnts)[4], "IndPlots%03d.pdf", sep = "")),
    width = 12, height = 4.5, onefile = FALSE, useDingbats = FALSE )

cowplot::plot_grid(p2, p1, nrow = 1, labels = c("A", "B"))

dev.off()

