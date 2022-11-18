## clearing the environment
rm(list = ls())  
gc()    

library(rstan)
library(loo)
library(tidyverse)
library(bayesplot)
####################################################################################

## model specific details that needs to be change for every run
modelName <- "rtem"

## Setting all the directories for opeartions
projectDir <- getwd()
scriptDir <- file.path(projectDir, "scripts")
modelDir <- file.path(projectDir, "stan_models")
dataDir <- file.path(projectDir, "datafiles")
toolsDir <- file.path(scriptDir, "tools")
outputDir <- file.path(projectDir, "output_fit")
saveDir <- file.path(projectDir, 'save_csv')
LooDir <- file.path('loo_fit') 

# loadiong the scr# loadiong the script that contains functions for plotting stan parameters
source(file.path(toolsDir, "stanTools.R"))                # save results in new folder

# compiling multiple stan objects together that ran on different nodes
stanfit1 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_c1", ".csv")))
stanfit2 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_c2",".csv")))

fit <- sflist2stanfit(list(stanfit1, stanfit2))

# finding the parameters used in the model 
# using the last parameter("sigma4") in the array to get the total number of parameters set in the model
num_pars <- which(fit@model_pars %in% "sigma_host_ki_thy")      # the variable "sigma4" will change depdending on the data used
parametersToPlot <- fit@model_pars[1:num_pars]

# number of post-burnin samples that are used for plotting 
nPost <- nrow(fit)

################################################################################################
################################################################################################

## loading required datasets for plotting
counts_file <- file.path("data", "Counts_naiTreg.csv")
counts_data <- read.csv(counts_file) %>% 
  arrange(age.at.S1K) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 56, 'agebin1',
                             ifelse(age.at.BMT <= 70, 'agebin2',
                                    ifelse(age.at.BMT <= 84, 'agebin3', 'agebin4')))) 

Nfd_file <- file.path("data", "Nfd_naiTreg.csv")
Nfd_data <- read.csv(Nfd_file) %>% 
  arrange(age.at.S1K)%>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 56, 'agebin1',
                             ifelse(age.at.BMT <= 70, 'agebin2',
                                    ifelse(age.at.BMT <= 84, 'agebin3', 'agebin4'))))

hostki_file <- file.path("data", "hostKi67_naiTreg.csv")
hostki_data <- read.csv(hostki_file) %>% 
  arrange(age.at.S1K) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 56, 'agebin1',
                             ifelse(age.at.BMT <= 70, 'agebin2',
                                    ifelse(age.at.BMT <= 84, 'agebin3', 'agebin4'))))

donorki_file <- file.path("data", "donorKi67_naiTreg.csv")
donorki_data <- read.csv(donorki_file) %>% 
  arrange(age.at.S1K) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 56, 'agebin1',
                             ifelse(age.at.BMT <= 70, 'agebin2',
                                    ifelse(age.at.BMT <= 84, 'agebin3', 'agebin4'))))

# ################################################################################################
# calculating PSIS-L00-CV for the fit
thycounts_loglik <- extract_log_lik(fit, parameter_name = "log_lik_counts_thy", merge_chains = TRUE)
thyfd_loglik <- extract_log_lik(fit, parameter_name = "log_lik_Nfd_thy", merge_chains = TRUE)
percounts_loglik <- extract_log_lik(fit, parameter_name = "log_lik_counts_per", merge_chains = TRUE)
perfd_loglik <- extract_log_lik(fit, parameter_name = "log_lik_Nfd_per", merge_chains = TRUE)
ki_donor_thy_loglik <- extract_log_lik(fit, parameter_name = "log_lik_ki_donor_thy", merge_chains = TRUE)
ki_host_thy_loglik <- extract_log_lik(fit, parameter_name = "log_lik_ki_host_thy", merge_chains = TRUE)
ki_donor_per_loglik <- extract_log_lik(fit, parameter_name = "log_lik_ki_donor_per", merge_chains = TRUE)
ki_host_per_loglik <- extract_log_lik(fit, parameter_name = "log_lik_ki_host_per", merge_chains = TRUE)

#combined_loglik <- extract_log_lik(fit, parameter_name = "log_lik", merge_chains = TRUE)
log_lik_comb <- cbind(thycounts_loglik, thyfd_loglik, percounts_loglik, perfd_loglik,
                      ki_donor_thy_loglik, ki_host_thy_loglik, ki_donor_per_loglik, ki_host_per_loglik)


# optional but recommended
ll_array <- extract_log_lik(fit,parameter_name = "log_lik_counts_thy", merge_chains = FALSE)
r_eff <- relative_eff(exp(ll_array))

# loo-ic values
loo_loglik <- loo(log_lik_comb, save_psis = FALSE, cores = 8)

# Widely applicable AIC
AICw_lok <- waic(cbind(thycounts_loglik, thyfd_loglik, percounts_loglik, perfd_loglik,
                       ki_donor_thy_loglik, ki_host_thy_loglik, ki_donor_per_loglik, ki_host_per_loglik))

# AIC from LLmax
#AIC_lok <-  -2 * max(combined_loglik)  + 2 * length(parametersToPlot)
ploocv <- data.frame("Model" = modelName,
                     "LooIC" = loo_loglik$estimates[3],
                     "SE" = loo_loglik$estimates[6], 
                     "PLoo" = loo_loglik$estimates[2])
ploocv

write.table(ploocv, file = file.path(outputDir, "stat_table.csv"),
            sep = ",", append = T, quote = FALSE,
            col.names = F, row.names = FALSE)

### posterior distributions of parameters
ptable <- monitor(as.array(fit, pars = parametersToPlot), warmup = 0, print = FALSE)
out_table <- ptable[1:num_pars, c(1, 3, 4, 8)]
out_table
write.csv(out_table, file = file.path(outputDir, paste0('params_', modelName, ".csv")))



################################################################################################
################################################################################################
## posterior predictive distributions

ggplot() +
  geom_ribbon(data = Counts_thy_pred, aes(x = timeseries, ymin = lb, ymax = ub, fill = ageBMT_bin), alpha = 0.2)+
  geom_line(data = Counts_thy_pred, aes(x = timeseries, y = median, color = ageBMT_bin), size=1.2) +
  geom_point(data = counts_data, aes(x = age.at.S1K, y = Thymus, color = ageBMT_bin), size=2) +
  labs(title=paste('Total counts of thymic naive Tregs'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(60, 450) , trans="log10", breaks=c(10, 30, 100, 300))+
  scale_y_continuous(limits = c(5e3, 5e5), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  guides(fill = 'none') + myTheme


ggplot() +
  geom_ribbon(data = Counts_per_pred, aes(x = timeseries, ymin = lb, ymax = ub, fill = ageBMT_bin), alpha = 0.2)+
  geom_line(data = Counts_per_pred, aes(x = timeseries, y = median, color = ageBMT_bin), size=1.2) +
  geom_point(data = counts_data, aes(x = age.at.S1K, y = Periphery, color = ageBMT_bin), size=2) +
  labs(title=paste('Total counts of peripheral naive Tregs'),  y=NULL, x= "Host age (days)") + 
  scale_x_continuous(limits = c(60, 450) , trans="log10", breaks=c(10, 30, 100, 300))+
  scale_y_continuous(limits = c(1e5, 1e7), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  guides(fill = 'none') + myTheme

# normalised donr fractions

ggplot() +
  geom_hline(aes(yintercept = 1), color = "#d11100", linetype = 2, size=1.2)+
  geom_ribbon(data = Y2pred, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#0099cc", alpha = 0.2)+
  geom_line(data = Y2pred, aes(x = timeseries, y = median), color = "#1e2366", size=1.2) +
  geom_point(data = NTregThy_Nfd, aes(x = time.post.BMT, y = Nfd), color = "#1e2366", size=2) +
  labs(x = "Days post BMT", y = NULL, title = "Chimerism in thymic naive T regs normalised to chimersim in thymic SP4 T cells") +
  scale_x_continuous(limits = c(0, 300), breaks = c(0,100,200,300))+
  scale_y_continuous(limits =c(0, 1.02), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0))+ 
  guides(color = FALSE)+ myTheme


Y3pred <- as.data.frame(fit, pars = "y3_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)

Y3pred <- Y3pred%>%
  filter(timeseries >= 10)%>% filter(timeseries <= 300)



ggplot() +
  # geom_ribbon(data = Cpred, aes(x = timeseries, ymin = lb, ymax=ub), fill = '#0099cc', alpha = 0.2) +
  geom_ribbon(data = Y3pred, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#0099cc", alpha = 0.2)+
  geom_line(data = Y3pred, aes(x = timeseries, y = median), color = "#1e2366", size=1.2) +
  geom_point(data = NTregPer_counts, aes(x = time.post.BMT, y = total_counts), col = 'darkblue', size=2) +
  labs(title=paste('Total counts of thymic naive Tregs'),  y=NULL, x= "Days post BMT") + 
  scale_x_continuous(limits = c(10, 300) , trans="log10", breaks=c(10, 30, 100, 300))+
  scale_y_continuous(limits = c(1e5, 5e6), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  guides(color = FALSE) + myTheme


# normalised donr fractions
Y4pred <- as.data.frame(fit, pars = "y4_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)

Y4pred <- Y4pred%>%
  filter(timeseries >= 0)%>% filter(timeseries <= 300)


ggplot() +
  geom_hline(aes(yintercept = 1), color = "#d11100", linetype = 2, size=1.2)+
  geom_ribbon(data = Y4pred, aes(x = timeseries, ymin = lb, ymax = ub), fill = "#0099cc", alpha = 0.2)+
  geom_line(data = Y4pred, aes(x = timeseries, y = median), color = "#1e2366", size=1.2) +
  geom_point(data = NTregPer_Nfd, aes(x = time.post.BMT, y = Nfd), color = "#1e2366", size=2) +
  labs(x = "Days post BMT", y = NULL, title = "Chimerism in peripheral naive T regs normalised to chimersim in thymic SP4 T cells") +
  scale_x_continuous(limits = c(0, 300), breaks = c(0,100,200,300))+
  scale_y_continuous(limits =c(0, 1.02), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0))+ 
  guides(color = FALSE)+ myTheme


dev.off()




################################################################################################
## open graphics device 
## saving  plots for quality control 
pdf(file = file.path(outputDir, paste(modelName,"Plots%03d.pdf", sep = "")),
    width = 8, height = 5, onefile = F)

pairs(fit, pars = parametersToPlot)

options(bayesplot.base_size = 15,
        bayesplot.base_family = "sans")
bayesplot::color_scheme_set(scheme = "viridis")
myTheme <- theme(text = element_text(size = 12), axis.text = element_text(size = 12), axis.title =  element_text(size = 12, face = "bold"),
                 plot.title = element_text(size=12, face = 'bold',  hjust = 0.5), legend.text = element_text(size=12),
                 legend.title = element_text(size = 12))

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

rhats <- rhat(fit, pars = parametersToPlot)
mcmc_rhat(rhats) + yaxis_text() + myTheme

ratios1 <- neff_ratio(fit, pars = parametersToPlot)
mcmc_neff(ratios1) + yaxis_text() + myTheme

posterior <- as.array(fit)
mcmc_acf(posterior, pars = parametersToPlot) + myTheme

mcmcHistory(fit, pars = parametersToPlot, nParPerPage = 4, myTheme = myTheme)

mcmc_dens_overlay(posterior, parametersToPlot)
mcmc_dens(posterior, parametersToPlot) + myTheme
