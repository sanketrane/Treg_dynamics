## clearing the environment
rm(list = ls())  
gc()    

library(rstan)
library(loo)
library(tidyverse)
library(bayesplot)
####################################################################################

## model specific details that needs to be change for every run
modelName <- "asm_deltavar_Ki2"

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
stanfit3 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_c3",".csv")))
stanfit4 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_c4",".csv")))

fit <- sflist2stanfit(list(stanfit1, stanfit4))

# finding the parameters used in the model 
# using the last parameter("sigma4") in the array to get the total number of parameters set in the model
num_pars <- which(fit@model_pars %in% "global_params") -1      # the variable "sigma4" will change depdending on the data used
parametersToPlot <- fit@model_pars[1:num_pars]

# number of post-burn-in samples that are used for plotting 
nPost <- nrow(fit)

################################################################################################
################################################################################################

## loading required datasets for plotting
counts_file <- file.path("data", "Counts_Treg.csv")
counts_data <- read.csv(counts_file) %>% 
  arrange(age.at.S1K) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 56, 'agebin1',
                             ifelse(age.at.BMT <= 70, 'agebin2',
                                    ifelse(age.at.BMT <= 84, 'agebin3', 'agebin4')))) %>%
  #gather(c(naive, memory), key='location', value = 'total_counts')
  select(-memory)

Nfd_file <- file.path("data", "Nfd_Treg.csv")
Nfd_data <- read.csv(Nfd_file) %>% 
  arrange(age.at.S1K)%>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 56, 'agebin1',
                             ifelse(age.at.BMT <= 70, 'agebin2',
                                    ifelse(age.at.BMT <= 84, 'agebin3', 'agebin4')))) 

hostki_file <- file.path("data", "hostKi67_Treg.csv")
hostki_data <- read.csv(hostki_file) %>% 
  arrange(age.at.S1K) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 56, 'agebin1',
                             ifelse(age.at.BMT <= 70, 'agebin2',
                                    ifelse(age.at.BMT <= 84, 'agebin3', 'agebin4'))),
         subcomp='Host') %>%
  #gather(c(naive, memory), key='location', value = 'prop_ki')
  select(-memory)

donorki_file <- file.path("data", "donorKi67_Treg.csv")
donorki_data <- read.csv(donorki_file) %>% 
  arrange(age.at.S1K) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 56, 'agebin1',
                             ifelse(age.at.BMT <= 70, 'agebin2',
                                    ifelse(age.at.BMT <= 84, 'agebin3', 'agebin4'))),
         subcomp='Donor')%>%
  #gather(c(naive, memory), key='location', value = 'prop_ki')
  select(-memory)


ki_data <- rbind(donorki_data, hostki_data)

# ################################################################################################
# # calculating PSIS-L00-CV for the fit
# naive_counts_loglik <- extract_log_lik(fit, parameter_name = "log_lik_counts_naive", merge_chains = TRUE)
# naive_fd_loglik <- extract_log_lik(fit, parameter_name = "log_lik_Nfd_naive", merge_chains = TRUE)
# ki_donor_naive_loglik <- extract_log_lik(fit, parameter_name = "log_lik_ki_donor_naive", merge_chains = TRUE)
# ki_host_naive_loglik <- extract_log_lik(fit, parameter_name = "log_lik_ki_host_naive", merge_chains = TRUE)
# 
# #combined_loglik <- extract_log_lik(fit, parameter_name = "log_lik", merge_chains = TRUE)
# log_lik_comb <- cbind(naive_counts_loglik, naive_fd_loglik,
#                       ki_donor_naive_loglik, ki_host_naive_loglik)
# 
# 
# # optional but recommended
# ll_array <- extract_log_lik(fit,parameter_name = "log_lik_counts_naive", merge_chains = FALSE)
# r_eff <- relative_eff(exp(ll_array))
# 
# # loo-ic values
# loo_loglik <- loo(log_lik_comb, save_psis = FALSE, cores = 8)
# 
# # Widely applicable AIC
# AICw_lok <- waic(cbind(naive_counts_loglik, naive_fd_loglik, 
#                       ki_donor_naive_loglik, ki_host_naive_loglik))
# 
# # AIC from LLmax
# #AIC_lok <-  -2 * max(combined_loglik)  + 2 * length(parametersToPlot)
# ploocv <- data.frame("Model" = modelName,
#                      "LooIC" = loo_loglik$estimates[3],
#                      "SE" = loo_loglik$estimates[6], 
#                      "PLoo" = loo_loglik$estimates[2])
# ploocv
# 
# write.table(ploocv, file = file.path(outputDir, "stat_table.csv"),
#             sep = ",", append = F, quote = FALSE,
#             col.names = T, row.names = FALSE)

### posterior distributions of parameters
ptable <- monitor(as.array(fit, pars = parametersToPlot), warmup = 0, print = FALSE)
out_table <- ptable[1:num_pars, c(1, 3, 4, 8)]
out_table
write.csv(out_table, file = file.path(outputDir, paste0('params_', modelName, ".csv")))



################################################################################################
################################################################################################
## posterior predictive distributions

source('scripts/stan_extract_forplotting.R')

legn_labels <- c('6-8', '8-10', '10-12', '12-25')

ggplot() +
  geom_ribbon(data = Counts_pred, aes(x = timeseries, ymin = lb, ymax = ub, fill = ageBMT_bin), alpha = 0.2)+
  geom_ribbon(data = Counts_pred, aes(x = timeseries, ymin = lb, ymax = ub, fill = ageBMT_bin), alpha = 0.2)+
  geom_ribbon(data = Counts_naive_withsigma, aes(x = timeseries, ymin = lb, ymax = ub, fill = ageBMT_bin), alpha = 0.1)+
  geom_line(data = Counts_pred, aes(x = timeseries, y = median, color = ageBMT_bin)) +
  #geom_errorbar(data = Counts_naive_sigma_obs, aes(x = timeseries, ymin = lb, ymax = ub, col=ageBMT_bin),
  #              alpha = 0.25, width=0.02)+
  #geom_errorbar(data = Counts_per_sigma_obs, aes(x = timeseries, ymin = lb, ymax = ub, col=ageBMT_bin),
  #              alpha = 0.15, width=0.02)+
  geom_point(data = counts_data, aes(x = age.at.S1K, y = naive, color = ageBMT_bin), size=2) +
  labs(title=paste('Total counts of naive Tregs'),  y=NULL, x= "Host age (days)") + 
  scale_color_discrete(name="Host age at \n BMT (Wks)", labels=legn_labels)+
  scale_x_continuous(limits = c(60, 450) , trans="log10", breaks=c(10, 30, 100, 300)) + #scale_y_log10() +
  scale_y_continuous(limits = c(5e4, 2e7), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  #facet_wrap(~ factor(location, levels = c('naive', "memory")))+
  guides(fill = 'none') + myTheme 

ggsave(filename = file.path(outputDir, paste0(modelName, "P1.pdf")), last_plot(),
       device = "pdf", height = 4.5, width = 6)


# normalised donor fractions

ggplot() +
  geom_ribbon(data = Nfd_pred, aes(x = timeseries, ymin = lb, ymax = ub, fill = ageBMT_bin), alpha = 0.15)+
  geom_line(data = Nfd_pred, aes(x = timeseries, y = median, color = ageBMT_bin)) +
  geom_point(data = Nfd_data, aes(x = age.at.S1K, y = naive_DP1, color = ageBMT_bin), size=2) +
  labs(x = "Host age (days)", y = NULL, title = "Normalised Chimerism in naive Tregs") +
  scale_color_discrete(name="Host age at \n BMT (Wks)", labels=legn_labels)+
  scale_x_continuous(limits = c(1, 450), breaks = c(0,100,200,300, 400, 500))+
  scale_y_continuous(limits =c(0, 1.02), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) + 
 # facet_wrap(~ factor(location, levels = c('naive', "memory")))+
  guides(fill='none')+ myTheme


ggsave(filename = file.path(outputDir, paste0(modelName, "P2.pdf")), last_plot(),
       device = "pdf", height = 4.5, width = 6)

# thymic Ki67 fractions

fac_labels <- c(`agebin1`= '6-8 weeks', `agebin2`= '8-10 weeks', `agebin3`= '10-12 weeks', `agebin4`= '12-25 weeks')

ggplot() +
  geom_ribbon(data = ki_naive_pred, aes(x = timeseries, ymin = lb*100, ymax = ub*100, fill = subcomp), alpha = 0.15)+
  geom_line(data = ki_naive_pred, aes(x = timeseries, y = median*100, color = subcomp)) +
  geom_point(data = ki_data, aes(x = age.at.S1K, y = naive*100, color = subcomp), size=1.5) +
  labs(x = "Host age (days)", y = NULL, title = "% Ki67hi in naive Tregs") +
  scale_x_continuous(limits = c(60, 450), breaks = c(0,100,200,300, 400, 500))+
  scale_y_continuous(limits =c(0, 50), breaks = c(0, 10, 20, 30, 40, 50))+ 
  facet_wrap(~ ageBMT_bin, scales = 'free', labeller = as_labeller(fac_labels))+
  guides(fill='none') + myTheme + theme(legend.title = element_blank())


ggsave(filename = file.path(outputDir, paste0(modelName, "P3.pdf")), last_plot(),
       device = "pdf", height = 6, width = 8.5)


################################################################################################
## open graphics device 
## saving  plots for quality control 
pdf(file = file.path(outputDir, paste(modelName,"Plots%03d.pdf", sep = "")),
    width = 12, height = 5, onefile = F)

pairs(fit, pars = parametersToPlot)

options(bayesplot.base_size = 15,
        bayesplot.base_family = "sans")
bayesplot::color_scheme_set(scheme = "viridis")

rhats <- rhat(fit, pars = parametersToPlot)
mcmc_rhat(rhats) + yaxis_text() + myTheme

ratios1 <- neff_ratio(fit, pars = parametersToPlot)
mcmc_neff(ratios1) + yaxis_text() + myTheme

posterior <- as.array(fit)
mcmc_acf(posterior, pars = parametersToPlot) + myTheme

mcmcHistory(fit, pars = parametersToPlot, nParPerPage = 4, myTheme = myTheme)

mcmc_dens_overlay(posterior, parametersToPlot)
mcmc_dens(posterior, parametersToPlot) + myTheme


dev.off()
