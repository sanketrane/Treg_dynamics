## Simple linear ODE model -- Homogeneous model 
## clearing the environment
rm(list = ls()); gc()  
## loading required libraries
require(cmdstanr)
require(parallel)
require(loo)
require(tidyverse)

# model specific details
modelName <- "Neutral_SP4source"

# variables for data files
data_derived1 <- "counts_THnaiTreg.csv"
data_derived2 <- "Nfd_THnaiTreg.csv"
data_derived3 <- "counts_naiTreg.csv"
data_derived4 <- "Nfd_naiTreg.csv"

## Setting all the directories for opeartions
projectDir <- getwd()
scriptDir <- file.path(projectDir, "scripts")
modelDir <- file.path(projectDir, "models")
dataDir <- file.path(projectDir, "data")
outputDir <- file.path(projectDir, "output")
saveDir <- file.path(outputDir, modelName)

###############################################################################################
################################################################################################
NTregThy_counts <- read_csv(file.path(dataDir, data_derived1))%>% arrange(age.at.S1K)
NTregThy_Nfd <- read_csv(file.path(dataDir, data_derived2))%>% arrange(age.at.S1K)
NTregPer_counts <- read_csv(file.path(dataDir, data_derived3))%>% arrange(age.at.S1K)
NTregPer_Nfd <- read_csv(file.path(dataDir, data_derived4))%>% arrange(age.at.S1K)

## time-points in data
data_time <- NTregPer_Nfd$age.at.S1K
## ages of animals at their respective BMTs
ageAtBMT <- NTregPer_Nfd$age.at.BMT
## tb_time -- unique ageAtBMT for solving initial conditions 
tb_time <- ageAtBMT %>% unique() %>% sort()
##keep track of index of tb_time points in relation to data_AgeBMT
tb_index <- purrr::map_dbl(ageAtBMT, function(x) which(x == tb_time))

## In this analysis mice were binned in 3 groups dependeing on the age at BMT
## mean ages for each group were calculated and are --> 49, 72, 128
## predictions are made for each bin from the lowest timepoint in that group to d450.
ts_pred1 = seq(from = 66, to = 450, length.out = 300)
ts_pred2 = seq(from = 90, to = 450, length.out = 300)
ts_pred3 = seq(from = 180, to = 450, length.out = 300)


## create data set
data_list <- list(
  num_obs = nrow(NTregPer_Nfd),      ## number of observations
  num_tb = length(tb_time),          ## number of unique age at BMTs
  num_index = length(tb_index),      ## indexes of unique ageATBMT matching up with all ageAtBMTs
  data_time = data_time,
  ageAtBMT = ageAtBMT,
  tb_time = tb_time,
  tb_index = tb_index,               
  thy_counts = NTregThy_counts$total_counts,  
  thy_Nfd = NTregThy_Nfd$Nfd,                 
  per_counts = NTregPer_counts$total_counts,  
  per_Nfd = NTregPer_Nfd$Nfd,
  num_pred = 300,
  ts_pred1 = ts_pred1,
  ts_pred2 = ts_pred2,
  ts_pred3 = ts_pred3
  )


## create initial estimates
init_list <- function() list(
  alpha1 = exp(rnorm(1, log(0.5), 0.2)),
  alpha2 = exp(rnorm(1, log(0.5), 0.2)),
  mu = exp(rnorm(1,log(0.3), 1)),
  beta = exp(rnorm(1,log(0.1), 1)),
  lambda = exp(rnorm(1,log(0.05), 1)),
  
  thy_N0Log = rnorm(1, 10 , 0.5),
  per_N0Log = rnorm(1, 14 , 0.5),
  
  sigma1 = exp(rnorm(1,log(1.5), 1)),
  sigma2 = exp(rnorm(1,log(1.5), 1)),
  sigma3 = exp(rnorm(1,log(1.5), 1)),
  sigma4 = exp(rnorm(1,log(1.5), 1)))

## Specify the variables for which you want history and density plots
parametersToPlot <- c("alpha1", "alpha2", "mu", "beta", "lambda", "thy_N0Log",
                      "per_N0Log", "sigma1", "sigma2",  "sigma3", "sigma4")

## Additional variables to monitor
otherRVs <- c("log_lik1", "log_lik2", "log_lik3", "log_lik4",
              "y1_mean_pred1", "y2_mean_pred1", "y3_mean_pred1", "y4_mean_pred1",
              "y1_mean_pred2", "y2_mean_pred2", "y3_mean_pred2", "y4_mean_pred2",
              "y1_mean_pred3", "y2_mean_pred3", "y3_mean_pred3", "y4_mean_pred3")

parameters <- c(parametersToPlot, otherRVs)

################################################################################################
## run Stan
#### compile model
model_file <- file.path(modelDir, 'Neutral_SP4Source.stan')
mod <- cmdstan_model(model_file)

### sampling params to fit model to data
stan_fit <- mod$sample(data = data_list, init = init_list, 
                       iter_warmup = 300, iter_sampling = 300,
                       chains = 2, parallel_chains = 2,
                       save_warmup = T, refresh = 100, 
                       adapt_delta = 0.9)

# saving the stan fit object 
## dir to save output
if (!file.exists(saveDir)){
  dir.create(saveDir)
}
stanfit_save <- rstan::read_stan_csv(stan_fit$output_files())
output_filename <- file.path(saveDir, paste0(modelName, ".rds"))
write_rds(stanfit_save, file = file.path(output_filename))


### parameters table
num_pars <- length(parametersToPlot)
ptable <- monitor(as.array(stanfit_save, pars = parameters_to_plot), warmup = 0, print = FALSE)
out_table <- ptable[1:num_pars, c(1, 3, 4, 8)]
write.csv(out_table, file = file.path(saveDir, paste0('params_', modelName, ".csv")))


# loo-ic values
loo_loglik1 <- loo::extract_log_lik(stanfit_save, parameter_name = "log_lik1", merge_chains = TRUE)
loo_loglik2 <- loo::extract_log_lik(stanfit_save, parameter_name = "log_lik2", merge_chains = TRUE)
loo_loglik3 <- loo::extract_log_lik(stanfit_save, parameter_name = "log_lik3", merge_chains = TRUE)
loo_loglik4 <- loo::extract_log_lik(stanfit_save, parameter_name = "log_lik4", merge_chains = TRUE)
combined_loglik <- cbind(loo_loglik1, loo_loglik2, loo_loglik3, loo_loglik4)

loo_ic <- loo::loo(combined_loglik,  save_psis = FALSE, cores = 4)
ploocv <- loo_ic$estimates
write.csv(ploocv, file = file.path(saveDir, paste0('stats_', modelName, ".csv")))


print("DONE!!!")

