## Simple linear ODE model -- Homogeneous model 
## clearing the environment
rm(list = ls())  
gc()    
#setwd("Desktop/GIt_repos/Treg_dynamics")

####################################################################################
## Installing r-stan pachage on the go:
if(!("rstan" %in% rownames(installed.packages())) ){
  install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies=TRUE)
}

## Installing loo:
if(!("loo" %in% rownames(installed.packages())) ){
  install.packages("loo", repos = "https://cloud.r-project.org/", dependencies=TRUE)
}

## Installing tidyverse:
if(!("tidyverse" %in% rownames(installed.packages())) ){
  install.packages("tidyverse", repos = "https://cloud.r-project.org/", dependencies=TRUE)
}
####################################################################################

require(rstan)
require(parallel)
require(loo)
require(tidyverse)

# model specific details
modelName <- "SH_model"
source_pop <- "DP1"

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
toolsDir <- file.path(scriptDir, "tools")
outputDir <- file.path(projectDir, "output")
saveDir <- file.path(outputDir, paste(modelName, "_", substr(data_derived1, 1,2), sep=""))

# loadiong the scr# loadiong the script that contains functions for plotting stan parameters
source(file.path(toolsDir, "stanTools.R"))

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# set different seeds for different tasks
i <- as.numeric(Sys.getenv("SLURM_PROCID"))
seed <- 4010 + i
################################################################################################
################################################################################################
NTregThy_counts <- read_csv(file.path(dataDir, data_derived1))%>% arrange(age.at.S1K)
NTregThy_Nfd <- read_csv(file.path(dataDir, data_derived2))%>% arrange(age.at.S1K)
NTregPer_counts <- read_csv(file.path(dataDir, data_derived3))%>% arrange(age.at.S1K)
NTregPer_Nfd <- read_csv(file.path(dataDir, data_derived4))%>% arrange(age.at.S1K)

## solving ODEs only for unique timepoints
## stan ode solver throws an erro when time points are repeated
unique_times_df <- NTregPer_Nfd %>% distinct(age.at.S1K, .keep_all = TRUE)

## actual timepoints in data
data_time <- NTregPer_Nfd$age.at.S1K
## unique timepoints for which ODEs will be solved
solve_time <- unique_times_df$age.at.S1K  
#keep track of index of solve_time points in relation to data_time
time_index <- purrr::map_dbl(data_time, function(x) which(x == solve_time))
  
## delay time for age correction
## For each host with different age at BMT -- time zero is different which is another varibale in the ODEs
## tb_time -- time at BMT -- for solving initial conditions at each tb we use the earliest host age at BMT to calculate the initial conidtions at each tb
data_AgeBMT <- unique_times_df$age.at.BMT
tb_time <- data_AgeBMT %>% unique() %>% sort()
##keep track of index of tb_time points in relation to data_AgeBMT
tb_index <- purrr::map_dbl(data_AgeBMT, function(x) which(x == tb_time))

## time sequence for predictions specific to age bins within the data
## In this analysis 3 age bins were selected with mean ages for each bin calculated as follows:
Age_BMT_sorted <- NTregPer_Nfd %>% arrange(age.at.BMT)
Mean_agebin1 <- round(mean(Age_BMT_sorted$age.at.BMT[1:11]))
Mean_agebin2 <- round(mean(Age_BMT_sorted$age.at.BMT[12:34]))
Mean_agebin3 <- round(mean(Age_BMT_sorted$age.at.BMT[35:43]))

## predictions will be made for following time points within each bin
ts_pred1 = seq(from = Mean_agebin1, to = 450, by = 1)
ts_pred2 = seq(from = Mean_agebin2, to = 450, by = 1)
ts_pred3 = seq(from = Mean_agebin3, to = 450, by = 1)

## predictions will be for following host ages at which BM was transferred within each bin
tb_pred1 <- rep(Mean_agebin1, length(ts_pred1))
tb_pred2 <- rep(Mean_agebin2, length(ts_pred2))
tb_pred3 <- rep(Mean_agebin3, length(ts_pred3))

## in order to get inital conditions for each bin the vcator of prediction tb_time points
tb_time_pred1 <- c(min(Age_BMT_sorted), Mean_agebin1)
tb_time_pred2 <- c(min(Age_BMT_sorted), Mean_agebin2)
tb_time_pred3 <- c(min(Age_BMT_sorted), Mean_agebin3)


### data inputs depdending on the precursor population used
## these value dictate how splines for source counts and chimerism change with time
if (grepl("SP4", source_pop) == TRUE){
  source_list = data.frame("nu" = 0.00396, "theta0" = 15.29, "chiEst" = 0.82, "qEst" = 0.047)
} else {
  ## deafult precursor pop is Thymic DP1
  source_list = data.frame("nu" = 0.00293, "theta0" = 18.06, "chiEst" = 0.84, "qEst" = 0.057)
} 

## create data set
data <- list(
  numObs = nrow(NTregPer_Nfd),      ## total number of observations within data
  num_index = length(solve_time),   ## to decalre dimensions of solve_time array in stan
  solve_time = solve_time,          ## unique time points for ODE solver
  time_index = time_index,
  num_tb = length(tb_time),         ## to decalre dimensions of tb_time array in stan
  tb_time = tb_time,
  tb_index = tb_index,              ## dimensions of tb_index are as same as solve_time
  ageBMT = data_AgeBMT,
  dpBMT = NTregPer_Nfd$time.post.BMT,  ## for spline of chimerism of the source pop
  thy_counts = NTregThy_counts$total_counts,  ## data to fit
  thy_Nfd = NTregThy_Nfd$Nfd,                 ## data to fit
  per_counts = NTregPer_counts$total_counts,  ## data to fit
  per_Nfd = NTregPer_Nfd$Nfd,                 ## data to fit
  numPred1 = length(ts_pred1),      ## to decalre dimensions of ts_pred1 array in stan
  numPred2 = length(ts_pred2),      ## to decalre dimensions of ts_pred2 array in stan
  numPred3 = length(ts_pred3),      ## to decalre dimensions of ts_pred3 array in stan
  ts_pred1 = ts_pred1,
  ts_pred2 = ts_pred2,
  ts_pred3 = ts_pred3,
  tb_pred1 = tb_pred1,
  tb_pred2 = tb_pred2,
  tb_pred3 = tb_pred3,
  tb_time_pred1 = tb_time_pred1,
  tb_time_pred2 = tb_time_pred2,
  tb_time_pred3 = tb_time_pred3,
  thy_Nd0 = 0,                      ## initial condition for donor pool in thymus
  per_Nd0 = 0,                      ## initial condition for donor pool in periphery
  theta0 = exp(source_list$theta0),
  nu = source_list$nu,
  chiEst = source_list$chiEst,
  qEst = source_list$qEst
)


## create initial estimates
init <- function() list(
  psi = exp(rnorm(1, log(0.01), 0.2)),
  alpha = exp(rnorm(1, log(0.01), 0.2)),
  lambda_thy = exp(rnorm(1,log(0.02), 1)),
  lambda_per = exp(rnorm(1,log(0.02), 1)),
  
  thy_y0Log = rnorm(1, 13 , 0.1),
  thy_y0Log = rnorm(1, 11 , 0.1),
  
  sigma1 = exp(rnorm(1,log(1.5), 1)),
  sigma2 = exp(rnorm(1,log(1.5), 1)),
  sigma3 = exp(rnorm(1,log(1.5), 1)),
  sigma4 = exp(rnorm(1,log(1.5), 1)))

## Specify the variables for which you want history and density plots
parametersToPlot <- c("psi", "alpha", "lambda_thy", "lambda_per", "thy_y0Log", "sigma1", "sigma2",  "sigma3", "sigma4")

## Additional variables to monitor
otherRVs <- c("y1_mean_pred_age1", "y1_mean_pred_age2", "y1_mean_pred_age3", "thy_countspred_age1", "thy_countspred_age2", "thy_countspred_age3", 
              "y2_mean_pred_age1", "y2_mean_pred_age2", "y2_mean_pred_age3", "thy_fdpred_age1", "thy_fdpred_age2", "thy_fdpred_age3",
              "y3_mean_pred_age1", "y3_mean_pred_age2", "y3_mean_pred_age3", "per_countspred_age1", "per_countspred_age2", "per_countspred_age3",
              "y4_mean_pred_age1",  "y4_mean_pred_age2",  "y4_mean_pred_age3",  "per_fdpred_age1", "per_fdpred_age3", "per_fdpred_age3",
              "log_lik", "lambdaThy_inv" , "lambdaPer_inv")

parameters <- c(parametersToPlot, otherRVs)

################################################################################################
## run Stan

# parameters for running fits
nChains <- 2
nPost <- 500 ## Number of post-burn-in samples per chain after thinning
nBurn <- 500 ## Number of burn-in samples per chain after thinning
nThin <- 1

nIter <- (nPost + nBurn) * nThin
nBurnin <- nBurn * nThin

fit <- stan(file = file.path(modelDir, paste(modelName, ".stan", sep = "")),
            data = data,
            pars = parameters,
            iter = nIter,
            warmup = nBurnin,
            thin = nThin, 
            init = init,
            chains = nChains,
            control = list(adapt_delta = 0.9))

################################################################################################
# save results in output directory
if (!file.exists(saveDir)){
  dir.create(saveDir)
}

# saving output file as 
output_filename=paste(Sys.getenv("SLURM_JOB_ID"), "_P", Sys.getenv("SLURM_PROCID"), "_",  modelName, ".rds",  sep="")

# saving the stan fit object for individual runs as 
saveRDS(fit, file = file.path(saveDir, output_filename))

## calculating an output from individual runs for the validation of a sucessful run
loo_loglik <- extract_log_lik(fit, parameter_name = "log_lik", merge_chains = TRUE)
loo_ic <- loo(loo_loglik)

## this value is added to the job specific run_info file which is gets written in the save directory after each job copmpletion

# information about the cluster run
cluster_run_info <- paste("Time = ", timestamp(),
                          "|| Job ID = ", Sys.getenv("SLURM_JOB_ID"),
                          "|| Modelname = ", modelName,
                          "|| crude_loo_ic = ", loo_ic$estimates[3],
                          "|| Host = ", Sys.getenv("SLURM_SUBMIT_HOST"),
                          "|| Nodes = ",  Sys.info()["nodename"],
                          "|| Process_ID = P", Sys.getenv("SLURM_PROCID"))

# writing the cluster run info as a data file in the out dir 
write.table(cluster_run_info, file = file.path(saveDir, paste("job_", Sys.getenv("SLURM_JOB_ID"), "_run_info.txt")), append = TRUE)
