## Simple linear ODE model -- Homogeneous model 
## clearing the environment
rm(list = ls())  
gc()    
#setwd("/opt/mesh/eigg/sanket/Treg_dynamics")

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

modelName <- "INC_model"
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


arcsinsqrt_func <- function(x){
  return(asin(sqrt(x)))
}
################################################################################################
################################################################################################
#source_sorted <- read_csv(file.path(dataDir, data_derived1))%>% arrange(time.post.BMT)
NTregThy_counts <- read_csv(file.path(dataDir, data_derived1))%>% arrange(age.at.S1K)
NTregThy_Nfd <- read_csv(file.path(dataDir, data_derived2))%>% arrange(age.at.S1K)
NTregPer_counts <- read_csv(file.path(dataDir, data_derived3))%>% arrange(age.at.S1K)
NTregPer_Nfd <- read_csv(file.path(dataDir, data_derived4))%>% arrange(age.at.S1K)

### solving ODEs only for unique timepoints 
### also Stan ODE solver throws an error for repetitive time-points 
# unique time points in data for odes solver
unique_times_df <- NTregPer_counts %>% distinct(age.at.S1K, .keep_all = TRUE)

data_time <- NTregPer_Nfd$age.at.S1K                                          # data and solver time 
solve_time <- unique_times_df$age.at.S1K                                       # unique time points to solve ode
time_index <- purrr::map_dbl(data_time, function(x) which(x == solve_time))    # keeping track of index of time point in relation to solve_time


# delay time for age correction 
# for each host with different age at BMT -- time zero is different which is becomes another variable in the ODEs -- is it delayed differential equation now?
# tb_time -- for solving for initial conditions at each tb --  we use the earlist host age at BMT to calculate the initial conditions at each tb
solve_ageAtBMT <- unique_times_df$age.at.BMT   
tb_time <- solve_ageAtBMT %>% unique() %>% sort()                              # unique age at BMT points to solve ode
tb_index <- purrr::map_dbl(solve_ageAtBMT, function(x) which(x == tb_time))    #keeping track of index of  delay time in relation to tb_time


# time sequence for predictions specific to age bins within the data
# In this analysis 3 age bins were selected with mean ages for each bin as 48, 72 and 120 respectively.
counts_agebmt_sorted <- NTregPer_counts %>% arrange(age.at.BMT)
Mean_agebin1 <- round(mean(counts_agebmt_sorted$age.at.BMT[1:10]))
Mean_agebin2 <- round(mean(counts_agebmt_sorted$age.at.BMT[11:34]))
Mean_agebin3 <- round(mean(counts_agebmt_sorted$age.at.BMT[35:43]))

ts_pred1 <- seq(Mean_agebin1, to = 450, by = 1)
tb_pred1 <- rep(Mean_agebin1, length(ts_pred1))
ts_pred2 <- seq(Mean_agebin2, to = 450, by = 1)
tb_pred2 <- rep(Mean_agebin2, length(ts_pred2))
ts_pred3 <- seq(Mean_agebin3, to = 450, by = 1)
tb_pred3 <- rep(Mean_agebin3, length(ts_pred3))

# tb_time_pred -- for solving for initial conditions at each tb (here the mean age of each age bin)
tb_time_pred1 <- c(min(NTregPer_counts$age.at.BMT), Mean_agebin1)
tb_time_pred2 <- c(min(NTregPer_counts$age.at.BMT), Mean_agebin2)
tb_time_pred3 <- c(min(NTregPer_counts$age.at.BMT), Mean_agebin3)


# setting data input according to the precurosr pop used
## these values dictate how splines for counts and chimerism change with time.
## obtained by fitting splines separately to counts and chimerism in the respective precurosr pop
source_list = data.frame("nu" = 0.003951, "theta0" = 15.29, "chiEst" = 0.7935, "qEst" = 0.002259)

## create data set
data <- list(
  numObs = nrow(NTregPer_counts),              # Total number of observations within data
  solve_time = solve_time,                   # unique observations for ode solver
  num_index = length(solve_time),            # to declare dimensions of the array in Stan
  num_tb = length(tb_time),                 # to declare dimensions of the array in Stan
  time_index = time_index,
  tb_index = tb_index,
  ageAtBMT = solve_ageAtBMT,
  tb_time = tb_time,
  dpBMT = NTregPer_counts$time.post.BMT,
  thy_counts = NTregThy_counts$total_counts,  ## data to fit
  thy_Nfd = NTregThy_Nfd$Nfd,                 ## data to fit
  per_counts = NTregPer_counts$total_counts,  ## data to fit
  per_Nfd = NTregPer_Nfd$Nfd,                 ## data to fit
  numPred1 = length(ts_pred1),                # to declare dimensions of the array in Stan
  numPred2 = length(ts_pred2),                # to declare dimensions of the array in Stan
  numPred3 = length(ts_pred3),                # to declare dimensions of the array in Stan
  ts_pred1 = ts_pred1,
  ts_pred2 = ts_pred2,
  ts_pred3 = ts_pred3,
  tb_pred1 = tb_pred1,
  tb_pred2 = tb_pred2,
  tb_pred3 = tb_pred3,
  tb_time_pred1 = tb_time_pred1,
  tb_time_pred2 = tb_time_pred2,
  tb_time_pred3 =tb_time_pred3,
  thy_Nd0 = 0,                      ## initial condition for donor pool in thymus
  per_Nd0 = 0,                      ## initial condition for donor pool in periphery
  theta0 = exp(source_list$theta0),
  nu = source_list$nu,
  chiEst = source_list$chiEst,
  qEst = source_list$qEst
)

## create initial estimates
init <- function() list(
  psi = exp(rnorm(1, log(0.3), 0.5)),
  delta = exp(rnorm(1,log(0.05), 0.5)),
  lambda = exp(rnorm(1,log(0.03), 0.5)),
  mu = exp(rnorm(1, log(0.3), 0.5)),
  f_inc = exp(rnorm(1, log(0.2), 0.5)),
  
  X0Log = rnorm(1, 13 , 0.1),
  Y0Log = rnorm(1, 10 , 0.1),
  I0Log = rnorm(1, 7 , 0.1),
  
  sigma1 = exp(rnorm(1,log(1.5), 1)),
  sigma2 = exp(rnorm(1,log(1.5), 1)),
  sigma3 = exp(rnorm(1,log(1.5), 1)),
  sigma4 = exp(rnorm(1,log(1.5), 1)))

## Specify the variables for which you want history and density plots
parametersToPlot <- c("psi", "mu", "f_inc", "lambda", "delta", "X0Log", "Y0Log", "I0Log", "sigma1", "sigma2",  "sigma3", "sigma4")
## Additional variables to monitor
## Additional variables to monitor
otherRVs <- c("y1_mean_pred_age1", "y1_mean_pred_age2", "y1_mean_pred_age3", 
              "y2_mean_pred_age1", "y2_mean_pred_age2", "y2_mean_pred_age3",
              "y3_mean_pred_age1",  "y3_mean_pred_age2",  "y3_mean_pred_age3", 
              "y4_mean_pred_age1",  "y4_mean_pred_age2", "y4_mean_pred_age3", "log_lik",
              "log_lik1", "log_lik2", "log_lik3", "log_lik4")

parameters <- c(parametersToPlot, otherRVs)

################################################################################################
## run Sta

# parameters for running fits
nChains <- 6
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
            control = list(adapt_delta = 0.9),
            chains = nChains)

################################################################################################
# save results in output directory
if (!file.exists(saveDir)){
  dir.create(saveDir)
}

# saving output file as 
#output_filename=paste(Sys.getenv("SLURM_JOB_ID"), "_P", Sys.getenv("SLURM_PROCID"), "_",  modelName, ".rds",  sep="")
output_filename=paste0(modelName, ".rds")

# saving the stan fit object for individual runs as 
saveRDS(fit, file = file.path(saveDir, output_filename))

## calculating an output from individual runs for the validation of a sucessful run
loo_loglik1 <- extract_log_lik(fit, parameter_name = "log_lik1", merge_chains = TRUE)
loo_loglik2 <- extract_log_lik(fit, parameter_name = "log_lik1", merge_chains = TRUE)
loo_loglik3 <- extract_log_lik(fit, parameter_name = "log_lik1", merge_chains = TRUE)
loo_loglik4 <- extract_log_lik(fit, parameter_name = "log_lik1", merge_chains = TRUE)
loo_ic <- loo(cbind(loo_loglik1, loo_loglik2, loo_loglik3, loo_loglik4))
print(loo_ic)

### posterior distributions of parameters
ptable <- monitor(as.array(fit, pars = parametersToPlot), warmup = 0, print = FALSE)
print(ptable[1:10, c(1,4,8)])
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
