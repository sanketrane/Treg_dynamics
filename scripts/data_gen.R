### data wrangling for Marginal Zone compartmnet
rm(list = ls()); gc();

library(tidyverse)

## importing data to be fitted 
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


## Unique time points with indices to map
unique_times_counts <- counts_data %>% distinct(age.at.S1K, .keep_all = TRUE) 
data_time_counts <- counts_data$age.at.S1K 
solve_time_counts <- unique_times_counts$age.at.S1K  ## unique time points in the data
unique_times_chi <- Nfd_data %>% distinct(age.at.S1K, .keep_all = TRUE) 
data_time_chi <- Nfd_data$age.at.S1K 
solve_time_chi <- unique_times_chi$age.at.S1K  ## unique time points in the data
unique_times_donorki <- donorki_data %>% distinct(age.at.S1K, .keep_all = TRUE) 
data_time_donorki <- donorki_data$age.at.S1K 
solve_time_donorki <- unique_times_donorki$age.at.S1K  ## unique time points in the data
unique_times_hostki <- hostki_data %>% distinct(age.at.S1K, .keep_all = TRUE) 
data_time_hostki <- hostki_data$age.at.S1K 
solve_time_hostki <- unique_times_hostki$age.at.S1K  ## unique time points in the data

###
solve_time <- solve_time_counts   ### since all other solvetimes are subset of solve_time_counts

## Map of the unique time points on all the timepoints
time_index_counts <- purrr::map_dbl(data_time_counts, function(x) which(x == solve_time))    # keeping track of index of time point in relation to solve_time
time_index_chi <- purrr::map_dbl(data_time_chi, function(x) which(x == solve_time))    # keeping track of index of time point in relation to solve_time
time_index_donorki <- purrr::map_dbl(data_time_donorki, function(x) which(x == solve_time))    # keeping track of index of time point in relation to solve_time
time_index_hostki <- purrr::map_dbl(data_time_hostki, function(x) which(x == solve_time))    # keeping track of index of time point in relation to solve_time

## Data to import in Stan
numObs1 <- length(data_time_counts)
numObs2 <- length(data_time_chi)
numObs3 <- length(data_time_donorki)
numObs4 <- length(data_time_hostki)
n_solve <- length(solve_time)
n_shards <- n_solve
dpBMT <- unique_times_counts$age.at.S1K - unique_times_counts$age.at.BMT
ageAtBMT <- unique_times_counts$age.at.BMT 
counts_per <- counts_data$Periphery
counts_thy <- counts_data$Thymus
Nfd_per <- Nfd_data$Periphery + 1e-4
Nfd_thy <- Nfd_data$Thymus + 1e-4
ki_donor_per <- donorki_data$Periphery
ki_donor_thy <- donorki_data$Thymus
ki_host_per <- hostki_data$Periphery
ki_host_thy <- hostki_data$Thymus

logit_transf <- function(x){log(x/(1-x))}
asinsq_transf <- function(x){asin(sqrt(x))}
<<<<<<< HEAD
logit_transf(Nfd_thy)
=======
logit_transf(Nfd_per)
>>>>>>> 4e8f318e85797612de579e3e2ea67b5219212a34
plot(log(counts_thy)~data_time_counts, ylim = c(1, 16))
points(log(counts_per)~data_time_counts, col=2)

## defining the function to calculate mode of a vector series
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

counts_data %>%
  group_by(ageBMT_bin) %>%
  summarise('modevalue' = getmode(age.at.BMT),
            'ub' = max(age.at.BMT),
            'lb' = min(age.at.BMT),
            'counts' = n())


# time sequence for predictions specific to age bins within the data
ts_pred1 <- 10^seq(log10(66), log10(450), length.out = 300)
ts_pred2 <- 10^seq(log10(91), log10(450), length.out = 300)
ts_pred3 <- 10^seq(log10(90), log10(450), length.out = 300)
ts_pred4 <- 10^seq(log10(174), log10(450), length.out = 300)
tb_pred1 <- rep(45, 300)
tb_pred2 <- rep(66, 300)
tb_pred3 <- rep(76, 300)
tb_pred4 <- rep(118, 300)
numPred <- length(ts_pred1)


rstan::stan_rdump(c("numObs1", "numObs2", "numObs3", "numObs4", "n_solve", "n_shards", "solve_time", "dpBMT", "ageAtBMT", 
             "time_index_counts", "time_index_chi", "time_index_donorki", "time_index_hostki", 
             "counts_per",  "Nfd_per", "ki_donor_per", "ki_host_per",
             "counts_thy",  "Nfd_thy", "ki_donor_thy", "ki_host_thy",
             "ts_pred1", "ts_pred2", "ts_pred3", "ts_pred4",
             "tb_pred1", "tb_pred2", "tb_pred3",  "tb_pred4",  "numPred"),
           file = file.path('data', paste0('Treg_data_S', n_shards,".Rdump")))











