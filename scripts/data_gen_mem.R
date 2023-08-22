## memory Tregs
####### data generation for stan fitting procedure 
rm(list = ls()); gc();
## laoding libraries
library(tidyverse)

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


## importing data to be fitted 
Nfd_file <- file.path("data", "Treg_memory_Nfd.csv")
Treg_memory_Nfd <- read.csv(Nfd_file) %>% 
  arrange(age.at.S1K) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 56, 'agebin1',
                             ifelse(age.at.BMT <= 70, 'agebin2',
                                    ifelse(age.at.BMT <= 84, 'agebin3', 'agebin4'))))

legn_labels <- c('6-8', '8-10', '10-12', '12-25')
ggplot() +
  geom_point(data = Treg_memory_Nfd, aes(x = age.at.S1K, y = total_counts, color = ageBMT_bin), size=2) +
  labs(title=paste('Total counts of memory Tregs'),  y=NULL, x= "Host age (days)") + 
  scale_color_discrete(name="Host age at \n BMT (Wks)", labels=legn_labels)+
  scale_x_continuous(limits = c(60, 450) , trans="log10", breaks=c(10, 30, 100, 300)) + #scale_y_log10() +
  scale_y_continuous(limits = c(5e4, 5e6), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  guides(fill = 'none') + myTheme 

ggplot() +
  geom_point(data = Treg_memory_Nfd, aes(x = age.at.S1K, y = Nfd, color = ageBMT_bin), size=2) +
  labs(x = "Host age (days)", y = NULL, title = "Normalised Chimerism in memory Tregs") +
  scale_color_discrete(name="Host age at \n BMT (Wks)", labels=legn_labels)+
  scale_x_continuous(limits = c(1, 450), breaks = c(0,100,200,300, 400, 500))+
  scale_y_continuous(limits =c(0, 1.02), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) + 
  guides(fill='none')+ myTheme



ki_file <- file.path("data", "Treg_memory_ki.csv")
Treg_memory_ki <- read.csv(ki_file) %>% 
  arrange(age.at.S1K) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 56, 'agebin1',
                             ifelse(age.at.BMT <= 70, 'agebin2',
                                    ifelse(age.at.BMT <= 84, 'agebin3', 'agebin4')))) 


fac_labels <- c(`agebin1`= '6-8 weeks', `agebin2`= '8-10 weeks', `agebin3`= '10-12 weeks', `agebin4`= '12-25 weeks')
Treg_memory_ki %>%
  select(-Popln) %>% rename(Donor = ki_donor, Host = ki_host) %>%
  gather(c(Donor, Host), key = "subcomp", value = "ki_prop") %>%
  ggplot() +
  geom_point(aes(x = age.at.S1K, y = ki_prop * 100, color = subcomp), size=1.5) +
  labs(x = "Host age (days)", y = NULL, title = "% Ki67hi in memory Tregs") +
  scale_x_continuous(limits = c(60, 450), breaks = c(0,100,200,300, 400, 500))+
  scale_y_continuous(limits =c(0, 50), breaks = c(0, 10, 20, 30, 40, 50))+ 
  facet_wrap(~ ageBMT_bin, scales = 'free', labeller = as_labeller(fac_labels))+
  guides(fill='none') + myTheme + theme(legend.title = element_blank())


## Unique time points with indices to map
unique_times_df <- Treg_memory_Nfd %>% distinct(age.at.S1K, .keep_all = TRUE) 
data_time_Nfd <- Treg_memory_Nfd$age.at.S1K 
solve_time <- unique_times_df$age.at.S1K  ## unique time points in the data
unique_times_ki <- Treg_memory_ki %>% distinct(age.at.S1K, .keep_all = TRUE) 
data_time_ki <- Treg_memory_ki$age.at.S1K 

## Map of the unique time points on all the timepoints
time_index_chi <- purrr::map_dbl(data_time_Nfd, function(x) which(x == solve_time))    # keeping track of index of time point in relation to solve_time
time_index_ki <- purrr::map_dbl(data_time_ki, function(x) which(x == solve_time))    # keeping track of index of time point in relation to solve_time

## Data to import in Stan
numObs1 <- length(data_time_Nfd)
numObs2 <- length(data_time_ki)
n_solve <- length(solve_time)
n_shards <- n_solve
dpBMT <- unique_times_df$age.at.S1K - unique_times_df$age.at.BMT
ageAtBMT <- unique_times_df$age.at.BMT 
counts <- Treg_memory_Nfd$total_counts
Nfd <- Treg_memory_Nfd$Nfd + 1e-4
ki_donor <- Treg_memory_ki$ki_donor
ki_host <- Treg_memory_ki$ki_host


## defining the function to calculate mode of a vector series
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

Treg_memory_Nfd %>%
  group_by(ageBMT_bin) %>%
  summarise('modevalue' = getmode(age.at.BMT),
            'ub' = max(age.at.S1K),
            'lb' = min(age.at.S1K),
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

rstan::stan_rdump(c("numObs1", "numObs2", "n_solve", "n_shards", 
                    "solve_time", "dpBMT", "ageAtBMT", "time_index_chi", "time_index_ki", 
                    "counts",  "Nfd", "ki_donor", "ki_host",
                    "ts_pred1", "ts_pred2", "ts_pred3", "ts_pred4",
                    "tb_pred1", "tb_pred2", "tb_pred3",  "tb_pred4",  "numPred"),
                  file = file.path('data', paste0('Treg_memory_shards', n_shards,".Rdump")))

