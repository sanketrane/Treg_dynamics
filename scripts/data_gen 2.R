library(tidyverse)
library(rstan)

#######
Population <-  'cd4'


## importing data to be fitted 
chimera_file <- file.path("datafiles/original_data", paste0(Population, "_data.csv"))  
chimera_data <- read.csv(chimera_file) %>% 
  arrange(age.at.S1K) %>%
  filter(Nfd <= 1.2)


chi_Nfd_combine <- chimera_data %>%
  select('mouse.ID', "age.at.S1K", "age.at.BMT",  'total_counts', 'Nfd') %>% na.omit() 

#### importing data for GFP positive proportions
chi_data_df <- readxl::read_excel(file.path("datafiles", "pbio.2003949.s003.xlsx"), sheet = 1) %>%
  transform(DP1_donor = as.numeric(DP1_donor),
            DP1_total = as.numeric(DP1_total),
            NCD4_donor = as.numeric(NCD4_donor),
            NCD4_total = as.numeric(NCD4_total))%>%
  mutate(thy_chi = DP1_donor/DP1_total,
         Nfd = NCD4_donor/(NCD4_total * thy_chi)) %>%
  rename(expt_location = expt.location.x.x,
         age.at.S1K = mouse.age.days,
         total_counts = NCD4_total,
         age.at.BMT = age.at.bmt) %>%
  select('mouse.ID', contains('age'), 'total_counts', contains('Nfd'), -contains('group')) %>% na.omit()  %>%
  full_join(chi_Nfd_combine, by=c('mouse.ID', "age.at.S1K", "age.at.BMT"))

 

ggplot()+
  geom_point(data = chimera_data, aes(x=age.at.S1K, y=Nfd), col=2, size=2)+
  geom_point(data = chi_data_df, aes(x=age.at.S1K, y=Nfd_cd4), col=4, size=1)

ontogeny_file <- file.path("datafiles/original_data", paste0(Population, "_ln.csv"))  
ontogeny_data <- read.csv(ontogeny_file) %>%
  rename(age.at.S1K = time, 
         total_counts = counts,
         total_kiprop = ki67) %>%
  mutate(dataset = rep('ontogeny', 34)) %>% 
  arrange(age.at.S1K) 

#unique time points in data for odes solver
unique_times_ont <- ontogeny_data %>% distinct(age.at.S1K, .keep_all = TRUE) 
data_time_ont <- ontogeny_data$age.at.S1K                        # timepoints in observations 
solve_time_ont <- unique_times_ont$age.at.S1K
time_index_ont <- purrr::map_dbl(data_time_ont, function(x) which(x == solve_time_ont))    # keeping track of index of time point in relation to solve_time


unique_times_chi <- chimera_data %>% distinct(age.at.S1K, .keep_all = TRUE) 
data_time_chi <- chimera_data$age.at.S1K 
solve_time_chi <- unique_times_chi$age.at.S1K
time_index_chi <- purrr::map_dbl(data_time_chi, function(x) which(x == solve_time_chi))    # keeping track of index of time point in relation to solve_time


data_fit <- chimera_data %>% 
  select(contains('S1K'), contains('counts'), contains('kiprop')) %>%
  mutate(dataset = rep('chimera', nrow(chimera_data))) %>%
  bind_rows(ontogeny_data) %>% 
  arrange(age.at.S1K) 

## ki67 data
kidata_file <- file.path("datafiles/original_data", paste0(Population, "_donor_host.csv"))  
ki_data <- read.csv(kidata_file) %>%
  right_join(chimera_data) %>% 
  select(-contains("total"), -contains("Nfd"))
donor_ki_df <- ki_data %>% filter(subpop == "donor_ki")
host_ki_df <- ki_data %>% filter(subpop == "host_ki")


## time points for solving the PDE
data_time_ont <- ontogeny_data$age.at.S1K                        # timepoints in observations 
data_time_chi <- chimera_data$age.at.S1K 
tb_chi <- chimera_data$age.at.BMT
ont_counts <- ontogeny_data$total_counts
ont_ki <- ontogeny_data$total_kiprop
chi_counts <- chimera_data$total_counts
N_donor_fraction <- chimera_data$Nfd
donor_ki <- donor_ki_df$ki_prop
host_ki <- host_ki_df$ki_prop
ts_pred_ont <- 10^seq(log10(5), log10(450), length.out = 300)
ts_pred_chi1 <- 10^seq(log10(58), log10(450), length.out = 300)
ts_pred_chi2 <- 10^seq(log10(75), log10(450), length.out = 300)
ts_pred_chi3 <- 10^seq(log10(101), log10(450), length.out = 300)
tb_pred1 <- rep(54, 300)
tb_pred2 <- rep(71, 300)
tb_pred3 <- rep(97, 300)
numPred <- length(ts_pred_chi1)
numOnt <- length(data_time_ont)
numChi <- length(data_time_chi)
dat_t0 <- c(tb_chi)
dat_time <- c(data_time_chi)
numObs <- length(dat_time)

for (n_shards in c(145)) {
  stan_rdump(c("numOnt",  "ont_ki", "ont_counts", 
               "numChi", "chi_counts",  "N_donor_fraction", "donor_ki", "host_ki",
               "ts_pred_ont", "ts_pred_chi1", "ts_pred_chi2", "ts_pred_chi3",
               "tb_pred1", "tb_pred2", "tb_pred3", "numPred", "n_shards", "dat_time", "dat_t0"),
             file = file.path('datafiles', paste0(Population, '_data_s', n_shards,".Rdump")))
} 


#######
Population <-  'cd8'

## importing data to be fitted 
chimera_file <- file.path("datafiles/original_data", paste0(Population, "_data.csv"))  
chimera_data <- read.csv(chimera_file) 

ontogeny_file <- file.path("datafiles/original_data", paste0(Population, "_ln.csv"))  
ontogeny_data <- read.csv(ontogeny_file) %>%
  rename(age.at.S1K = time, 
         total_counts = counts,
         total_kiprop = ki67) %>%
  mutate(dataset = rep('ontogeny', 34))

data_fit <- chimera_data %>% 
  select(contains('S1K'), contains('counts'), contains('kiprop')) %>%
  mutate(dataset = rep('chimera', nrow(chimera_data))) %>%
  bind_rows(ontogeny_data) %>% arrange(age.at.S1K) 

## ki67 data
kidata_file <- file.path("datafiles/original_data", paste0(Population, "_donor_host.csv"))  
ki_data <- read.csv(kidata_file)
donor_ki_df <- ki_data %>% filter(subpop == "donor_ki")
host_ki_df <- ki_data %>% filter(subpop == "host_ki")

## time points for solving the PDE
time_ont <- data_fit$age.at.S1K                        # timepoints in observations 
time_chi <- chimera_data$age.at.S1K 
tb_chi <- chimera_data$age.at.BMT
total_counts <- data_fit$total_counts
total_ki <- data_fit$total_kiprop
N_donor_fraction <- chimera_data$Nfd
donor_ki <- donor_ki_df$ki_prop
host_ki <- host_ki_df$ki_pro
ts_pred_chi1 <- 10^seq(log10(58), log10(450), length.out = 300)
ts_pred_chi2 <- 10^seq(log10(75), log10(450), length.out = 300)
ts_pred_chi3 <- 10^seq(log10(101), log10(450), length.out = 300)
tb_pred1 <- rep(54, 300)
tb_pred2 <- rep(71, 300)
tb_pred3 <- rep(97, 300)
numPred <- length(ts_pred_chi1)
numOnt <- length(time_ont)
numChi <- length(time_chi)

for (n_shards in c(26, 52)) {
  stan_rdump(c("numOnt", "time_ont", "to_counts","total_ki", 
               "numChi", "time_chi", "tb_chi","N_donor_fraction", "donor_ki", "host_ki",
               "ts_pred_ont", "ts_pred_chi1", "ts_pred_chi2", "ts_pred_chi3",
               "numPred", "tb_pred1", "tb_pred2", "tb_pred3",  "n_shards"),
             file = file.path('datafiles', paste0(Population, '_data', numOnt,
                                                  "_s", n_shards,".Rdump")))
} 

## model and data specific details
Population <-  "cd8"
ModelName <- "asm_deltavar"
WorkinDir <- getwd()
DataDir <- file.path(WorkinDir, "datafiles")
ModelDir <- file.path(WorkinDir, "stan_models", "only_chimera", paste0("MAP_", ModelName, "_", Population, ".stan"))
OutputDir <- file.path(WorkinDir, "out_fit", "only_chimera", Population, ModelName)

rstan::expose_stan_functions(ModelDir)
#Importing fitted parameters
ParamsFile <- read.csv(file.path(OutputDir, paste0("params_", Population, "_", ModelName, ".csv")))
params_imported <- ParamsFile$mean[1:4]

theta <- c(0)
x_i <- c(dat_t0[100])
x_r <- c(dat_time[100])
math_reduce(params_imported, theta, x_r, x_i)

dat_time1 <- c(ts_pred_chi1)
dat_time2 <- c(ts_pred_chi2)
dat_time3 <- c(ts_pred_chi3)
dat_t01 <- c(tb_pred1)
dat_t02 <- c(tb_pred2)
dat_t03 <- c(tb_pred3)
numOnt <- length(ts_pred_ont)
numChi <- length(ts_pred_chi1)

logit <- function(x){log(x/(1-x))}

artf_dat1 <- data.frame()
artf_dat1 <- data.frame("x"= math_reduce(params_imported, theta, dat_time1[1], dat_t01[1]))
for (i in 1:length(dat_time1)) {
  artf_dat1[, i] = data.frame("x"=math_reduce(params_imported, theta, dat_time1[i], dat_t01[i]))
}
artf_df1 <- data.frame(t(artf_dat1)) #%>%
  #mutate(resid_counts = log(X1) - log(chimera_data$total_counts),
  #       resid_nfd = logit(X2) - logit(chimera_data$Nfd),
  #       host_age = chimera_data$age.at.S1K
  #)

artf_dat2 <- data.frame()
artf_dat2 <- data.frame("x"= math_reduce(params_imported, theta, dat_time2[1], dat_t02[1]))
for (i in 1:length(dat_time2)) {
  artf_dat2[, i] = data.frame("x"=math_reduce(params_imported, theta, dat_time2[i], dat_t02[i]))
}
artf_df2 <- data.frame(t(artf_dat2))

artf_dat3 <- data.frame()
artf_dat3 <- data.frame("x"= math_reduce(params_imported, theta, dat_time3[1], dat_t03[1]))
for (i in 1:length(dat_time3)) {
  artf_dat3[, i] = data.frame("x"=math_reduce(params_imported, theta, dat_time3[i], dat_t03[i]))
}
artf_df3 <- data.frame(t(artf_dat3))

y3_mean1 <- c(); y3_mean2 <- c(); y3_mean3 <- c(); 
y4_mean1 <- c(); y4_mean2 <- c(); y4_mean3 <- c()
y5_mean1 <- c(); y5_mean2 <- c(); y5_mean3 <- c();
y6_mean1 <- c(); y6_mean2 <- c(); y6_mean3 <- c();

set.seed(1357)

for (i in 1:numChi){
  y3_mean1[i] = artf_df1[i, 1] #+ rnorm(1, 5e5, 5e5)
  y3_mean2[i] = artf_df2[i, 1] #+ rnorm(1, 5e5, 5e5)
  y3_mean3[i] = artf_df3[i, 1] #+ rnorm(1, 5e5, 5e5)
  y4_mean1[i] = artf_df1[ i, 2] #+ rnorm(1, 0.01, 0.01)
  y4_mean2[i] = artf_df2[ i, 2] #+ rnorm(1, 0.01, 0.01)
  y4_mean3[i] = artf_df3[ i, 2] #+ rnorm(1, 0.01, 0.01)
  y5_mean1[i] = artf_df1[i, 3] #+ rnorm(1, 0.02, 0.02)
  y5_mean2[i] = artf_df2[i, 3] #+ rnorm(1, 0.02, 0.02)
  y5_mean3[i] = artf_df3[i, 3] #+ rnorm(1, 0.02, 0.02)
  y6_mean1[i] = artf_df1[ i, 4] #+ rnorm(1, 0.04, 0.04)
  y6_mean2[i] = artf_df2[ i, 4] #+ rnorm(1, 0.04, 0.04)
  y6_mean3[i] = artf_df3[ i, 4] #+ rnorm(1, 0.04, 0.04)
}

## model and data specific details
Population <-  "cd8"
ModelName <- "asm_rhovar"
WorkinDir <- getwd()
DataDir <- file.path(WorkinDir, "datafiles")
ModelDir <- file.path(WorkinDir, "stan_models", "only_chimera", paste0("MAP_", ModelName, "_", Population, ".stan"))
OutputDir <- file.path(WorkinDir, "out_fit", "only_chimera", Population, ModelName)

rstan::expose_stan_functions(ModelDir)
#Importing fitted parameters
ParamsFile <- read.csv(file.path(OutputDir, paste0("params_", Population, "_", ModelName, ".csv")))
params_imported <- ParamsFile$mean[1:4]

artf_dat4 <- data.frame()
artf_dat4 <- data.frame("x"= math_reduce(params_imported, theta, dat_time1[1], dat_t01[1]))
for (i in 1:length(dat_time1)) {
  artf_dat4[, i] = data.frame("x"=math_reduce(params_imported, theta, dat_time1[i], dat_t01[i]))
}
artf_df4 <- data.frame(t(artf_dat4)) 

artf_dat5 <- data.frame()
artf_dat5 <- data.frame("x"= math_reduce(params_imported, theta, dat_time2[1], dat_t02[1]))
for (i in 1:length(dat_time2)) {
  artf_dat5[, i] = data.frame("x"=math_reduce(params_imported, theta, dat_time2[i], dat_t02[i]))
}
artf_df5 <- data.frame(t(artf_dat5))

artf_dat6 <- data.frame()
artf_dat6 <- data.frame("x"= math_reduce(params_imported, theta, dat_time3[1], dat_t03[1]))
for (i in 1:length(dat_time3)) {
  artf_dat6[, i] = data.frame("x"=math_reduce(params_imported, theta, dat_time3[i], dat_t03[i]))
}
artf_df6 <- data.frame(t(artf_dat6))

y3_mean4 <- c(); y3_mean5 <- c(); y3_mean6 <- c(); 
y4_mean4 <- c(); y4_mean5 <- c(); y4_mean6 <- c()
y5_mean4 <- c(); y5_mean5 <- c(); y5_mean6 <- c();
y6_mean4 <- c(); y6_mean5 <- c(); y6_mean6 <- c();


for (i in 1:numChi){
  y3_mean4[i] = artf_df4[i, 1] #+ rnorm(1, 5e5, 5e5)
  y3_mean5[i] = artf_df5[i, 1] #+ rnorm(1, 5e5, 5e5)
  y3_mean6[i] = artf_df6[i, 1] #+ rnorm(1, 5e5, 5e5)
  y4_mean4[i] = artf_df4[ i, 2] #+ rnorm(1, 0.01, 0.01)
  y4_mean5[i] = artf_df5[ i, 2] #+ rnorm(1, 0.01, 0.01)
  y4_mean6[i] = artf_df6[ i, 2] #+ rnorm(1, 0.01, 0.01)
  y5_mean4[i] = artf_df4[i, 3] #+ rnorm(1, 0.02, 0.02)
  y5_mean5[i] = artf_df5[i, 3] #+ rnorm(1, 0.02, 0.02)
  y5_mean6[i] = artf_df6[i, 3] #+ rnorm(1, 0.02, 0.02)
  y6_mean4[i] = artf_df4[ i, 4] #+ rnorm(1, 0.04, 0.04)
  y6_mean5[i] = artf_df5[ i, 4] #+ rnorm(1, 0.04, 0.04)
  y6_mean6[i] = artf_df6[ i, 4] #+ rnorm(1, 0.04, 0.04)
}

counts_df <- data.frame("ageBMT_group1" = y3_mean1,
                         "ageBMT_group2" = y3_mean2,
                         "ageBMT_group3" = y3_mean3) %>%
  gather(key = ageBMT_bin, value = counts) %>%
  bind_cols(age.at.S1K = c(ts_pred_chi1, ts_pred_chi2, ts_pred_chi3))

counts_df2 <- data.frame("ageBMT_group1" = y3_mean4,
                     "ageBMT_group2" = y3_mean5,
                     "ageBMT_group3" = y3_mean6) %>%
  gather(key = ageBMT_bin, value = counts) %>%
  bind_cols(age.at.S1K = c(ts_pred_chi1, ts_pred_chi2, ts_pred_chi3))

Nfd_df <- data.frame("ageBMT_group1" = y4_mean1,
                     "ageBMT_group2" = y4_mean2,
                     "ageBMT_group3" = y4_mean3) %>%
  gather(key = ageBMT_bin, value = Nfd) %>%
  bind_cols(age.at.S1K = c(ts_pred_chi1, ts_pred_chi2, ts_pred_chi3))

Nfd_df2 <- data.frame("ageBMT_group1" = y4_mean4,
                      "ageBMT_group2" = y4_mean5,
                      "ageBMT_group3" = y4_mean6) %>%
  gather(key = ageBMT_bin, value = Nfd) %>%
  bind_cols(age.at.S1K = c(ts_pred_chi1, ts_pred_chi2, ts_pred_chi3))

kidonor_df <- data.frame("ageBMT_group1" = y6_mean1,
                         "ageBMT_group2" = y6_mean2,
                         "ageBMT_group3" = y6_mean3) %>%
  gather(key = ageBMT_bin, value = donor_ki) 

kihost_df <- data.frame("ageBMT_group1" = y5_mean1, 
                        "ageBMT_group2" = y5_mean2, 
                        "ageBMT_group3" = y5_mean3) %>%
  gather(key = ageBMT_bin, value = host_ki) 

kipred_median <- data.frame(kidonor_df, kihost_df) %>%
  select(-"ageBMT_bin.1") %>%
  bind_cols(age.at.S1K = c(ts_pred_chi1, ts_pred_chi2, ts_pred_chi3)) %>%
  gather(-c(age.at.S1K, ageBMT_bin), key = subpop, value = ki_prop) 

kidonor_df2 <- data.frame("ageBMT_group1" = y6_mean4,
                         "ageBMT_group2" = y6_mean5,
                         "ageBMT_group3" = y6_mean6) %>%
  gather(key = ageBMT_bin, value = donor_ki) 

kihost_df2 <- data.frame("ageBMT_group1" = y5_mean4, 
                        "ageBMT_group2" = y5_mean5, 
                        "ageBMT_group3" = y5_mean6) %>%
  gather(key = ageBMT_bin, value = host_ki) 

kipred_median2 <- data.frame(kidonor_df2, kihost_df2) %>%
  select(-"ageBMT_bin.1") %>%
  bind_cols(age.at.S1K = c(ts_pred_chi1, ts_pred_chi2, ts_pred_chi3)) %>%
  gather(-c(age.at.S1K, ageBMT_bin), key = subpop, value = ki_prop) 



##### main plots
ageBMT_names <- c(`ageBMT_group1` = "7-9wks",
                  `ageBMT_group2` = "9-11wks", 
                  `ageBMT_group3` = "11-25wks")
ggplot() +
  geom_point(data = chimera_data, aes(x = age.at.S1K, y = total_counts, col=ageBMT_bin),  alpha=0.8,  size=2) +
  geom_line(data = counts_df, aes(x =age.at.S1K, y = counts, col=ageBMT_bin))+
  geom_line(data = counts_df2, aes(x =age.at.S1K, y = counts, col=ageBMT_bin), linetype=2)+
  scale_color_manual(name=NULL, values=c("#004CA3", "#CD5E99", "#FFA47E")) +
  labs(title= paste0('Counts of naive ', toupper(Population), ' T cells'),  y=NULL, x= "Host age (days)") + 
  facet_grid(. ~ ageBMT_bin, labeller = as_labeller(ageBMT_names)) +
  myTheme + theme(legend.position = c(0.85, 0.15)) + guides(col='none')+
  scale_y_continuous(limits = c(1e6, 1e8), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  myTheme + guides(col="none")


ggsave(filename = file.path(OutputDir, "combined_counts.pdf"), last_plot(), device = "pdf", width = 12, height = 3.5)

ggplot()+
  geom_point(data = chimera_data, aes(x = age.at.S1K, y = Nfd, col = ageBMT_bin),  alpha=0.8,  size =2) + 
  geom_line(data = Nfd_df, aes(x =age.at.S1K, y = Nfd, col=ageBMT_bin))+
  geom_line(data = Nfd_df2, aes(x =age.at.S1K, y = Nfd, col=ageBMT_bin), linetype=2)+
  scale_color_manual(name=NULL, values=c("#004CA3", "#CD5E99", "#FFA47E")) +
  labs(x = 'Host age', y = NULL, title = paste0("Normalised donor fraction in naive ", toupper(Population), " T cells")) +
  scale_y_continuous(limits = c(0, 1.05), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) + 
  #scale_x_log10(limits = c(50, 1450), breaks = c(10, 30, 100, 300, 450)) +
  facet_grid(. ~ ageBMT_bin, labeller = as_labeller(ageBMT_names)) +
  myTheme + theme(legend.position = c(0.85, 0.15)) + guides(col='none')+
  geom_hline(yintercept = 1.0, col='magenta', linetype=2)


ggsave(filename = file.path(OutputDir, "combined_Nfd.pdf"), last_plot(), device = "pdf", width = 12, height = 3.50)

ggplot()+
  geom_line(data = kipred_median, aes(x =age.at.S1K, y = ki_prop*100, col=subpop))+
  geom_line(data = kipred_median2, aes(x =age.at.S1K, y = ki_prop*100, col=subpop), size=1.0, linetype=2)+
  geom_point(data = ki_data, aes(x = age.at.S1K, y = ki_prop*100, col = subpop), alpha=0.8, size =2) + 
  scale_color_manual(name = NULL, values = c("#a81f3c", "#1fa2a8"), labels = c('Donor', "Host")) +
  labs(x = 'Host age', y = NULL, title = paste0("% Ki67high cells in naive ", toupper(Population), " subset")) +
  scale_y_log10(limits = c(0.2, 100), breaks = c(0.1, 1, 10, 100)) +
  facet_grid(. ~ ageBMT_bin, labeller = as_labeller(ageBMT_names)) +
  myTheme + theme(legend.position = c(0.85, 0.15)) + guides(col='none')

ggsave(filename = file.path(OutputDir, "combined_ki.pdf"), last_plot(), device = "pdf", width = 12, height = 3.50)

