## clearing the environment
rm(list = ls())  
gc()    

library(rstan)
library(tidyverse)
####################################################################################

#### plotting style
myTheme <- theme(text = element_text(size = 12), axis.text = element_text(size = 12),
                 axis.title =  element_text(size = 12, face = "bold"),
                 plot.title = element_text(size=12,  hjust = 0.5, face = "bold"),
                 legend.background = element_blank(), legend.key = element_blank())

# setting ggplot theme for rest fo the plots
theme_set(theme_bw())


####### plotting
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

log10minorbreaks=as.numeric(1:10 %o% 10^(3:8))

## Setting all the directories for opeartions
projectDir <- getwd()
dataDir <- file.path(projectDir, "datafiles")
outputDir <- file.path(projectDir, "output_fit")
saveDir <- file.path(projectDir, 'save_csv')


## model specific details that needs to be change for every run
M1 <- "Incumbent_logit"
M2 <- "Incumbent"
M3 <- "Incumbent_simple"

# compiling multiple stan objects together that ran on different nodes
stanfitM1.1 <- read_stan_csv(file.path(saveDir, paste0(M1, "_c1", ".csv")))
stanfitM1.2 <- read_stan_csv(file.path(saveDir, paste0(M1, "_c2",".csv")))
stanfitM1.3 <- read_stan_csv(file.path(saveDir, paste0(M1, "_c3",".csv")))

fit1 <- sflist2stanfit(list(stanfitM1.3))

# finding the parameters used in the model 
# using the last parameter("sigma4") in the array to get the total number of parameters set in the model
num_pars_M1 <- which(fit1@model_pars %in% "sigma_counts_per") -1      
params_M1 <- fit1@model_pars[1:num_pars_M1]


# compiling multiple stan objects together that ran on different nodes
stanfitM2.1 <- read_stan_csv(file.path(saveDir, paste0(M2, "_c1", ".csv")))
stanfitM2.2 <- read_stan_csv(file.path(saveDir, paste0(M2, "_c2",".csv")))
stanfitM2.3 <- read_stan_csv(file.path(saveDir, paste0(M2, "_c3",".csv")))

fit2 <- sflist2stanfit(list(stanfitM2.1, stanfitM2.2, stanfitM2.3))

# finding the parameters used in the model 
# using the last parameter("sigma4") in the array to get the total number of parameters set in the model
num_pars_M2 <- which(fit2@model_pars %in% "sigma_counts_per") -1      
params_M2 <- fit2@model_pars[1:num_pars_M2]


# compiling multiple stan objects together that ran on different nodes
stanfitM3.1 <- read_stan_csv(file.path(saveDir, paste0(M3, "_c1", ".csv")))
stanfitM3.2 <- read_stan_csv(file.path(saveDir, paste0(M3, "_c2",".csv")))
stanfitM3.3 <- read_stan_csv(file.path(saveDir, paste0(M3, "_c3",".csv")))

fit3 <- sflist2stanfit(list(stanfitM3.3, stanfitM3.1))

# finding the parameters used in the model 
# using the last parameter("sigma4") in the array to get the total number of parameters set in the model
num_pars_M3 <- which(fit3@model_pars %in% "sigma_counts_per") -1      
params_M3 <- fit3@model_pars[1:num_pars_M3]


ptable1 <- monitor(as.array(fit1, pars = params_M1), warmup = 0, print = FALSE)
out_table1 <- data.frame(ptable1[1:num_pars_M1, c(1, 4, 8)])
names(out_table1) <- c('Estimates', 'par_lb', 'par_ub')

ptable2 <- monitor(as.array(fit2, pars = params_M2), warmup = 0, print = FALSE)
out_table2 <- data.frame(ptable2[1:num_pars_M2, c(1, 4, 8)])
names(out_table2) <- c('Estimates', 'par_lb', 'par_ub')


ptable3 <- monitor(as.array(fit3, pars = params_M3), warmup = 0, print = FALSE)
out_table3 <- data.frame(ptable3[1:num_pars_M3, c(1, 4, 8)])
names(out_table3) <- c('Estimates', 'par_lb', 'par_ub')


df_pars <- rbind(out_table1, out_table2) %>%
  mutate(parname = c(row.names(out_table1), row.names(out_table2)),
         Model = c(rep('M1', num_pars_M1), rep('M2', num_pars_M2))) %>%
  filter(!grepl("_0", parname))

blank_data <- data.frame(parname = rep(c("alpha", "delta", "sigma"), 3),
                         Param = rep(c("Rate of influx", "Loss rate", "Sigma"), 3),
                         Model = c(rep('M1', 3), rep('M2', 3), rep('GT', 3)),
                         Estimates = c(0.008,0.008,0.2,0.05,0.1,0.5, 0, 0, 0))

ggplot(df_pars, aes(y=Estimates, x=factor(Model), col=Model))+
  labs(y=NULL) +
  geom_errorbar(aes(y=Estimates, ymin=par_lb, ymax=par_ub, x=Model),
                width=0.2, linetype=1,  position=position_dodge(0.4)) +
  #geom_blank(data = blank_data)+
  geom_point(position=position_dodge(width=0.4), stat = "identity", size=4) + 
  facet_wrap(~ factor(parname), scales = "free") + 
  expand_limits(y = 0)  +
  myTheme + theme(axis.text.x=element_blank(),
                  axis.title.x=element_blank())




# time sequence for predictions specific to age bins within the data
ts_pred1 <- 10^seq(log10(66), log10(450), length.out = 300)
ts_pred2 <- 10^seq(log10(91), log10(450), length.out = 300)
ts_pred3 <- 10^seq(log10(90), log10(450), length.out = 300)
ts_pred4 <- 10^seq(log10(174), log10(450), length.out = 300)
tb_pred1 <- rep(45, 300)
tb_pred2 <- rep(66, 300)
tb_pred3 <- rep(76, 300)
tb_pred4 <- rep(118, 300)

# naive Treg counts in the thymus with 90% envelopes
Counts_thy_pred1 <- as.data.frame(fit1, pars = c("counts_thy_mean_pred1", "counts_thy_mean_pred2",
                                               "counts_thy_mean_pred3", "counts_thy_mean_pred4")) %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols("timeseries" = c(ts_pred1, ts_pred2, ts_pred3, ts_pred4))%>%
  mutate(ageBMT_bin = ifelse(grepl("pred1", key),"agebin1",
                             ifelse(grepl("pred2", key), "agebin2",
                                    ifelse(grepl("pred3", key), "agebin3", "agebin4"))),
         location = "Thymus",
         Model = "M1") 



# naive Treg counts in the periphery  with 90% envelopes
Counts_per_pred1 <- as.data.frame(fit1, pars = c("counts_per_mean_pred1", "counts_per_mean_pred2",
                                               "counts_per_mean_pred3", "counts_per_mean_pred4")) %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95))  %>%
  bind_cols("timeseries" = c(ts_pred1, ts_pred2, ts_pred3, ts_pred4))%>%
  mutate(ageBMT_bin = ifelse(grepl("pred1", key),"agebin1",
                             ifelse(grepl("pred2", key), "agebin2",
                                    ifelse(grepl("pred3", key), "agebin3", "agebin4"))),
         location = "Periphery",
         Model = "M1")  


# naive Treg counts in the thymus with 90% envelopes
Counts_thy_pred2 <- as.data.frame(fit2, pars = c("counts_thy_mean_pred1", "counts_thy_mean_pred2",
                                                 "counts_thy_mean_pred3", "counts_thy_mean_pred4")) %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols("timeseries" = c(ts_pred1, ts_pred2, ts_pred3, ts_pred4))%>%
  mutate(ageBMT_bin = ifelse(grepl("pred1", key),"agebin1",
                             ifelse(grepl("pred2", key), "agebin2",
                                    ifelse(grepl("pred3", key), "agebin3", "agebin4"))),
         location = "Thymus",
         Model = "M2") 


# naive Treg counts in the periphery  with 90% envelopes
Counts_per_pred2 <- as.data.frame(fit2, pars = c("counts_per_mean_pred1", "counts_per_mean_pred2",
                                                 "counts_per_mean_pred3", "counts_per_mean_pred4")) %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95))  %>%
  bind_cols("timeseries" = c(ts_pred1, ts_pred2, ts_pred3, ts_pred4))%>%
  mutate(ageBMT_bin = ifelse(grepl("pred1", key),"agebin1",
                             ifelse(grepl("pred2", key), "agebin2",
                                    ifelse(grepl("pred3", key), "agebin3", "agebin4"))),
         location = "Periphery",
         Model = "M2")  



Counts_pred <- rbind(Counts_thy_pred1, Counts_thy_pred2, Counts_per_pred1, Counts_per_pred2)

legn_labels <- c('6-8', '8-10', '10-12', '12-25')

ggplot() +
  geom_ribbon(data = Counts_pred, aes(x = timeseries, ymin = lb, ymax = ub, fill = ageBMT_bin), alpha = 0.15)+
  geom_line(data = Counts_pred, aes(x = timeseries, y = median, color = ageBMT_bin)) +
  geom_point(data = counts_data, aes(x = age.at.S1K, y = total_counts, color = ageBMT_bin), size=2) +
  labs(title=paste('Total counts of naive Tregs'),  y=NULL, x= "Host age (days)") + 
  scale_color_discrete(name="Host age at \n BMT (Wks)", labels=legn_labels)+
  scale_x_continuous(limits = c(60, 450) , trans="log10", breaks=c(10, 30, 100, 300))+
  scale_y_continuous(limits = c(5e3, 5e6), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  facet_grid(factor(Model)~factor(location, levels =c('Thymus', "Periphery")))+
  guides(fill = 'none') + myTheme 










