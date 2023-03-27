### extracting predcitions from stanfit object

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
Counts_thy_pred <- as.data.frame(fit, pars = c("counts_thy_mean_pred1", "counts_thy_mean_pred2",
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
         location = "Thymus") 


#Counts_pred <- rbind(Counts_thy_pred)


# Nfd in thymic naive Tregs with 90% envelopes
Nfd_thy_pred <- as.data.frame(fit, pars = c("Nfd_thy_mean_pred1", "Nfd_thy_mean_pred2",
                                            "Nfd_thy_mean_pred3", "Nfd_thy_mean_pred4")) %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95))  %>%
  bind_cols("timeseries" = c(ts_pred1, ts_pred2, 
                             ts_pred3, ts_pred4))%>%
  mutate(ageBMT_bin = ifelse(grepl("pred1", key),"agebin1",
                             ifelse(grepl("pred2", key), "agebin2",
                                    ifelse(grepl("pred3", key), "agebin3", "agebin4"))),
         location = "Thymus")   

#Nfd_pred <- rbind(Nfd_thy_pred)




# naive Treg counts in the periphery  with 90% envelopes
Counts_per_pred <- as.data.frame(fit, pars = c("counts_per_mean_pred1", "counts_per_mean_pred2",
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
         location = "Periphery")  

# Nfd in peripheral naive Tregs with 90% envelopes
Nfd_per_pred <- as.data.frame(fit, pars = c("Nfd_per_mean_pred1", "Nfd_per_mean_pred2",
                                             "Nfd_per_mean_pred3", "Nfd_per_mean_pred4")) %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols("timeseries" = c(ts_pred1, ts_pred2, 
                             ts_pred3, ts_pred4))%>%
  mutate(ageBMT_bin = ifelse(grepl("pred1", key),"agebin1",
                             ifelse(grepl("pred2", key), "agebin2",
                                    ifelse(grepl("pred3", key), "agebin3", "agebin4"))),
         location = "Periphery")  




Counts_pred <- rbind(Counts_thy_pred, Counts_per_pred)
Nfd_pred <- rbind(Nfd_thy_pred, Nfd_per_pred)


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





# Proportion of Ki67hi in donor thymic naive Tregs with 90% envelopes
ki_donor_thy_pred <- as.data.frame(fit, pars = c("ki_donor_thy_mean_pred1", "ki_donor_thy_mean_pred2",
                                                 "ki_donor_thy_mean_pred3", "ki_donor_thy_mean_pred4")) %>%
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
         subcomp = 'Donor')   

# Proportion of Ki67hi in host thymic naive Tregs with 90% envelopes
ki_host_thy_pred <- as.data.frame(fit, pars = c("ki_host_thy_mean_pred1", "ki_host_thy_mean_pred2",
                                                "ki_host_thy_mean_pred3", "ki_host_thy_mean_pred4")) %>%
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
         subcomp = 'Host')    


ki_thy_pred <- rbind(ki_donor_thy_pred, ki_host_thy_pred)


# Proportion of Ki67hi in donor peripheral naive Tregs with 90% envelopes
ki_donor_per_pred <- as.data.frame(fit, pars = c("ki_donor_per_mean_pred1", "ki_donor_per_mean_pred2",
                                                 "ki_donor_per_mean_pred3", "ki_donor_per_mean_pred4")) %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols("timeseries" = c(ts_pred1, ts_pred2, ts_pred3, ts_pred4))%>%
  mutate(ageBMT_bin = ifelse(grepl("pred1", key),"agebin1",
                             ifelse(grepl("pred2", key), "agebin2",
                                    ifelse(grepl("pred3", key), "agebin3", "agebin4"))),
         location = "Periphery", 
         subcomp = 'Donor')   


# Proportion of Ki67hi in host peripheral naive Tregs with 90% envelopes
ki_host_per_pred <- as.data.frame(fit, pars = c("ki_host_per_mean_pred1", "ki_host_per_mean_pred2",
                                                "ki_host_per_mean_pred3", "ki_host_per_mean_pred4")) %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols("timeseries" = c(ts_pred1, ts_pred2, ts_pred3, ts_pred4))%>%
  mutate(ageBMT_bin = ifelse(grepl("pred1", key),"agebin1",
                             ifelse(grepl("pred2", key), "agebin2",
                                    ifelse(grepl("pred3", key), "agebin3", "agebin4"))),
         location = "Periphery", 
         subcomp = 'Host')    

ki_per_pred <- rbind(ki_donor_per_pred, ki_host_per_pred)
