### extracting predcitions from stanfit object

# time sequence for predictions specific to age bins within the data
ts_pred1 <- 10^seq(log10(66), log10(200), length.out = 300)
ts_pred2 <- 10^seq(log10(91), log10(330), length.out = 300)
ts_pred3 <- 10^seq(log10(90), log10(350), length.out = 300)
ts_pred4 <- 10^seq(log10(174), log10(450), length.out = 300)


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
                                    ifelse(grepl("pred3", key), "agebin3", "agebin4")))) 



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
                                    ifelse(grepl("pred3", key), "agebin3", "agebin4")))) 


# Nfd in thymic naive Tregs with 90% envelopes
Nfd_thy_pred <- as.data.frame(fit, pars = c("Nfd_thy_mean_pred1", "Nfd_thy_mean_pred2",
                                             "Nfd_thy_mean_pred3", "Nfd_thy_mean_pred4")) %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95))  %>%
  bind_cols("timeseries" = c(ts_pred1, ts_pred2, ts_pred3, ts_pred4))%>%
  mutate(ageBMT_bin = ifelse(grepl("pred1", key),"agebin1",
                             ifelse(grepl("pred2", key), "agebin2",
                                    ifelse(grepl("pred3", key), "agebin3", "agebin4")))) 

# Nfd in peripheral naive Tregs with 90% envelopes
Nfd_per_pred <- as.data.frame(fit, pars = c("Nfd_per_mean_pred1", "Nfd_per_mean_pred2",
                                             "Nfd_per_mean_pred3", "Nfd_per_mean_pred4")) %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols("timeseries" = c(ts_pred1, ts_pred2, ts_pred3, ts_pred4))%>%
  mutate(ageBMT_bin = ifelse(grepl("pred1", key),"agebin1",
                             ifelse(grepl("pred2", key), "agebin2",
                                    ifelse(grepl("pred3", key), "agebin3", "agebin4")))) 


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
                                    ifelse(grepl("pred3", key), "agebin3", "agebin4")))) 

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
                                    ifelse(grepl("pred3", key), "agebin3", "agebin4")))) 

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
                                    ifelse(grepl("pred3", key), "agebin3", "agebin4")))) 

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
                                    ifelse(grepl("pred3", key), "agebin3", "agebin4")))) 


