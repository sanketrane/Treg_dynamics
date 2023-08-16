### extracting predcitions from stanfit object

# time sequence for predictions specific to age bins within the data
ts_pred1 <- 10^seq(log10(66), log10(450), length.out = 300)
ts_pred2 <- 10^seq(log10(91), log10(450), length.out = 300)
ts_pred3 <- 10^seq(log10(90), log10(450), length.out = 300)
ts_pred4 <- 10^seq(log10(174), log10(450), length.out = 300)
tb_pred1 <- rep(49, 300)
tb_pred2 <- rep(66, 300)
tb_pred3 <- rep(77, 300)
tb_pred4 <- rep(128, 300)

# naive Treg counts in the thymus with 90% envelopes
Counts_naive_pred <- as.data.frame(fit, pars = c("counts_naive_mean_pred1", "counts_naive_mean_pred2",
                                                "counts_naive_mean_pred3", "counts_naive_mean_pred4")) %>%
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

# naive Treg counts in the thymus with 90% envelopes
Counts_naive_withsigma <- as.data.frame(fit, pars = c("counts_naive_pred1", "counts_naive_pred2",
                                               "counts_naive_pred3", "counts_naive_pred4")) %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols("timeseries" = c(ts_pred1, ts_pred2, ts_pred3, ts_pred4)) %>%
  mutate(ageBMT_bin = ifelse(grepl("pred1", key),"agebin1",
                             ifelse(grepl("pred2", key), "agebin2",
                                    ifelse(grepl("pred3", key), "agebin3", "agebin4"))),
         location = "Thymus")   

### selecting rows closest to observed timpoints
find_nearest_time <- function(test, target_vec){
  target.index <- which(abs(target_vec - test) == min(abs(target_vec - test)))
  #target_vec[target.index]
}

obstimes_index <- sapply(counts_data$age.at.S1K, find_nearest_time, target_vec = ts_pred1)

Counts_naive_sigma_obs <- Counts_naive_withsigma[obstimes_index, ]%>%
  mutate(ageBMT_bin = counts_data$ageBMT_bin)
#Counts_pred <- rbind(Counts_naive_pred)


# Nfd in thymic naive Tregs with 90% envelopes
Nfd_naive_pred <- as.data.frame(fit, pars = c("Nfd_naive_mean_pred1", "Nfd_naive_mean_pred2",
                                            "Nfd_naive_mean_pred3", "Nfd_naive_mean_pred4")) %>%
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

#Nfd_pred <- rbind(Nfd_naive_pred)


Counts_pred <- rbind(Counts_naive_pred)
Counts_withsigma <- rbind(Counts_naive_withsigma)
Nfd_pred <- rbind(Nfd_naive_pred)


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
ki_donor_naive_pred <- as.data.frame(fit, pars = c("ki_donor_naive_mean_pred1", "ki_donor_naive_mean_pred2",
                                                 "ki_donor_naive_mean_pred3", "ki_donor_naive_mean_pred4")) %>%
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
ki_host_naive_pred <- as.data.frame(fit, pars = c("ki_host_naive_mean_pred1", "ki_host_naive_mean_pred2",
                                                "ki_host_naive_mean_pred3", "ki_host_naive_mean_pred4")) %>%
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


ki_naive_pred <- rbind(ki_donor_naive_pred, ki_host_naive_pred)

