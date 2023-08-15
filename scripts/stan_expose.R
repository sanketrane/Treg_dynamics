rm(list = ls()); gc()

library(rstan)
library(tidyverse)

## model
modelName <- "ASM_deltavar_FP3"
outputDir <- file.path("output_fit")
rstan::expose_stan_functions(file.path("stan_models", paste0(modelName, ".stan")))
#params <- c(1.363198e+05, 0.01, 0.04, 0.02, 0.002)
params <- c(1096284, 0.00624, 0.03, 0.017, 0.0012)
par_inc <- c(0.3, 0.01, 0.05, 0.02, 10, 11, 9, 10)
theta <- c(0)
x_i <- c(60)
x_r <- c(100)
math_reduce(params, theta, x_r, x_i)

## ts
ts_pred1 <- 10^seq(log10(66), log10(450), length.out = 300)
ts_pred2 <- 10^seq(log10(91), log10(450), length.out = 300)
ts_pred3 <- 10^seq(log10(90), log10(450), length.out = 300)
ts_pred4 <- 10^seq(log10(174), log10(450), length.out = 300)
tb_pred1 <- rep(45, 300)
tb_pred2 <- rep(66, 300)
tb_pred3 <- rep(76, 300)
tb_pred4 <- rep(118, 300)

kseq <- seq(0, 1, 0.01)
ki_dist_i <- sapply(kseq, ki_dist_theta, 2)
ggplot()+
  geom_point(aes(x=kseq, y=ki_dist_i))


ki_d <- function(ki){
  r_ki_init = 5.69
  exp(-ki * r_ki_init)/((1 - exp(-r_ki_init))/r_ki_init)
}

integrate(ki_d, exp(-1), 1)


chi_vec1 <- sapply(ts_pred1 - tb_pred1, Chi_spline)
chi_vec2 <- sapply(ts_pred2 - tb_pred2, Chi_spline)
chi_vec3 <- sapply(ts_pred3 - tb_pred3, Chi_spline)
chi_vec4 <- sapply(ts_pred4 - tb_pred4, Chi_spline)

ggplot()+
  geom_line(aes(x=ts_pred1 - tb_pred1, y=chi_vec1))

total_counts1 <- N_pooled_time(ts_pred1, tb_pred1, params)
total_counts2 <- N_pooled_time(ts_pred2, tb_pred2, params)
total_counts3 <- N_pooled_time(ts_pred3, tb_pred3, params)
total_counts4 <- N_pooled_time(ts_pred4, tb_pred4, params)

donor_counts1 <- N_donor_time(ts_pred1, tb_pred1, params)
donor_counts2 <- N_donor_time(ts_pred2, tb_pred2, params)
donor_counts3 <- N_donor_time(ts_pred3, tb_pred3, params)
donor_counts4 <- N_donor_time(ts_pred4, tb_pred4, params)

Nfd_pred1 <- donor_counts1/(total_counts1 * chi_vec1)
Nfd_pred2 <- donor_counts2/(total_counts2 * chi_vec2)
Nfd_pred3 <- donor_counts3/(total_counts3 * chi_vec3)
Nfd_pred4 <- donor_counts4/(total_counts4 * chi_vec4)

host_counts_mean1 = N_host_time(ts_pred1, tb_pred1, params);
host_counts_mean2 = N_host_time(ts_pred2, tb_pred2, params);
host_counts_mean3 = N_host_time(ts_pred3, tb_pred3, params);
host_counts_mean4 = N_host_time(ts_pred4, tb_pred4, params);

host_ki_counts1 = U_host_time(ts_pred1, tb_pred1, params);
host_ki_counts2 = U_host_time(ts_pred2, tb_pred2, params);
host_ki_counts3 = U_host_time(ts_pred3, tb_pred3, params);
host_ki_counts4 = U_host_time(ts_pred4, tb_pred4, params);

host_ki_mean1 = host_ki_counts1/host_counts_mean1;
host_ki_mean2 = host_ki_counts2/host_counts_mean2;
host_ki_mean3 = host_ki_counts3/host_counts_mean3;
host_ki_mean4 = host_ki_counts4/host_counts_mean4;

donor_ki_counts1 = U_donor_time(ts_pred1, tb_pred1, params);
donor_ki_counts2 = U_donor_time(ts_pred2, tb_pred2, params);
donor_ki_counts3 = U_donor_time(ts_pred3, tb_pred3, params);
donor_ki_counts4 = U_donor_time(ts_pred4, tb_pred4, params);

donor_ki_mean1 = donor_ki_counts1/donor_counts1;
donor_ki_mean2 = donor_ki_counts2/donor_counts2;
donor_ki_mean3 = donor_ki_counts3/donor_counts3;
donor_ki_mean4 = donor_ki_counts4/donor_counts4;


## compile pred
stan_pred <- data.frame(
  "age.at.S1K" = c(ts_pred1, ts_pred2, ts_pred3, ts_pred4),
  "counts" = c(total_counts1, total_counts2, total_counts3, total_counts4),
  "Nfd" = c(Nfd_pred1, Nfd_pred2, Nfd_pred3, Nfd_pred4),
  "host_ki" = c(host_ki_mean1, host_ki_mean2, host_ki_mean3, host_ki_mean4),
  "donor_ki" = c(donor_ki_mean1, donor_ki_mean2, donor_ki_mean3, donor_ki_mean4),
  "ageBMT_bin" = c(rep("agebin1", 300), rep("agebin2", 300), rep("agebin3", 300), rep("agebin4", 300))) %>%
  gather(c(counts, Nfd, host_ki, donor_ki), key="subdata", value = "naive")



## loading required datasets for plotting
counts_file <- file.path("data", "Counts_Treg.csv")
counts_data <- read.csv(counts_file) %>% 
  arrange(age.at.S1K) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 56, 'agebin1',
                             ifelse(age.at.BMT <= 70, 'agebin2',
                                    ifelse(age.at.BMT <= 84, 'agebin3', 'agebin4'))),
         subdata = "counts") %>%
  #gather(c(naive, memory), key='location', value = 'total_counts')
  select(-memory) 

Nfd_file <- file.path("data", "Nfd_Treg.csv")
Nfd_data <- read.csv(Nfd_file) %>% 
  arrange(age.at.S1K)%>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 56, 'agebin1',
                             ifelse(age.at.BMT <= 70, 'agebin2',
                                    ifelse(age.at.BMT <= 84, 'agebin3', 'agebin4'))),
         subdata = "Nfd") %>%
  select(-contains('memory')) %>%
  rename(naive = naive_DP1)
  

hostki_file <- file.path("data", "hostKi67_Treg.csv")
hostki_data <- read.csv(hostki_file) %>% 
  arrange(age.at.S1K) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 56, 'agebin1',
                             ifelse(age.at.BMT <= 70, 'agebin2',
                                    ifelse(age.at.BMT <= 84, 'agebin3', 'agebin4'))),
         subcomp='Host',
         subdata = "host_ki") %>%
  #gather(c(naive, memory), key='location', value = 'prop_ki')
  select(-memory)

donorki_file <- file.path("data", "donorKi67_Treg.csv")
donorki_data <- read.csv(donorki_file) %>% 
  arrange(age.at.S1K) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 56, 'agebin1',
                             ifelse(age.at.BMT <= 70, 'agebin2',
                                    ifelse(age.at.BMT <= 84, 'agebin3', 'agebin4'))),
         subcomp='Donor',
         subdata = "donor_ki")%>%
  #gather(c(naive, memory), key='location', value = 'prop_ki')
  select(-memory)


ki_data <- rbind(donorki_data, hostki_data)


getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

Nfd_data %>%
  group_by(ageBMT_bin) %>%
  summarise(md = getmode(age.at.BMT))

####plots 
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

legn_labels <- c('6-8', '8-10', '10-12', '12-25')

ggplot()+
  geom_point(data = counts_data, aes(x = age.at.S1K , y=naive, col=ageBMT_bin)) +
  geom_line(data= filter(stan_pred, subdata == "counts"),
            aes(x = age.at.S1K, y=naive, col=ageBMT_bin)) +
  labs(title=paste('Total counts of naive Tregs'),  y=NULL, x= "Host age (days)") + 
  scale_color_discrete(name="Host age at \n BMT (Wks)", labels=legn_labels)+
  scale_x_continuous(limits = c(60, 450) , trans="log10", breaks=c(10, 30, 100, 300)) + #scale_y_log10() +
  scale_y_continuous(limits = c(5e4, 2e7), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  #facet_wrap(~ factor(location, levels = c('naive', "memory")))+
  guides(fill = 'none') + myTheme 


ggsave(filename = file.path("output_fit", paste0(modelName, "_P1.pdf")), last_plot(),
       device = "pdf", height = 6, width = 8.5)


ggplot()+
  geom_point(data = Nfd_data, aes(x = age.at.S1K , y=naive, col=ageBMT_bin)) +
  geom_line(data= filter(stan_pred, subdata == "Nfd"),
            aes(x = age.at.S1K, y=naive, col=ageBMT_bin)) +
  labs(x = "Host age (days)", y = NULL, title = "Normalised Chimerism in naive Tregs") +
  scale_color_discrete(name="Host age at \n BMT (Wks)", labels=legn_labels)+
  scale_x_continuous(limits = c(1, 450), breaks = c(0,100,200,300, 400, 500))+
  scale_y_continuous(limits =c(0, 1.02), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) + 
  # facet_wrap(~ factor(location, levels = c('naive', "memory")))+
  guides(fill='none')+ myTheme


ggsave(filename = file.path("output_fit", paste0(modelName, "_P2.pdf")), last_plot(),
       device = "pdf", height = 6, width = 8.5)
  
ki_pred <- stan_pred %>%
  filter(subdata != "counts")%>%
  filter(subdata != "Nfd") %>%
  mutate(subcomp = ifelse(subdata == "host_ki", "Host", "Donor"))

fac_labels <- c(`agebin1`= '6-8 weeks', `agebin2`= '8-10 weeks', `agebin3`= '10-12 weeks', `agebin4`= '12-25 weeks')

ggplot()+
  geom_point(data = ki_data, aes(x = age.at.S1K , y=naive*100, col=subcomp)) +
  geom_line(data= ki_pred,
            aes(x = age.at.S1K, y=naive*100, col=subcomp))+
  labs(x = "Host age (days)", y = NULL, title = "% Ki67hi in naive Tregs") +
  scale_x_continuous(limits = c(60, 450), breaks = c(0,100,200,300, 400, 500))+
  scale_y_continuous(limits =c(0, 50), breaks = c(0, 10, 20, 30, 40, 50))+ 
  facet_wrap(~ ageBMT_bin, scales = 'free', labeller = as_labeller(fac_labels))+
  guides(fill='none') + myTheme + theme(legend.title = element_blank())


ggsave(filename = file.path("output_fit", paste0(modelName, "_P3.pdf")), last_plot(),
       device = "pdf", height = 6, width = 8.5)

ageseq <- seq(0, 40, 0.01)
gvec <- sapply(ageseq, g_age, parms = params)
ggplot()+
  geom_point(aes(x=ageseq, y=gvec))+
  scale_y_log10(limits=c(1e4, 1e5))

kiseq <- seq(0, 1, 0.01)
kvec <- sapply(kiseq, ki_dist_init)
ggplot()+
  geom_point(aes(x=kiseq, y=kvec)) 


ktheta_tot <- sapply(kiseq, ki_dist_theta, subpop=1)
ktheta_donor <- sapply(kiseq, ki_dist_theta, subpop=2)
ktheta_host <- sapply(kiseq, ki_dist_theta, subpop=3)
ggplot() +
  geom_point(aes(x=kiseq, y=ktheta_tot), col=2) +
  geom_point(aes(x=kiseq, y=ktheta_donor), col=3) +
  geom_point(aes(x=kiseq, y=ktheta_host), col=4) 





