## Fitting splines to SP4 data to estimate the means

## cleaning environment and cache
rm(list = ls()); gc()

#Loading required libraries
library(tidyverse)

# importing data ---------------------------------
### We consider thymic FoxP3 negative SP4 as the source for Tregs (naive Tregs) 
### for all purposes henceforth -- SP4 = FoxP3 negative SP4

## loading pre-processed data and munging
### different data frames because the number of viable observations per dataset are different
# thynaiveTregs_counts <- read_csv("data/Counts_Treg.csv") 
# thynaiveTregs_chimerism <- read_csv("data/Nfd_Treg.csv")
# thynaiveTregs_donorki <- read_csv("data/donorKi67_Treg.csv") 
# thynaiveTregs_hostki <- read_csv("data/hostKi67_Treg.csv") 


### Ki67 Proportions
source_donorKi67 <- readxl::read_excel(path = "data/master_doc.xlsx", sheet = 7) %>%
  select(contains("mouse.ID"), contains("time"), contains("age"), contains("DP1"), contains("SP4"), contains("Fox25"))%>% 
  na.omit() %>% unique() 

source_hostKi67 <- readxl::read_excel(path = "data/master_doc.xlsx", sheet = 8) %>%
  select(contains("mouse.ID"), contains("time"), contains("age"), contains("DP1"), contains("SP4"), contains("Fox25"))%>% 
  na.omit() %>% unique() 

### Total counts and donor fractions for the source population
source_donor <- readxl::read_excel(path = "data/master_doc.xlsx", sheet = 2) %>%
  select(contains("mouse.ID"), contains("time"), contains("age"), contains("DP1"), contains("SP4"), contains("Fox25"))%>% 
  left_join(source_donorKi67, by = c("mouse.ID", "time.post.BMT", "age.at.S1K", "age.at.BMT"), suffix= c("_counts", "_ki")) %>%
  na.omit() %>% unique() %>%
  mutate(### 1st and 4th quadrants of FOXP3 (x-axis) and CD25 (y-axis) Boolean gate, FoxP3- SP4 = Treg free SP4
    FoxP3_neg_SP4_counts = TH.Fox25Q1_counts + TH.Fox25Q4_counts,
    ### 2nd and 3rd quadrants of FOXP3 (x-axis) and CD25 (y-axis) Boolean gate, FoxP3+ SP4 = SP4 Tregs
    FoxP3_pos_SP4_counts = TH.Fox25Q2_counts + TH.Fox25Q3_counts,
    FoxP3_neg_SP4_ki =  (((TH.Fox25Q1_ki/100) * TH.Fox25Q1_counts) + ((TH.Fox25Q4_ki/100) * TH.Fox25Q4_counts))/(TH.Fox25Q1_counts + TH.Fox25Q4_counts),
    FoxP3_pos_SP4_ki =  (((TH.Fox25Q2_ki/100) * TH.Fox25Q2_counts) + ((TH.Fox25Q3_ki/100) * TH.Fox25Q3_counts))/(TH.Fox25Q2_counts + TH.Fox25Q3_counts))%>%   
  select(-contains("Fox25"))

source_host <- readxl::read_excel(path = "data/master_doc.xlsx", sheet = 3) %>%
  select(contains("mouse.ID"), contains("time"), contains("age"), contains("DP1"), contains("SP4"), contains("Fox25"))%>% 
  left_join(source_hostKi67, by = c("mouse.ID", "time.post.BMT", "age.at.S1K", "age.at.BMT"), suffix= c("_counts", "_ki")) %>%
  na.omit() %>% unique() %>%
  mutate(### 1st and 4th quadrants of FOXP3 (x-axis) and CD25 (y-axis) Boolean gate, FoxP3- SP4 = Treg free SP4
    FoxP3_neg_SP4_counts = TH.Fox25Q1_counts + TH.Fox25Q4_counts,
    ### 2nd and 3rd quadrants of FOXP3 (x-axis) and CD25 (y-axis) Boolean gate, FoxP3+ SP4 = SP4 Tregs
    FoxP3_pos_SP4_counts = TH.Fox25Q2_counts + TH.Fox25Q3_counts,
    FoxP3_neg_SP4_ki =  (((TH.Fox25Q1_ki/100) * TH.Fox25Q1_counts) + ((TH.Fox25Q4_ki/100) * TH.Fox25Q4_counts))/(TH.Fox25Q1_counts + TH.Fox25Q4_counts),
    FoxP3_pos_SP4_ki =  (((TH.Fox25Q2_ki/100) * TH.Fox25Q2_counts) + ((TH.Fox25Q3_ki/100) * TH.Fox25Q3_counts))/(TH.Fox25Q2_counts + TH.Fox25Q3_counts))%>%   
  select(-contains("Fox25"))

# merging total counts for host and donor compartments
# calculating total counts, donor fractions
source_join <- full_join(source_host, source_donor, by = c("mouse.ID", "time.post.BMT", "age.at.S1K", "age.at.BMT"), suffix= c(".host", ".donor")) %>%
  mutate(# total = donor + host
    total_DP1 = TH.DP1_counts.host + TH.DP1_counts.donor,
    total_SP4 = TH.SP4_counts.host + TH.SP4_counts.donor,
    total_FoxP3_neg_SP4 = FoxP3_neg_SP4_counts.host + FoxP3_neg_SP4_counts.donor,
    total_FoxP3_pos_SP4 = FoxP3_pos_SP4_counts.host + FoxP3_pos_SP4_counts.donor,
    ## fd = donor fraction
    fd_DP1 = TH.DP1_counts.donor/total_DP1,
    fd_SP4 = TH.SP4_counts.donor/total_SP4,
    fd_FoxP3_neg_SP4 = FoxP3_neg_SP4_counts.donor/total_FoxP3_neg_SP4,
    fd_FoxP3_pos_SP4 = FoxP3_pos_SP4_counts.donor/total_FoxP3_pos_SP4) %>%
  select(-contains("counts"), -contains("ki"))%>%
  na.omit() %>% unique()


## vectors depicting timeseq for predictions 
timeseq <- seq(40, 450)  ## host age
dpt_seq <- seq(14, 300)   ## days post BMT
  

## fd_fit
## phenomenological function
chi_spline <- function(Time, chiEst, qEst){
  Chi = ifelse(Time-10 < 0, 0,
               chiEst * (1-exp(-qEst * (Time-10))))
  return(Chi)
}

             
## LL function
logit_transf <- function(x){asin(sqrt(x))}
Thy_chi_nlm <- nls((fd_DP1) ~ (chi_spline(time.post.BMT, chiEst, qEst)),
                  data =  source_join,
                  start = list(chiEst=0.84, qEst=0.01))
Thy_chi_pars <- coef(Thy_chi_nlm)

# prediction
Thy_chi_fit <- data.frame(dpt_seq, "y_sp" = chi_spline(dpt_seq, Thy_chi_pars[1], Thy_chi_pars[2]))
Thy_chi_fit_m <- data.frame(dpt_seq, "y_sp" = chi_spline(dpt_seq, 0.8, 0.05))
ggplot() + 
  geom_point(data= source_join, aes(x=age.at.S1K-age.at.BMT, y=fd_DP1), size =2) +
  geom_line(data = Thy_chi_fit, aes(x = dpt_seq , y = y_sp), col=4, size =1) + 
  geom_line(data = Thy_chi_fit_m, aes(x = dpt_seq , y = y_sp), col=2, size =1) + 
  ylim(0,1) + #scale_x_log10(limits=c(1, 350)) +
  labs(title = 'Donor fraction in FoxP3 positive naive SP4 cells',  y=NULL,  x = 'Days post BMT') 

## Counts fit
## phenomenological function
counts_spline <- function(age.at.S1K, b1, b0, nu){
  Time = age.at.S1K   ###
  t0 = 40
  return(10^b0 + (10^b1/(1+ ((Time-t0)/nu)^2)))
}

thy_counts_nlm <- nls(log(total_SP4) ~ log(counts_spline(age.at.S1K, b0, b1, nu)),
                   data =  source_join,
                   start = list(b0=5, b1=6, nu=80))
thy_counts_pars <- coef(thy_counts_nlm)

thycounts_fit <- data.frame(timeseq, "y_sp" = counts_spline(timeseq, thy_counts_pars[1], thy_counts_pars[2], thy_counts_pars[3]))
#thycounts_fit_m <- data.frame(timeseq, "y_sp" = counts_spline(timeseq,4.3, 5.1, 40))

ggplot() + 
  geom_point(data=source_join, aes(x=age.at.S1K, y=total_SP4), col=4, size =2) +
  geom_line(data = thycounts_fit, aes(x = timeseq , y = y_sp), col=4, size =1) + 
  #geom_line(data = thycounts_fit_m, aes(x = timeseq , y = y_sp), col=4, size =1) + 
  #scale_y_log10(limits=c(1e3, 2e5)) + xlim(40, 450)+
  labs(title = 'Counts of thymic SP4 cells',  y=NULL,  x = 'Host age (days)') 


source_donorKi67$TH.SP4 %>% mean()
source_hostKi67$TH.SP4 %>% mean()

## phenomenological function
ki_spline <- function(age.at.S1K, b0, eps_f){
  Time = age.at.S1K   ###
  return(b0 * (1 + exp(-eps_f * (Time-t0))))
}

timeseq <- seq(0, 300)
ki_spline <- function(Time, b0, b1, eps_f){
  return(b0 + (b1/(1 + (Time/eps_f)^4)))
}


ki_nlm <- nls((TH.SP4/100) ~ (ki_spline(age.at.S1K - age.at.BMT, b0, b1, eps_f)),
                      data =  source_donorKi67,
                      start = list(b0=0.4, b1=0.6,  eps_f=60))
ki_pars <- coef(ki_nlm)

ki_fit <- data.frame(timeseq, "y_sp" = ki_spline(timeseq, ki_pars[1], ki_pars[2], ki_pars[3]))
ki_fit_m <- data.frame(timeseq, "y_sp" = ki_spline(timeseq, 0.37, 0.63, 20))

ggplot() + 
  geom_point(data=source_donorKi67, aes(x=age.at.S1K-age.at.BMT, y=TH.SP4), col=4, size =2) +
  #geom_line(data = ki_fit, aes(x = timeseq , y = y_sp*100), col=4, size =1) + 
  geom_line(data = ki_fit_m, aes(x = timeseq , y = y_sp*100), col=4, size =1) + 
  scale_y_log10(limits=c(10, 100)) + xlim(0, 300)+
  labs(title = 'Ki poportions in donor SP4 cells',  y=NULL,  x = 'Host age (days)') 


naiTreg_host_ki <- source_hostKi67$TH.SP4 %>% mean()

ggplot() + 
  #geom_point(data=source_donorKi67, aes(x=age.at.S1K, y=TH.SP4/100), col=4, size =2) +
  geom_point(data=source_hostKi67, aes(x=age.at.S1K, y=TH.SP4/100), col=2, size =2) +
  geom_line(aes(x = dpt_seq , y = naiTreg_host_ki), col=2, size =1) + 
  xlim(40, 450) + ylim(0,1)+
  labs(title = 'Total counts of thymic FoxP3 positive naive SP4 cells',  y=NULL,  x = 'Host age (days)') 


dev.off()

