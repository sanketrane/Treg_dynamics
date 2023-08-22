### data wrangling for T regulatory cells from Busulfan chimera experiments

## claning environment and cache
rm(list = ls()); gc()

#Loading required libraries
library(tidyverse)

myTheme <-  theme(axis.text = element_text(size = 14),
                  axis.title =  element_text(size = 14, face = "bold"),
                  plot.title = element_text(size = 14, face = "bold",  hjust = 0.5),
                  legend.text = element_text(size = 12),
                  legend.title = element_text(size = 12, face = "bold"),
                  strip.text = element_text(size = 14))


#### Putative precursors for memory Tregs: naive Tregs, conventional naive and total cd4 (naive + memory) 

### Correcting data 

### memory T cells
counts_memory_donor <- readxl::read_excel(path = "data/cleaned_doc.xlsx", sheet = 2) %>%
  select(mouse.ID, contains("age"), TH.DP1, contains("cm"), contains("em"), - contains("+")) 

df_memory_donor <- readxl::read_excel(path = "data/cleaned_doc.xlsx", sheet = 4) %>%
  select(mouse.ID, contains("age"), contains("cm"), contains("em"), - contains("+")) %>%
  left_join(counts_memory_donor, by = c("mouse.ID", "age.at.BMT", "age.at.S1K"), suffix=c(".ki", ".counts")) %>%
  mutate(
    counts_conv_memory = SP.4cm.counts + SP.4em.counts + LN.4cm.counts + LN.4em.counts,
    ki_mem = ((SP.4cm.ki * SP.4cm.counts) + (SP.4em.ki * SP.4em.counts) + (LN.4cm.ki * LN.4cm.counts) + (LN.4em.ki * LN.4em.counts))/counts_conv_memory,
    ki_conv_memory = ifelse(is.na(ki_mem) , SP.4mem.ki, ki_mem),
    counts_Treg_memory = SP.memTreg.counts + LN.memTreg.counts,
    ki_Treg_memory = ((SP.memTreg.ki * SP.memTreg.counts) + (LN.memTreg.ki * LN.memTreg.counts))/counts_Treg_memory) %>%
  select(mouse.ID, contains("age"), TH.DP1, contains("memory")) 

counts_memory_host <- readxl::read_excel(path = "data/cleaned_doc.xlsx", sheet = 3) %>%
  select(mouse.ID, contains("age"), TH.DP1, contains("cm"), contains("em"), - contains("+")) 

df_memory_host <- readxl::read_excel(path = "data/cleaned_doc.xlsx", sheet = 5) %>%
  select(mouse.ID, contains("age"), contains("cm"), contains("em"), - contains("+")) %>%
  left_join(counts_memory_host, by = c("mouse.ID", "age.at.BMT", "age.at.S1K"), suffix=c(".ki", ".counts")) %>%
  mutate(
    counts_conv_memory = SP.4cm.counts + SP.4em.counts + LN.4cm.counts + LN.4em.counts,
    ki_mem = ((SP.4cm.ki * SP.4cm.counts) + (SP.4em.ki * SP.4em.counts) + (LN.4cm.ki * LN.4cm.counts) + (LN.4em.ki * LN.4em.counts))/counts_conv_memory,
    ki_conv_memory = ifelse(is.na(ki_mem) , SP.4mem.ki, ki_mem),
    counts_Treg_memory = SP.memTreg.counts + LN.memTreg.counts,
    ki_Treg_memory = ((SP.memTreg.ki * SP.memTreg.counts) + (LN.memTreg.ki * LN.memTreg.counts))/counts_Treg_memory) %>%
  select(mouse.ID, contains("age"), TH.DP1, contains("memory")) 


### memory T cells
counts_naive_donor <- readxl::read_excel(path = "data/cleaned_doc.xlsx", sheet = 2) %>%
  select(mouse.ID, contains("age"), TH.DP1, contains("nai"), - contains("+")) 

df_naive_donor <- readxl::read_excel(path = "data/cleaned_doc.xlsx", sheet = 4) %>%
  select(mouse.ID, contains("age"),  contains("nai"), - contains("+")) %>%
  left_join(counts_naive_donor, by = c("mouse.ID", "age.at.BMT", "age.at.S1K"), suffix=c(".ki", ".counts")) %>%
  mutate(
    counts_conv_naive = SP.4nai.counts + LN.4nai.counts,
    ki_conv_naive = ((SP.4nai.ki * SP.4nai.counts) + (LN.4nai.ki * LN.4nai.counts))/counts_conv_naive,
    counts_Treg_naive = SP.naiTreg.counts + LN.naiTreg.counts,
    ki_Treg_naive = ((SP.naiTreg.ki * SP.naiTreg.counts) + (LN.naiTreg.ki * LN.naiTreg.counts))/counts_Treg_naive) %>%
  select(mouse.ID, contains("age"), TH.DP1, contains("naive")) 

counts_naive_host <- readxl::read_excel(path = "data/cleaned_doc.xlsx", sheet = 3) %>%
  select(mouse.ID, contains("age"), TH.DP1, contains("nai"), - contains("+")) 

df_naive_host <- readxl::read_excel(path = "data/cleaned_doc.xlsx", sheet = 5) %>%
  select(mouse.ID, contains("age"),  contains("nai"), - contains("+")) %>%
  left_join(counts_naive_host, by = c("mouse.ID", "age.at.BMT", "age.at.S1K"), suffix=c(".ki", ".counts")) %>%
  mutate(
    counts_conv_naive = SP.4nai.counts + LN.4nai.counts,
    ki_conv_naive = ((SP.4nai.ki * SP.4nai.counts) + (LN.4nai.ki * LN.4nai.counts))/counts_conv_naive,
    counts_Treg_naive = SP.naiTreg.counts + LN.naiTreg.counts,
    ki_Treg_naive = ((SP.naiTreg.ki * SP.naiTreg.counts) + (LN.naiTreg.ki * LN.naiTreg.counts))/counts_Treg_naive) %>%
  select(mouse.ID, contains("age"), TH.DP1, contains("naive")) 

#####
source_naive_df <- df_naive_donor %>%
  left_join(df_naive_host, by = c("mouse.ID", "age.at.BMT", "age.at.S1K"), suffix=c(".donor", ".host")) %>%
  select(mouse.ID, contains("age"), contains("conv_naive"), contains("TH.DP1")) %>%
  mutate(
    fd_DP1 = TH.DP1.donor/(TH.DP1.donor + TH.DP1.host),
    total_counts = counts_conv_naive.host + counts_conv_naive.donor,
    donor_ki = (ki_conv_naive.donor/100) * counts_conv_naive.donor,
    host_ki = (ki_conv_naive.host/100) * counts_conv_naive.host,
    donor_counts = counts_conv_naive.donor) %>%
  select(mouse.ID, contains("age"), fd_DP1, total_counts, donor_counts, donor_ki, host_ki)

source_memory_df <- df_memory_donor %>%
  left_join(df_memory_host, by = c("mouse.ID", "age.at.BMT", "age.at.S1K"), suffix=c(".donor", ".host")) %>%
  select(mouse.ID, contains("age"), contains("conv_memory")) %>%
  mutate(
    total_counts = counts_conv_memory.host + counts_conv_memory.donor,
    donor_ki = (ki_conv_memory.donor/100) * counts_conv_memory.donor,
    host_ki = (ki_conv_memory.host/100) * counts_conv_memory.host,
    donor_counts = counts_conv_memory.donor) %>%
  select(mouse.ID, contains("age"),  total_counts, donor_counts, donor_ki, host_ki)


## total T cells
total_Tcell_df <- source_naive_df %>%
  full_join(source_memory_df, by = c("mouse.ID", "age.at.BMT", "age.at.S1K"), suffix=c(".nai", ".mem")) %>%
  mutate(total_counts = total_counts.nai + total_counts.mem,
         donor_counts = donor_counts.nai + donor_counts.mem,
         fd = donor_counts/total_counts,
         Nfd = fd/fd_DP1,
         ki_donor = (donor_ki.nai + donor_ki.mem)/donor_counts,
         ki_host = (host_ki.nai + host_ki.mem)/(total_counts - donor_counts)) %>%
  select(mouse.ID, contains("age"), total_counts, fd, Nfd, ki_donor, ki_host)

total_Tcell_Nfd <- total_Tcell_df %>%
  select(mouse.ID, age.at.BMT, age.at.S1K, total_counts, fd, Nfd) %>% 
  mutate(Popln = "Total_Tcell") %>%
  na.omit() %>% filter(Nfd <=1.25) %>% arrange(age.at.S1K)

total_Tcell_ki <- total_Tcell_df %>% 
  mutate(Popln = "Total_Tcell") %>%
  select(mouse.ID, age.at.BMT, age.at.S1K, ki_donor, ki_host, Popln) %>% na.omit() %>% arrange(age.at.S1K)

## naive T cells
conv_naive_Nfd <- source_naive_df %>%
  mutate(fd = donor_counts/total_counts,
         Nfd = fd/fd_DP1,
         Popln = "naive_conv") %>%
  select(mouse.ID, contains("age"), total_counts, fd, Nfd, Popln) %>% 
  na.omit()  %>% filter(Nfd <=1.25) %>% arrange(age.at.S1K)

conv_naive_ki <- source_naive_df %>%
  mutate(ki_donor = donor_ki/donor_counts,
         ki_host = host_ki/(total_counts - donor_counts),
         Popln = "naive_conv") %>%
  select(mouse.ID, contains("age"), ki_donor, ki_host, Popln) %>% na.omit() %>% arrange(age.at.S1K)


### naive Tregs
naive_Treg_df <- df_naive_donor %>%
  left_join(df_naive_host, by = c("mouse.ID", "age.at.BMT", "age.at.S1K"), suffix=c(".donor", ".host")) %>%
  select(mouse.ID, contains("age"), contains("Treg_naive"), contains("TH.DP1")) %>%
  mutate(
    fd_DP1 = TH.DP1.donor/(TH.DP1.donor + TH.DP1.host),
    total_counts = counts_Treg_naive.host + counts_Treg_naive.donor,
    ki_donor = (ki_Treg_naive.donor/100),
    ki_host = (ki_Treg_naive.host/100),
    fd = counts_Treg_naive.donor/total_counts,
    Nfd = fd/fd_DP1) %>%
  select(mouse.ID, contains("age"), total_counts, ki_donor, ki_host, fd, Nfd)

Treg_naive_Nfd <- naive_Treg_df %>%
  select(mouse.ID, contains("age"), total_counts, fd, Nfd) %>%
  mutate(Popln = "naive_Treg") %>%
  na.omit() %>% filter(Nfd <=1.25) %>% arrange(age.at.S1K)

Treg_naive_ki <- naive_Treg_df %>%
  mutate(Popln = "naive_Treg") %>%
  select(mouse.ID, contains("age"), ki_donor, ki_host, Popln) %>% na.omit()%>% arrange(age.at.S1K)

### memory Tregs
memory_Treg_df <- df_memory_donor %>%
  left_join(df_memory_host, by = c("mouse.ID", "age.at.BMT", "age.at.S1K"), suffix=c(".donor", ".host")) %>%
  select(mouse.ID, contains("age"), contains("Treg_memory"), contains("TH.DP1")) %>%
  mutate(
    fd_DP1 = TH.DP1.donor/(TH.DP1.donor + TH.DP1.host),
    total_counts = counts_Treg_memory.host + counts_Treg_memory.donor,
    ki_donor = (ki_Treg_memory.donor/100),
    ki_host = (ki_Treg_memory.host/100),
    fd = counts_Treg_memory.donor/total_counts,
    Nfd = fd/fd_DP1) %>%
  select(mouse.ID, contains("age"), total_counts, ki_donor, ki_host, fd, Nfd)

Treg_memory_Nfd <- memory_Treg_df %>%
  select(mouse.ID, contains("age"), total_counts, fd, Nfd) %>% 
  mutate(Popln = "memory_Treg") %>%
  na.omit() %>% filter(Nfd <=1.25) %>% arrange(age.at.S1K)

Treg_memory_ki <- memory_Treg_df %>%
  mutate(Popln = "memory_Treg") %>%
  select(mouse.ID, contains("age"), ki_donor, ki_host, Popln) %>% na.omit() %>% arrange(age.at.S1K)


Nfd_all_df <- rbind(total_Tcell_Nfd, conv_naive_Nfd, Treg_naive_Nfd, Treg_memory_Nfd) 

ggplot(Nfd_all_df) +
  geom_point(aes(x=age.at.S1K-age.at.BMT, y=Nfd, col=Popln))+
  ylim(0, 1.05) + labs(x="Days post BMT", y=NULL, title="fd normalized to chimerism in DP1")

ggplot(Nfd_all_df)+
  geom_point(aes(x=age.at.S1K, y=total_counts, col=Popln))+
  labs(x="Host age (days)", y=NULL, title="Cell counts") + scale_y_log10(limits=c(1e5, 8e7)) +
  facet_wrap(.~Popln)


Ki_all_df <- rbind(total_Tcell_ki, conv_naive_ki, Treg_naive_ki, Treg_memory_ki) %>%
  gather(-c(mouse.ID, age.at.S1K, age.at.BMT, Popln), key = "Subpop", value = "Ki_prop")

ggplot(Ki_all_df)+
  geom_point(aes(x=age.at.S1K-age.at.BMT, y=Ki_prop, col=Subpop))+
  ylim(0, 0.5) + labs(x="Days post BMT", y=NULL, title="Proportion of Ki67+ cells ")+
  facet_wrap(.~Popln)


ggplot(Ki_all_df)+
  geom_point(aes(x=age.at.S1K, y=Ki_prop, col=Subpop))+
  ylim(0, 0.5) + labs(x="Host age (days)", y=NULL, title="Proportion of Ki67+ cells ")+
  facet_wrap(.~Popln)

# 
# ### Exporting data files as .csv for model fitting process #######
# Treg_memory_Nfd %>%
#   write.csv(file = "data/Treg_memory_Nfd.csv", row.names = FALSE)
# 
# Treg_memory_ki %>%
#   write.csv(file = "data/Treg_memory_ki.csv", row.names = FALSE)


### splines
## fd_fit
## phenomenological function
chi_spline <- function(Time, chiEst, qEst){
  Chi = ifelse((Time - 10) < 0, 0,
               chiEst * (1 - exp(-qEst * (Time - 10))))
  return(Chi)
}

# naive T cells
dpt_seq <- seq(0, 400)
asq_transf <- function(x){asin(sqrt(x)/sqrt(1.2))}
naive_chi_nlm <- nls(asq_transf(Nfd) ~ asq_transf(chi_spline(age.at.S1K - age.at.BMT, chiEst, qEst)),
                   data =  conv_naive_Nfd,
                   start = list(chiEst=0.84, qEst=0.01))
naive_chi_pars <- coef(naive_chi_nlm)

naive_chi_fit <- data.frame(dpt_seq, "y_sp" = chi_spline(dpt_seq, naive_chi_pars[1], naive_chi_pars[2]))
ggplot() + 
  geom_point(data=conv_naive_Nfd, aes(x=age.at.S1K-age.at.BMT, y=Nfd, col=age.at.BMT), size =2) +
  geom_line(data = naive_chi_fit, aes(x = dpt_seq , y = y_sp), col=4, size =1) + 
  ylim(0,1) + #scale_x_log10(limits=c(10, 750)) +
  labs(title = 'Donor fraction in naive CD4 T cells',  y=NULL,  x = 'Days post BMT') 


# total T cells
Tcell_chi_nlm <- nls(asq_transf(Nfd) ~ asq_transf(chi_spline(age.at.S1K - age.at.BMT, chiEst, qEst)),
                     data =  total_Tcell_Nfd,
                     start = list(chiEst=0.84, qEst=0.01))
Tcell_chi_pars <- coef(Tcell_chi_nlm)
Tcell_chi_fit <- data.frame(dpt_seq, "y_sp" = chi_spline(dpt_seq, Tcell_chi_pars[1], Tcell_chi_pars[2]))
ggplot() + 
  geom_point(data=total_Tcell_Nfd, aes(x=age.at.S1K-age.at.BMT, y=Nfd, col=age.at.BMT), size =2) +
  geom_line(data = Tcell_chi_fit, aes(x = dpt_seq , y = y_sp), col=4, size =1) + 
  ylim(0,1) + #scale_x_log10(limits=c(10, 750)) +
  labs(title = 'Donor fraction in naive CD4 T cells',  y=NULL,  x = 'Days post BMT') 

## naive Tregs
naiTreg_chi_nlm <- nls(asq_transf(Nfd) ~ asq_transf(chi_spline(age.at.S1K - age.at.BMT, chiEst, qEst)),
                     data =  Treg_naive_Nfd,
                     start = list(chiEst=0.84, qEst=0.01))
naiTreg_chi_pars <- coef(naiTreg_chi_nlm)

naiTreg_chi_fit <- data.frame(dpt_seq, "y_sp" = chi_spline(dpt_seq, naiTreg_chi_pars[1], naiTreg_chi_pars[2]))
ggplot() + 
  geom_point(data= Treg_naive_Nfd, aes(x=age.at.S1K-age.at.BMT, y=Nfd, col=age.at.BMT), size =2) +
  geom_line(data = naiTreg_chi_fit, aes(x = dpt_seq , y = y_sp), col=4, size =1) + 
  ylim(0,1) + #scale_x_log10(limits=c(10, 750)) +
  labs(title = 'Donor fraction in naive CD4 T cells',  y=NULL,  x = 'Days post BMT') 


## Counts fit
## phenomenological function
timeseq <- seq(40, 450)
counts_spline <- function(age.at.S1K, basl, nu){
  Time = age.at.S1K   ###
  t0 = 40
  return(10^basl * exp(-nu * (Time - t0)))
}

## naive T cells
naive_counts_nlm <- nls(log(total_counts) ~ log(counts_spline(age.at.S1K, basl, nu)),
                      data =  conv_naive_Nfd,
                      start = list(basl=6, nu=0.03))
naive_counts_pars <- coef(naive_counts_nlm)

naive_counts_fit <- data.frame(timeseq, "y_sp" = counts_spline(timeseq, naive_counts_pars[1], naive_counts_pars[2]))

ggplot() + 
  geom_point(data= conv_naive_Nfd, aes(x=age.at.S1K, y= total_counts), col=4, size =2) +
  geom_line(data = naive_counts_fit, aes(x = timeseq , y = y_sp), col=4, size =1) + 
  scale_y_log10(limits=c(1e5, 8e7)) +
  labs(title = 'Cell counts of naive T cells',  y=NULL,  x = 'Host age (days)') 


### total T cells
Tcell_counts_nlm <- nls(log(total_counts) ~ log(counts_spline(age.at.S1K, basl, nu)),
                        data =  total_Tcell_Nfd,
                        start = list(basl=6, nu=0.03))
Tcell_counts_pars <- coef(Tcell_counts_nlm)

Tcell_counts_fit <- data.frame(timeseq, "y_sp" = counts_spline(timeseq, Tcell_counts_pars[1], Tcell_counts_pars[2]))

ggplot() + 
  geom_point(data= total_Tcell_Nfd, aes(x=age.at.S1K, y= total_counts), col=4, size =2) +
  geom_line(data = Tcell_counts_fit, aes(x = timeseq , y = y_sp), col=4, size =1) + 
  scale_y_log10(limits=c(1e5, 8e7)) +
  labs(title = 'Cell counts total T cells',  y=NULL,  x = 'Host age (days)') 


### naive Tregs
naiTreg_counts_nlm <- nls(log(total_counts) ~ log(counts_spline(age.at.S1K, basl, nu)),
                        data =  Treg_naive_Nfd,
                        start = list(basl=6, nu=0.03))
naiTreg_counts_pars <- coef(naiTreg_counts_nlm)

naiTreg_counts_fit <- data.frame(timeseq, "y_sp" = counts_spline(timeseq, naiTreg_counts_pars[1], naiTreg_counts_pars[2]))

ggplot() + 
  geom_point(data= Treg_naive_Nfd, aes(x=age.at.S1K, y= total_counts), col=4, size =2) +
  geom_line(data = naiTreg_counts_fit, aes(x = timeseq , y = y_sp), col=4, size =1) + 
  scale_y_log10(limits=c(1e5, 8e7)) +
  labs(title = 'Total counts of FoxP3 negative SP4 cells',  y=NULL,  x = 'Host age (days)') 


### Donor Ki67 fit
## phenomenological function
ki_spline <- function(Time, b0, b1, eps_f){
  #Time = age.at.S1K - age.at.BMT
  return(b0 + (b1/(1 + (Time/eps_f)^4)))
}

dpt_seq <- seq(0, 300)

## naive CD4 T cells
naive_ki_vec <- sapply(dpt_seq, ki_spline, b0=0.02, b1=0.98, 14)

ggplot() + 
  geom_point(data= conv_naive_ki, aes(x=age.at.S1K-age.at.BMT, y=ki_donor), col=4, size =2) +
  geom_line(aes(x = dpt_seq , y = naive_ki_vec), col=4, size =1) + 
  scale_y_log10(limits=c(0.003,1)) +
  labs(title = 'Ki67Hi fraction in donor naive CD4 T cells',  y=NULL,  x = 'Days post BMT') 


naive_host_ki <- conv_naive_ki$ki_host %>% mean()

ggplot() + 
  geom_point(data= conv_naive_ki, aes(x=age.at.S1K-age.at.BMT, y=ki_host), col=2, size =2) +
  geom_line(aes(x = dpt_seq , y = naive_host_ki), col=2, size =1) + 
  scale_y_log10(limits=c(0.003,1)) +
  labs(title = 'Ki67Hi fraction in host naive CD4 T cells',  y=NULL,  x = 'Days post BMT') 


## Total CD4
Tcell_ki_vec <- sapply(dpt_seq, ki_spline, b0=0.05, b1=0.95, 14)

ggplot() + 
  geom_point(data= total_Tcell_ki, aes(x=age.at.S1K-age.at.BMT, y=ki_donor), col=4, size =2) +
  geom_line(aes(x = dpt_seq, y = Tcell_ki_vec), col=4, size =1) + 
  scale_y_log10(limits=c(0.003,1)) +
  labs(title = 'Ki67Hi fraction in donor CD4 T cells',  y=NULL,  x = 'Days post BMT') 

Tcell_host_ki <- total_Tcell_ki$ki_host %>% mean()

ggplot() + 
  geom_point(data= total_Tcell_ki, aes(x=age.at.S1K-age.at.BMT, y=ki_host), col=2, size =2) +
  geom_line(aes(x = dpt_seq , y = Tcell_host_ki), col=2, size =1) + 
  scale_y_log10(limits=c(0.003,1)) +
  labs(title = 'Ki67Hi fraction in host CD4 T cells',  y=NULL,  x = 'Days post BMT') 


### naive Tregs
ki_spline2 <- function(Time, b0, eps_f){
  #Time = age.at.S1K - age.at.BMT
  return(b0 * exp(eps_f * Time))
}
naiTreg_ki_vec <- sapply(dpt_seq, ki_spline2, b0=0.085, eps_f = 0.001)

ggplot() + 
  geom_point(data= Treg_naive_ki, aes(x=age.at.S1K-age.at.BMT, y=ki_donor), col=4, size =2) +
  geom_line(aes(x = dpt_seq , y = naiTreg_ki_vec), col=4, size =1) + 
  scale_y_log10(limits=c(0.003,1)) +
  labs(title = 'Ki67Hi fraction in donor naive Treg cells',  y=NULL,  x = 'Days post BMT') 

naiTreg_host_ki <- Treg_naive_ki$ki_host %>% mean()

ggplot() + 
  geom_point(data= Treg_naive_ki, aes(x=age.at.S1K-age.at.BMT, y=ki_host), col=2, size =2) +
  geom_line(aes(x = dpt_seq , y = naiTreg_host_ki), col=2, size =1) + 
  scale_y_log10(limits=c(0.003,1)) +
  labs(title = 'Ki67Hi fraction in host naive Tregs',  y=NULL,  x = 'Days post BMT') 




