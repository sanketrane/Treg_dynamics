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

### Total counts and donor fractions for the source population
df_donor <- readxl::read_excel(path = "data/cleaned_doc.xlsx", sheet = 2) %>%
  mutate(
    periph_conv_naive = `SP+LN 4nai`,
    periph_conv_memory = `4mem (em+cm)`,
    Spleen_Treg_naive = SP.naiTreg,
    LN_Treg_naive = LN.naiTreg,
    periph_Treg_naive = SP.naiTreg + LN.naiTreg,
    Spleen_Treg_memory = SP.memTreg,
    LN_Treg_memory = LN.memTreg,
    periph_Treg_memory = SP.memTreg + LN.memTreg) %>%   
  select(mouse.ID, contains("age"), TH.DP1, contains("naive"), contains("memory"))

df_host <- readxl::read_excel(path = "data/cleaned_doc.xlsx", sheet = 3) %>%
  mutate(
    periph_conv_naive = `SP+LN 4nai`,
    periph_conv_memory = `4mem (em+cm)`,
    Spleen_Treg_naive = SP.naiTreg,
    LN_Treg_naive = LN.naiTreg,
    periph_Treg_naive = SP.naiTreg + LN.naiTreg,
    Spleen_Treg_memory = SP.memTreg,
    LN_Treg_memory = LN.memTreg,
    periph_Treg_memory = SP.memTreg + LN.memTreg) %>%   
  select(mouse.ID, contains("age"), TH.DP1, contains("naive"), contains("memory"))


# merging total counts for host and donor compartments
# calculating total counts, donor fractions
df_join <- full_join(df_host, df_donor, by = c("mouse.ID", "age.at.S1K", "age.at.BMT"), suffix= c(".host", ".donor")) %>%
  mutate(# total = donor + host
    Thymic_DP1 = TH.DP1.host + TH.DP1.donor,
    counts_conv_naive = periph_conv_naive.host + periph_conv_naive.donor,
    counts_conv_memory = periph_conv_memory.host + periph_conv_memory.donor,
    counts_Treg_naive = periph_Treg_naive.host + periph_Treg_naive.donor,
    counts_Treg_memory = periph_Treg_memory.host + periph_Treg_memory.donor,
    Spleen_Treg_naive = Spleen_Treg_naive.host + Spleen_Treg_naive.donor,
    LN_Treg_naive = LN_Treg_naive.host + LN_Treg_naive.donor,
    Spleen_Treg_memory = Spleen_Treg_memory.host + Spleen_Treg_memory.donor,
    LN_Treg_memory = LN_Treg_memory.host + LN_Treg_memory.donor,
    ## fd = donor fraction
    fd_DP1 = TH.DP1.donor/Thymic_DP1,
    fd_conv_naive = periph_conv_naive.donor/counts_conv_naive,
    fd_conv_memory = periph_conv_memory.donor/counts_conv_memory,
    fd_Treg_naive = periph_Treg_naive.donor/counts_Treg_naive,
    fd_Treg_memory = periph_Treg_memory.donor/counts_Treg_memory,
    fd_spl_Treg_naive = Spleen_Treg_naive.donor/Spleen_Treg_naive,
    fd_ln_Treg_naive = LN_Treg_naive.donor/LN_Treg_naive,
    fd_spl_Treg_memory = Spleen_Treg_memory.donor/Spleen_Treg_memory,
    fd_ln_Treg_memory = LN_Treg_memory.donor/LN_Treg_memory
    ) %>%
  select(mouse.ID, contains("age"), Thymic_DP1, contains("counts"), contains("fd"))

# df_join %>% na.omit() %>%
#   summarise("cov_naive" = cov(fd_spl_Treg_naive, fd_ln_Treg_naive),
#             "cov_memory" = cov(fd_spl_Treg_memory, fd_ln_Treg_memory),
#             "cor_naive" = cor(fd_spl_Treg_naive, fd_ln_Treg_naive),
#             "cor_memory" = cor(fd_spl_Treg_memory, fd_ln_Treg_memory))
# 
# ### correlation plots
# ggplot() +
#   geom_point(data=df_join, aes(x=fd_spl_Treg_naive, y=fd_ln_Treg_naive,  col=age.at.S1K-age.at.BMT), size=2) +
#   geom_line(aes(x=seq(0,0.8, 0.1), y=seq(0,0.8,0.1)))+
#   scale_color_viridis_c(name = NULL) +
#   labs(x="fd Spleen naive Tregs", y = "fd LN naive Tregs")
# 
# ggplot() +
#   geom_point(data=df_join, aes(x=fd_spl_Treg_memory, y=fd_ln_Treg_memory,  col=age.at.S1K-age.at.BMT), size=2) +
#   geom_line(aes(x=seq(0,0.6, 0.1), y=seq(0,0.6,0.1)))+
#   scale_color_viridis_c(name = NULL) +
#   labs(x="fd Spleen memory Tregs", y = "fd LN memory Tregs")
#   
#   
# total counts of differrnt subpopulations combined from spleen and LN
df_counts <- df_join %>%
  select("mouse.ID", "age.at.S1K", "age.at.BMT", contains('counts')) 

# count_labs <- c(`counts_conv_naive` = "Naive CD4", `counts_conv_memory` = "Memory CD4",
#                 `counts_Treg_naive` = "Naive Treg", `counts_Treg_memory` = "Memory Treg")
# 
# df_counts %>%
#   gather(c(counts_conv_naive, counts_conv_memory, counts_Treg_naive, counts_Treg_memory),
#          key = "Subpop", value="Counts") %>%
#   ggplot() +
#   geom_point(aes(x=age.at.S1K-age.at.BMT, y=Counts,  col=Subpop), size=1.5)  +
#   labs(x="Time post BMT (days)", y  = NULL, title= "Total counts") +
#   scale_y_log10() + xlim(0, 400)+
#   scale_color_discrete(name = NULL, labels = count_labs, 
#                        breaks = c("counts_conv_naive", "counts_Treg_naive",
#                                   "counts_conv_memory", "counts_Treg_memory")) +
#   myTheme + theme(legend.key = element_blank())+
#   facet_wrap(.~Subpop, scales= "free_y", labeller = as_labeller(count_labs)) 
#   

# normalising donor fraction in splenic naive tregs by dividing with the donor fractions in the source compartment
df_fd_Norm <- df_join %>%
  select("mouse.ID", "age.at.S1K", "age.at.BMT", contains('fd')) %>%
  mutate(## fd normalized to fd in DP1 
    Nfd_conv_naive = fd_conv_naive/fd_DP1,
    Nfd_conv_memory = fd_conv_memory/fd_DP1,
    Nfd_Treg_naive = fd_Treg_naive/fd_DP1,
    Nfd_Treg_memory = fd_Treg_memory/fd_DP1) %>% 
  filter(Nfd_Treg_memory <=1.2) %>%
  select(contains("mouse.ID"), contains("age"), contains("Nfd")) %>% na.omit() 
# 
# 
# 
# 
# df_fd_Norm %>%
#   gather(c(Nfd_conv_naive, Nfd_conv_memory, Nfd_Treg_naive, Nfd_Treg_memory),
#          key = "Subpop", value="Nfd") %>%
#   ggplot() +
#   geom_point(aes(x=age.at.S1K-age.at.BMT, y=Nfd,  col=Subpop), size=1.5)  +
#   labs(x="Time post BMT (days)", y  = NULL, title= "fd normalized to chimerism in DP1") +
#   #scale_x_log10() +
#   ylim(0,1) + xlim(0,300)+
#   scale_color_discrete(name = NULL, labels = legen_labs, 
#                        breaks = c("Nfd_conv_naive", "Nfd_Treg_naive", "Nfd_conv_memory", "Nfd_Treg_memory")) +
#   myTheme + theme(legend.key = element_blank())+
#   guides(color = guide_legend(override.aes = list(size=1.7)))
#   

######### Ki67 proportions in donor ######### 
donor_counts <- readxl::read_excel(path = "data/cleaned_doc.xlsx", sheet = 2)  %>% 
  select(mouse.ID, contains("age"), SP.4nai, LN.4nai, SP.4mem, LN.4mem,
         SP.naiTreg, SP.memTreg, LN.naiTreg, LN.memTreg)

## Percent ki67 in donor naive Tregs 
df_donorKi67 <- readxl::read_excel(path = "data/cleaned_doc.xlsx", sheet = 4)  %>% 
  left_join(donor_counts, by = c("mouse.ID", "age.at.BMT", "age.at.S1K"), suffix=c(".ki", ".counts"))  %>%
  mutate(Ki67_conv_naive = ((SP.4nai.ki/100) * SP.4nai.counts + (LN.4nai.ki/100) * LN.4nai.counts)/
           (SP.4nai.counts + LN.4nai.counts),
         Ki67_conv_memory = ((SP.4mem.ki/100) * SP.4mem.counts + (LN.4mem.ki/100) * LN.4mem.counts)/
           (SP.4mem.counts + LN.4mem.counts),
         Ki67_Treg_naive = ((SP.naiTreg.ki/100) * SP.naiTreg.counts + (LN.naiTreg.ki/100) * LN.naiTreg.counts)/
           (SP.naiTreg.counts + LN.naiTreg.counts),
         Ki67_Treg_memory = ((SP.memTreg.ki/100) * SP.memTreg.counts + (LN.memTreg.ki/100) * LN.memTreg.counts)/
           (SP.memTreg.counts + LN.memTreg.counts)) %>%
  select(mouse.ID, contains("age"), contains("Ki67")) 

######### Ki67 proportions in host  ######### 
host_counts <- readxl::read_excel(path = "data/cleaned_doc.xlsx", sheet = 3)  %>% 
  select(mouse.ID, contains("age"), SP.4nai, LN.4nai, SP.4mem, LN.4mem,
         SP.naiTreg, SP.memTreg, LN.naiTreg, LN.memTreg)

## Percent ki67 in donor naive Tregs 
df_hostKi67 <- readxl::read_excel(path = "data/cleaned_doc.xlsx", sheet = 5)  %>% 
  left_join(donor_counts, by = c("mouse.ID", "age.at.BMT", "age.at.S1K"), suffix=c(".ki", ".counts"))  %>%
  mutate(Ki67_conv_naive = ((SP.4nai.ki/100) * SP.4nai.counts + (LN.4nai.ki/100) * LN.4nai.counts)/
           (SP.4nai.counts + LN.4nai.counts),
         Ki67_conv_memory = ((SP.4mem.ki/100) * SP.4mem.counts + (LN.4mem.ki/100) * LN.4mem.counts)/
           (SP.4mem.counts + LN.4mem.counts),
         Ki67_Treg_naive = ((SP.naiTreg.ki/100) * SP.naiTreg.counts + (LN.naiTreg.ki/100) * LN.naiTreg.counts)/
           (SP.naiTreg.counts + LN.naiTreg.counts),
         Ki67_Treg_memory = ((SP.memTreg.ki/100) * SP.memTreg.counts + (LN.memTreg.ki/100) * LN.memTreg.counts)/
           (SP.memTreg.counts + LN.memTreg.counts)) %>%
  select(mouse.ID, contains("age"), contains("Ki67")) 


df_ki67 <- df_donorKi67  %>% 
  left_join(df_hostKi67, by = c("mouse.ID", "age.at.BMT", "age.at.S1K"), suffix=c(".donor", ".host"))  
  
fac_labs <- c(`Ki67_conv_naive` = "Naive CD4", `Ki67_conv_memory` = "Memory CD4",
                `Ki67_Treg_naive` = "Naive Treg", `Ki67_Treg_memory` = "Memory Treg")


df_ki67 %>%
  gather(-c(mouse.ID, age.at.S1K, age.at.BMT),key = "popln", value="Ki67") %>%
  mutate(Subcomp = ifelse(grepl("host", popln), "Host", "Donor"),
         Subpop = ifelse(grepl("Ki67_conv_naive", popln), "Ki67_conv_naive",
                         ifelse(grepl("Ki67_conv_memory", popln), "Ki67_conv_memory",
                                ifelse(grepl("Ki67_Treg_naive", popln), "Ki67_Treg_naive",
                                       "Ki67_Treg_memory")))
         ) %>%
  ggplot() +
  geom_point(aes(x=age.at.S1K, y=Ki67,  col=Subcomp), size=1.5)  +
  labs(x="Time post BMT (days)", y  = NULL, title= "Proportions of Ki67+ cells") +
  #scale_x_log10() +
  ylim(0,1) + xlim(0,500) +
  myTheme + theme(legend.key = element_blank())+
  guides(color = guide_legend(override.aes = list(size=1.7)))+
  facet_wrap(.~Subpop, labeller = as_labeller(fac_labs))


### defining source populations
memory_donor <- df_donor %>%
  select(mouse.ID, contains("age"), contains("conv_memory")) %>% na.omit() %>%
  rename(donor_counts = periph_conv_memory)

memory_host <- df_host %>%
  select(mouse.ID, contains("age"), contains("conv_memory")) %>% na.omit() %>%
  rename(host_counts = periph_conv_memory)
  

memory_ki <- df_ki67 %>% 
  select(mouse.ID, contains("age"), contains("conv_memory")) %>% na.omit() %>%
  left_join(memory_donor, by = c("mouse.ID", "age.at.S1K", "age.at.BMT")) %>%
  left_join(memory_host, by = c("mouse.ID", "age.at.S1K", "age.at.BMT")) %>%
  mutate(donor_memory_ki = Ki67_conv_memory.donor * donor_counts,
         host_memory_ki = Ki67_conv_memory.host * host_counts) %>%
  select(mouse.ID, contains("age"), contains("memory_ki"))

naive_donor <- df_donor %>%
  select(mouse.ID, contains("age"), contains("conv_naive")) %>% na.omit() %>%
  rename(donor_counts = periph_conv_naive)

naive_host <- df_host %>%
  select(mouse.ID, contains("age"), contains("conv_naive")) %>% na.omit() %>%
  rename(host_counts = periph_conv_naive)


naive_ki <- df_ki67 %>% 
  select(mouse.ID, contains("age"), contains("conv_naive")) %>% na.omit() %>%
  left_join(naive_donor, by = c("mouse.ID", "age.at.S1K", "age.at.BMT")) %>%
  left_join(naive_host, by = c("mouse.ID", "age.at.S1K", "age.at.BMT")) %>%
  mutate(donor_naive_ki = Ki67_conv_naive.donor * donor_counts,
         host_naive_ki = Ki67_conv_naive.host * host_counts) %>%
  select(mouse.ID, contains("age"), contains("naive_ki"))


Total_Tcells_ki <- memory_ki %>%
  left_join(naive_ki, by = c("mouse.ID", "age.at.S1K", "age.at.BMT"), suffix = c(".n", ".m")) %>%
  mutate(donor_total_ki = donor_memory_ki + donor_naive_ki,
         host_total_ki = host_memory_ki + host_naive_ki) %>%
  select(mouse.ID, contains("age"), contains("total"))
  




####### PLots ########
#######Importing datafiles as .csv for model fitting precess #######
Treg_join %>%
  mutate(naive = round(total_naiveTregs_periph, 0),
         memory = round(total_memoryTregs_periph, 0)) %>% 
  select(-contains("total"), -contains('fd')) %>%
  write.csv(file = "data/Counts_Treg.csv", row.names = FALSE)

Treg_fd_Norm %>%
  mutate(naive = round(naiveTregs_Nfd_periph, 4),
         memory = round(memoryTregs_Nfd_periph, 4)) %>%
  select(-contains("Nfd")) %>%
  write.csv(file = "data/Nfd_Treg.csv", row.names = FALSE)

Treg_hostKi67 %>%
  mutate(naive = round(Ki67_naiveTregs_periph, 4),
         memory = round(Ki67_memoryTregs_periph, 4)) %>%
  select(-contains("Ki67"))  %>%
  write.csv(file = "data/hostKi67_Treg.csv", row.names = FALSE)

Treg_donorKi67 %>%
  mutate(naive = round(Ki67_naiveTregs_periph, 4),
         memory = round(Ki67_memoryTregs_periph, 4)) %>%
  select(-contains("Ki67"))  %>%
  write.csv(file = "data/donorKi67_Treg.csv", row.names = FALSE)
  
source_join %>%
  mutate(FoxP3_Neg_SP4 = round(total_FoxP3_neg_SP4, 0))%>%
  select(-contains("total"), -contains("fd")) %>%
  write.csv(file = "data/Counts_thymicSource.csv", row.names = FALSE)

source_join %>% 
  mutate(FoxP3_Neg_SP4 = round(fd_FoxP3_neg_SP4, 4)) %>%
  select(-contains("total"), contains("fd")) %>%
  write.csv(file = "data/Chimerism_thymicSource.csv", row.names = FALSE)


source_donor %>%
  select(contains("age"), time.post.BMT, contains("neg_SP4_ki")) %>%
  write.csv(file = "data/donorKi67_thymicSource.csv", row.names = FALSE)


source_host %>%
  select(contains("age"), time.post.BMT, contains("neg_SP4_ki")) %>%
  write.csv(file = "data/hostKi67_thymicSource.csv", row.names = FALSE)


