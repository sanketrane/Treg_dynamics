### data wrangling for T regulatory cells from Busulfan chimera experiments

## claning environment and cache
rm(list = ls()); gc()

#Loading required libraries
library(tidyverse)

## Precursor population dynamics

#### Putative precursors for naive Tregs: Thymic DP1, Thymic FoxP3negative SP4
#### Putative precursors for memory Tregs: Thymic DP1/Thymic FoxP3negative SP4 (representing RTE), naive Tregs

### Ki67 Proportions
source_donorKi67 <- readxl::read_excel(path = "data/master_doc.xlsx", sheet = 7) %>%
  select(contains("mouse.ID"), contains("time"), contains("age"), contains("DP1"), contains("Fox25"))%>% 
  na.omit() %>% unique() 

source_hostKi67 <- readxl::read_excel(path = "data/master_doc.xlsx", sheet = 8) %>%
  select(contains("mouse.ID"), contains("time"), contains("age"), contains("DP1"), contains("Fox25"))%>% 
  na.omit() %>% unique() 

### Total counts and donor fractions for the source population
source_donor <- readxl::read_excel(path = "data/master_doc.xlsx", sheet = 2) %>%
  select(contains("mouse.ID"), contains("time"), contains("age"), contains("DP1"), contains("Fox25"))%>% 
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
  select(contains("mouse.ID"), contains("time"), contains("age"), contains("DP1"), contains("Fox25"))%>% 
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
    total_FoxP3_neg_SP4 = FoxP3_neg_SP4_counts.host + FoxP3_neg_SP4_counts.donor,
    total_FoxP3_pos_SP4 = FoxP3_pos_SP4_counts.host + FoxP3_pos_SP4_counts.donor,
    ## fd = donor fraction
    fd_DP1 = TH.DP1_counts.donor/total_DP1,
    fd_FoxP3_neg_SP4 = FoxP3_neg_SP4_counts.donor/total_FoxP3_neg_SP4,
    fd_FoxP3_pos_SP4 = FoxP3_pos_SP4_counts.donor/total_FoxP3_pos_SP4) %>%
  select(-contains("counts"))%>%
  na.omit() %>% unique()


######### naive Tregs total counts and donor fraction ######### 
## Counts of donor naive Tregs 
Treg_donor <- readxl::read_excel(path = "data/master_doc.xlsx", sheet = 2) %>%
  select(contains("mouse.ID"), contains("time"), contains("age"), contains('Treg')) %>% 
  rename(Spleen_naiveTregs = SP.naiTreg,
         LN_naiveTregs = LN.naiTreg,
         Thymic_naiveTregs = TH.naiTreg,
         Spleen_memoryTregs = SP.memTreg,
         LN_memoryTregs = LN.memTreg,
         Thymic_memoryTregs = TH.memTreg) %>% 
  na.omit() %>% unique()


## Counts of host naive Tregs 
Treg_host <- readxl::read_excel(path = "data/master_doc.xlsx", sheet = 3) %>%
  select(contains("mouse.ID"), contains("time"), contains("age"), contains('Treg')) %>% 
  rename(Spleen_naiveTregs = SP.naiTreg,
         LN_naiveTregs = LN.naiTreg,
         Thymic_naiveTregs = TH.naiTreg,
         Spleen_memoryTregs = SP.memTreg,
         LN_memoryTregs = LN.memTreg,
         Thymic_memoryTregs = TH.memTreg) %>% 
  na.omit() %>% unique()


# merging total counts for host and donor compartments
# calculating total counts, donor fractions 
Treg_join <- full_join(Treg_host, Treg_donor, by = c("mouse.ID", "time.post.BMT", "age.at.S1K", "age.at.BMT"), suffix= c(".host", ".donor"))%>%
  mutate(# total = donor + host
         total_naiveTregs_spl = Spleen_naiveTregs.host + Spleen_naiveTregs.donor,
         total_naiveTregs_ln = LN_naiveTregs.host + LN_naiveTregs.donor,
         total_naiveTregs_periph = total_naiveTregs_spl + total_naiveTregs_ln,  ## periphery = spleen + LN
         total_naiveTregs_thy = Thymic_naiveTregs.host + Thymic_naiveTregs.donor,
         total_memoryTregs_spl = Spleen_memoryTregs.host + Spleen_memoryTregs.donor,
         total_memoryTregs_ln = LN_memoryTregs.host + LN_memoryTregs.donor,
         total_memoryTregs_periph = total_memoryTregs_spl + total_memoryTregs_ln,  ## periphery = spleen + LN
         total_memoryTregs_thy = Thymic_memoryTregs.host + Thymic_memoryTregs.donor,
         ## fd = donor fraction
         fd_naiveTregs_spl = Spleen_naiveTregs.donor/total_naiveTregs_spl,
         fd_naiveTregs_ln = LN_naiveTregs.donor/total_naiveTregs_ln,
         fd_naiveTregs_periph = (Spleen_naiveTregs.donor + LN_naiveTregs.donor)/total_naiveTregs_periph,
         fd_naiveTregs_thy = Thymic_naiveTregs.donor/total_naiveTregs_thy,
         fd_memoryTregs_spl = Spleen_memoryTregs.donor/total_memoryTregs_spl,
         fd_memoryTregs_ln = LN_memoryTregs.donor/total_memoryTregs_ln,
         fd_memoryTregs_periph = (Spleen_memoryTregs.donor + LN_memoryTregs.donor)/total_memoryTregs_periph,
         fd_memoryTregs_thy = Thymic_memoryTregs.donor/total_memoryTregs_thy)  %>%
  filter(mouse.ID != "314807") %>%  ## filtering weird datapoint!
  select(-contains(".host"), -contains(".donor")) 

# normalising donor fraction in splenic naive tregs by dividing with the donor fractions in the source compartment
Treg_fd_Norm <- source_join %>%
  select("mouse.ID", "time.post.BMT", "age.at.S1K", "age.at.BMT", contains('fd')) %>%
  full_join(Treg_join, by = c("mouse.ID", "time.post.BMT", "age.at.S1K", "age.at.BMT"))%>%
  mutate(naiveTregs_Nfd_spl = fd_naiveTregs_spl/ fd_DP1,
         naiveTregs_Nfd_ln = fd_naiveTregs_ln/ fd_DP1,
         naiveTregs_Nfd_periph = fd_naiveTregs_periph/ fd_DP1,
         naiveTregs_Nfd_thy = fd_naiveTregs_thy/ fd_DP1,
         memoryTregs_Chi_spl = fd_memoryTregs_spl/ fd_DP1,
         memoryTregs_Chi_ln = fd_memoryTregs_ln/ fd_DP1,
         memoryTregs_Chi_periph = fd_memoryTregs_periph/ fd_DP1,
         memoryTregs_Chi_thy = fd_memoryTregs_thy/fd_DP1,
         memoryTregs_Nfd_spl = fd_memoryTregs_spl/fd_naiveTregs_spl,
         memoryTregs_Nfd_ln = fd_memoryTregs_ln/fd_naiveTregs_ln,
         memoryTregs_Nfd_periph = fd_memoryTregs_periph/fd_naiveTregs_periph,
         memoryTregs_Nfd_thy = fd_memoryTregs_thy/fd_naiveTregs_thy) %>% 
  filter(memoryTregs_Nfd_periph <= 1.2)%>%
  select(contains("mouse.ID"), contains("time"), contains("age"), contains("Nfd"), contains("Chi")) %>% na.omit()


######### naive Tregs Ki67 proportions in host and donor ######### 
## Percent ki67 in donor naive Tregs 
Treg_donorKi67 <- readxl::read_excel(path = "data/master_doc.xlsx", sheet = 7)  %>% 
  mutate(Ki67_naiveTregs_spl = SP.naiTreg/100,
         Ki67_naiveTregs_ln = LN.naiTreg/100,
         Ki67_naiveTregs_thy = TH.naiTreg/100,
         Ki67_memoryTregs_spl = SP.memTreg/100,
         Ki67_memoryTregs_ln = LN.memTreg/100,
         Ki67_memoryTregs_thy = TH.memTreg/100) %>%
  select(contains("mouse.ID"), contains("time"), contains("age"), contains('Ki67')) %>%
  na.omit() %>% unique() %>%
  left_join(Treg_donor, by = c("mouse.ID", "time.post.BMT", "age.at.BMT", "age.at.S1K"))  %>%
  mutate(Ki67_naiveTregs_periph = (Ki67_naiveTregs_spl * Spleen_naiveTregs + Ki67_naiveTregs_ln * LN_naiveTregs)/
           (Spleen_naiveTregs + LN_naiveTregs),
         Ki67_memoryTregs_periph = (Ki67_memoryTregs_spl * Spleen_memoryTregs + Ki67_memoryTregs_ln * LN_memoryTregs)/
           (Spleen_memoryTregs + LN_memoryTregs))%>%
  select(-contains('fd'),  -contains('spl'), -contains('ln'), -contains('Thymic'))

## Percent ki67 in host naive Tregs 
Treg_hostKi67 <- readxl::read_excel(path = "data/master_doc.xlsx", sheet = 8)  %>% 
  mutate(Ki67_naiveTregs_spl = SP.naiTreg/100,
         Ki67_naiveTregs_ln = LN.naiTreg/100,
         Ki67_naiveTregs_thy = TH.naiTreg/100,
         Ki67_memoryTregs_spl = SP.memTreg/100,
         Ki67_memoryTregs_ln = LN.memTreg/100,
         Ki67_memoryTregs_thy = TH.memTreg/100) %>%
  select(contains("mouse.ID"), contains("time"), contains("age"), contains('Ki67')) %>%
  na.omit() %>% unique() %>%
  left_join(Treg_host, by = c("mouse.ID", "time.post.BMT", "age.at.BMT", "age.at.S1K"))  %>%
  mutate(Ki67_naiveTregs_periph = (Ki67_naiveTregs_spl * Spleen_naiveTregs + Ki67_naiveTregs_ln * LN_naiveTregs)/
           (Spleen_naiveTregs + LN_naiveTregs),
         Ki67_memoryTregs_periph = (Ki67_memoryTregs_spl * Spleen_memoryTregs + Ki67_memoryTregs_ln * LN_memoryTregs)/
           (Spleen_memoryTregs + LN_memoryTregs))%>%
  select(-contains('fd'),  -contains('spl'), -contains('ln'), -contains('Thymic'))

####### PLots ########

## open graphics device to save plots
pdf(file = file.path(getwd(), "figures", paste("Tregs_","Plots%03d.pdf", sep = "")),
    width = 10, height = 4.5, onefile = F)

myTheme <-  theme(axis.text = element_text(size = 14),
      axis.title =  element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 14, face = "bold",  hjust = 0.5),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      strip.text = element_text(size = 14))

Treg_fd_Norm %>%
  select(mouse.ID, time.post.BMT, contains("age"), contains("spl"), contains('ln')) %>%
  gather(-c("mouse.ID", "time.post.BMT", "age.at.S1K", "age.at.BMT"), key = "Pop_of_interest", value = "Nfd") %>%
  mutate(Tissue_location = ifelse(grepl("spl", Pop_of_interest), "Spleen",
                                  ifelse(grepl("ln", Pop_of_interest), "LN",
                                        ifelse(grepl("periph", Pop_of_interest), "Periphery", "Thymus"))),
         pop_nfd = ifelse(grepl("naive", Pop_of_interest), "Naive_normto_DP1",
                                  ifelse(grepl("memoryTregs_Nfd", Pop_of_interest), "Memory_normto_Naive", 
                                         "Memory_normto_DP1"))) %>%
  #filter(pop_nfd != "Naive_normto_DP1") %>%
  ggplot(aes(x = time.post.BMT, y = Nfd)) +
  geom_point(aes(col=Tissue_location), size=2) +
  geom_hline(yintercept = 1.00, linetype = 2, size =1, col=1) +
  #scale_color_viridis_d(name = "Host age at BMT") + ylim(0,1.1)+
  scale_color_manual(name=NULL, values=c(2,4,7,3)) + ylim(0,1.1)+
  #scale_x_log10(breaks=c(75, 150, 300)) + 
  labs(x = "Days post BMT", y = NULL, title = "Normalized donor fractions in Tregs") +
  facet_wrap(.~ pop_nfd) + #+ guides(col='none')+
  theme_bw() + myTheme + theme(legend.position = c(0.95, 0.15), legend.background = element_blank())
  

Treg_join %>%
  select(mouse.ID, time.post.BMT, contains("age"), contains("periph"), contains('thy'), -contains('fd')) %>%
  gather(-c("mouse.ID", "time.post.BMT", "age.at.S1K", "age.at.BMT"), key = "Pop_of_interest", value = "Counts") %>%
  mutate(Tissue_location = ifelse(grepl("spl", Pop_of_interest), "Spleen",
                                  ifelse(grepl("ln", Pop_of_interest), "LN",
                                         ifelse(grepl("periph", Pop_of_interest), "Periphery", "Thymus"))),
         subpop = ifelse(grepl("naive", Pop_of_interest), "Naive", "Memory")) %>%
  ggplot(aes(x = age.at.S1K, y = Counts)) +
  geom_point(aes(col=Tissue_location), size=2) +
  scale_color_manual(values = c(2, 4), name = NULL) +
  scale_y_log10(limits = c(5e3, 1e7), breaks = c(1e4, 1e5, 1e6,1e7)) +
  scale_x_continuous(limits = c(60, 600), trans = "log10", breaks = c(75, 150, 300, 600)) + 
  labs(x = "Host age (days)", y = NULL, title = "Total Treg counts") + 
  facet_wrap(.~ subpop) + #+ guides(col = 'none') +
  theme_bw() + myTheme #+ theme(legend.position = c(0.9, 0.77), legend.background = element_blank())
  
### source fd
source_join %>%
  select(-contains("total"), -contains('ki')) %>%
  rename(DP1 = fd_DP1,  FoxP3_neg_SP4 = fd_FoxP3_neg_SP4, FoxP3_pos_SP4 = fd_FoxP3_pos_SP4) %>%
  gather(-c("mouse.ID", "time.post.BMT", "age.at.S1K", "age.at.BMT"), key = "Pop_of_interest", value = "Chi") %>%
  ggplot(aes(x = time.post.BMT, y = Chi)) +
  geom_point(aes(col=age.at.BMT), size=2) +
  geom_hline(yintercept = 1.00, linetype = 2, size =1, col=1) + ylim(0,1.05) +
  scale_color_viridis_c(name = "Age at BMT") +
  labs(x = "Days post BMT", y = NULL, title = "Donor fractions in thymic populations") +
  facet_wrap(.~ Pop_of_interest, ncol = 3) +
  theme_bw() + myTheme #+ theme(legend.position = c(0.92, 0.87), legend.background = element_blank())

## Source Ki67
source_join %>%
  mutate(DP1_ki.host = TH.DP1_ki.host/100, DP1_ki.donor = TH.DP1_ki.donor/100) %>%
  select(-contains("total"), -contains('fd'), -contains('TH')) %>%
  gather(-c("mouse.ID", "time.post.BMT", "age.at.S1K", "age.at.BMT"), key = "Pop_of_interest", value = "propns_Ki67") %>%
  mutate(subcomp = ifelse(grepl("donor", Pop_of_interest), "Donor", "Host"),
         subpop = ifelse(grepl("DP1", Pop_of_interest), "DP1",
                         ifelse(grepl("neg", Pop_of_interest), "FoxP3_neg_SP4", "FoxP3_pos_SP4"))) %>%
  ggplot(aes(x = time.post.BMT, y = propns_Ki67)) +
  geom_point(aes(col=subcomp), size=2) +
  geom_hline(yintercept = 1.00, linetype = 2, size =1, col=1) +
  scale_color_manual(values = c(7, 2, 4), name = NULL) +
  labs(x = "Days post BMT", y = NULL, title = "Ki67 fractions in thymic populations") +
  facet_wrap(.~ subpop, ncol = 3) +
  theme_bw() + myTheme + theme(legend.position = c(0.92, 0.87), legend.background = element_blank()) 


## Source Counts
source_join %>%
  select(-contains("ki"), -contains('fd'))  %>%
  rename(DP1 = total_DP1,  FoxP3_neg_SP4 = total_FoxP3_neg_SP4, FoxP3_pos_SP4 = total_FoxP3_pos_SP4) %>%
  gather(-c("mouse.ID", "time.post.BMT", "age.at.S1K", "age.at.BMT"), key = "Pop_of_interest", value = "Counts") %>%
  ggplot(aes(x = age.at.S1K, y = Counts)) +
  geom_point(aes(col=age.at.BMT), size=2) +
  scale_color_viridis_c(name = "Age at BMT") +
  #scale_color_manual(values = c(7, 2, 4), name = NULL, labels=c("DP1", "FoxP3- SP4", "FoxP3+ SP4")) +
  scale_y_log10(limits = c(1e4, 2e8), breaks = c(1e4, 1e6, 1e8)) +
  scale_x_continuous(limits = c(60, 500), trans = "log10", breaks = c(75, 150, 300, 600)) + 
  labs(x = "Host age (days)", y = NULL, title = "Total cell counts") +
  facet_wrap(.~ Pop_of_interest) +
  theme_bw() + myTheme #+ theme(legend.position = c(0.12, 0.25), legend.background = element_blank())

Treg_donorKi67 %>%
  left_join(Treg_hostKi67, 
            by = c("mouse.ID", "time.post.BMT", "age.at.S1K", "age.at.BMT"), suffix= c( ".donor", ".host")) %>%
  gather(-c("mouse.ID", "time.post.BMT", "age.at.S1K", "age.at.BMT"), key = "Pop_of_interest", value = "Propn_Ki67") %>%
  mutate(Tissue_location = ifelse(grepl("periph", Pop_of_interest), "Periphery", "Thymus"),
         subcomp = ifelse(grepl("donor", Pop_of_interest), "Donor", "Host"),
         subpop = ifelse(grepl("naive", Pop_of_interest), "Naive", "Memory")) %>%
  ggplot(aes(x = time.post.BMT, y = Propn_Ki67)) +
  geom_point(aes(col=subcomp), size=2) + ylim(0, 0.5) +
  scale_color_manual(values = c(7, 2, 4), name = NULL) +
  labs(x = "Days post BMT", y = NULL, title = NULL) +
  facet_grid(Tissue_location~ subpop) +
  theme_bw() + myTheme + theme(legend.position = c(0.92, 0.9), legend.background = element_blank(),
                               strip.background.y = element_rect(fill="#F783A5"))

Treg_hostKi67 %>%
  rename(Periphery = Ki67_naiveTregs_periph,
         Thymus = Ki67_naiveTregs_thy) %>%
  gather(-c("mouse.ID", "time.post.BMT", "age.at.S1K", "age.at.BMT"), key = "Pop_of_interest", value = "Donor_Ki67") %>%
  ggplot(aes(x = time.post.BMT, y = Donor_Ki67)) +
  geom_point(aes(col=age.at.BMT), size=2) +
  scale_color_viridis_c(name = "Host age at BMT") + ylim(0, 0.4) +
  #scale_color_manual(values = c(7, 2, 4), name = NULL, labels=c("DP1", "FoxP3- SP4", "FoxP3+ SP4")) +
  labs(x = "Days post BMT", y = NULL, title = "Ki67 fractions in host naive Tregs") +
  facet_wrap(.~ Pop_of_interest) +
  theme_bw() + myTheme

dev.off()

#######Importing datafiles as .csv for model fitting precess #######
Treg_join %>%
  mutate(naive = round(total_naiveTregs_periph, 0),
         memory = round(total_memoryTregs_periph, 0)) %>% 
  select(-contains("total"), -contains('fd')) %>%
  write.csv(file = "data/Counts_Treg.csv", row.names = FALSE)

Treg_fd_Norm %>%
  mutate(naive = round(memoryTregs_Nfd_periph, 4),
         memory = round(memoryTregs_Nfd_periph, 4)) %>%
  select(-contains("Nfd"), -contains('memory')) %>%
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
  select(-contains("total"), -contains("fd")) %>%
  write.csv(file = "data/Chimerism_thymicSource.csv", row.names = FALSE)


source_donor %>%
  select(contains("age"), time.post.BMT, contains("neg_SP4_ki")) %>%
  write.csv(file = "data/donorKi67_thymicSource.csv", row.names = FALSE)


source_host %>%
  select(contains("age"), time.post.BMT, contains("neg_SP4_ki")) %>%
  write.csv(file = "data/hostKi67_thymicSource.csv", row.names = FALSE)


