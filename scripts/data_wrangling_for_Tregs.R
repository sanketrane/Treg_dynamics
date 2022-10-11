### data wrangling for T regulatory cells from Busulfan chimera experiments

## claning environment and cache
rm(list = ls()); gc()

#Loading required libraries
library(tidyverse)

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


# total counts and donor fractions for the source poppulation
source_donor <- readxl::read_excel(path = "data/master_doc.xlsx", sheet = 2) %>%
  select(contains("mouse.ID"), contains("time"), contains("age"), contains("DP1"), contains("Fox25"))%>% 
  na.omit() %>% unique() %>%
  mutate(FoxP3_pos_SP4 = TH.Fox25Q2 + TH.Fox25Q3,  ### 2nd and 3r quadrants of FOXP3 (x-axis) and CD25 (y-axis) Boolean gate
         FoxP3_neg_SP4 = TH.Fox25Q1 + TH.Fox25Q4,
         ratioQ1_Q4 =  TH.Fox25Q1/TH.Fox25Q4) %>%   
  select(-contains("Fox25"))

source_host <- readxl::read_excel(path = "data/master_doc.xlsx", sheet = 3) %>%
  select(contains("mouse.ID"), contains("time"), contains("age"), contains("DP1"), contains("Fox25"))%>% 
  na.omit() %>% unique() %>%
  mutate(FoxP3_pos_SP4 = TH.Fox25Q2 + TH.Fox25Q3, ### 2nd and 3r quadrants of FOXP3 (x-axis) and CD25 (y-axis) Boolean gate
         FoxP3_neg_SP4 = TH.Fox25Q1 + TH.Fox25Q4,
         ratioQ1_Q4 =  TH.Fox25Q1/TH.Fox25Q4) %>%
  select(-contains("Fox25"))

# merging total counts for host and donor compartments
# calculating total counts, donor fractions
source_join <- full_join(source_host, source_donor, by = c("mouse.ID", "time.post.BMT", "age.at.S1K", "age.at.BMT"), suffix= c(".host", ".donor")) %>%
  mutate(# total = donor + host
         total_DP1 = TH.DP1.host + TH.DP1.donor,
         total_FoxP3_pos_SP4 = FoxP3_pos_SP4.host + FoxP3_pos_SP4.donor,
         total_FoxP3_neg_SP4 = FoxP3_neg_SP4.host + FoxP3_neg_SP4.donor,
         ## fd = donor fraction
         fd_DP1 = TH.DP1.donor/total_DP1,
         fd_FoxP3_pos_SP4 = FoxP3_pos_SP4.donor/total_FoxP3_pos_SP4,
         fd_FoxP3_neg_SP4 = FoxP3_neg_SP4.donor/total_FoxP3_neg_SP4)%>%
  select(-contains(".host"), -contains(".donor"))


# normalising donor fraction in splenic naive tregs by dividing with the donor fractions in the source compartment
Treg_fd_Norm <- Treg_join %>%
  full_join(source_join, by = c("mouse.ID", "time.post.BMT", "age.at.S1K", "age.at.BMT"))%>%
  mutate(naiveTregs_Nfd_spl = fd_naiveTregs_spl/ fd_FoxP3_neg_SP4,
         naiveTregs_Nfd_ln = fd_naiveTregs_ln/ fd_FoxP3_neg_SP4,
         naiveTregs_Nfd_periph = fd_naiveTregs_periph/ fd_FoxP3_neg_SP4,
         naiveTregs_Nfd_thy = fd_naiveTregs_thy/ fd_FoxP3_neg_SP4,
         memoryTregs_Chi_spl = fd_memoryTregs_spl/ fd_FoxP3_neg_SP4,
         memoryTregs_Chi_ln = fd_memoryTregs_ln/ fd_FoxP3_neg_SP4,
         memoryTregs_Chi_periph = fd_memoryTregs_periph/ fd_FoxP3_neg_SP4,
         memoryTregs_Chi_thy = fd_memoryTregs_thy/fd_FoxP3_neg_SP4,
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
         Ki67_naiveTregs_thy = TH.naiTreg/100) %>%
  select(contains("mouse.ID"), contains("time"), contains("age"), contains('naive')) %>%
  na.omit() %>% unique() %>%
  left_join(Treg_donor, by = c("mouse.ID", "time.post.BMT", "age.at.BMT", "age.at.S1K"))  %>%
  mutate(Ki67_naiveTregs_periph = (Ki67_naiveTregs_spl * Spleen_naiveTregs + Ki67_naiveTregs_ln * LN_naiveTregs)/
           (Spleen_naiveTregs + LN_naiveTregs))%>%
  select(-contains('fd'),  -contains('spl'), -contains('ln'), -contains('Thymic'))

## Percent ki67 in host naive Tregs 
Treg_hostKi67 <- readxl::read_excel(path = "data/master_doc.xlsx", sheet = 8)  %>% 
  mutate(Ki67_naiveTregs_spl = SP.naiTreg/100,
         Ki67_naiveTregs_ln = LN.naiTreg/100,
         Ki67_naiveTregs_thy = TH.naiTreg/100) %>%
  select(contains("mouse.ID"), contains("time"), contains("age"), contains('naive')) %>%
  na.omit() %>% unique() %>%
  left_join(Treg_host, by = c("mouse.ID", "time.post.BMT", "age.at.BMT", "age.at.S1K"))  %>%
  mutate(Ki67_naiveTregs_periph = (Ki67_naiveTregs_spl * Spleen_naiveTregs + Ki67_naiveTregs_ln * LN_naiveTregs)/
           (Spleen_naiveTregs + LN_naiveTregs))%>%
  select(-contains('fd'),  -contains('spl'), -contains('ln'), -contains('Thymic'))

source_donorKi67 <- readxl::read_excel(path = "data/master_doc.xlsx", sheet = 7) %>%
  select(contains("mouse.ID"), contains("time"), contains("age"), contains("DP1"), contains("Fox25"))%>% 
  na.omit() %>% unique() %>%
  left_join(source_donor, by = c("mouse.ID", "time.post.BMT", "age.at.BMT", "age.at.S1K")) %>%
  select(contains("mouse.ID"), contains("time"), contains("age"), contains("Q1"), contains("Q4")) %>%
  mutate(FoxP3_Neg_SP4 = ((TH.Fox25Q1/100) * ratioQ1_Q4 + (TH.Fox25Q4/100))/(1 + ratioQ1_Q4))  %>%   
  select(-contains("Q1"), -contains("Q4"))


source_hostKi67 <- readxl::read_excel(path = "data/master_doc.xlsx", sheet = 8) %>%
  select(contains("mouse.ID"), contains("time"), contains("age"), contains("DP1"), contains("Fox25"))%>% 
  na.omit() %>% unique() %>%
  left_join(source_host, by = c("mouse.ID", "time.post.BMT", "age.at.BMT", "age.at.S1K")) %>%
  select(contains("mouse.ID"), contains("time"), contains("age"), contains("Q1"), contains("Q4")) %>%
  mutate(FoxP3_Neg_SP4 = ((TH.Fox25Q1/100) * ratioQ1_Q4 + (TH.Fox25Q4/100))/(1 + ratioQ1_Q4))  %>%   
  select(-contains("Q1"), -contains("Q4"))



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

Treg_fd_Norm  %>%
  select(mouse.ID, time.post.BMT, contains("age"), contains("periph"), contains('thy')) %>%
  gather(-c("mouse.ID", "time.post.BMT", "age.at.S1K", "age.at.BMT"), key = "Pop_of_interest", value = "Nfd") %>%
  mutate(Tissue_location = ifelse(grepl("spl", Pop_of_interest), "Spleen",
                                  ifelse(grepl("ln", Pop_of_interest), "LN",
                                        ifelse(grepl("periph", Pop_of_interest), "Periphery", "Thymus"))),
         pop_nfd = ifelse(grepl("naive", Pop_of_interest), "Naive",
                                  ifelse(grepl("memoryTregs_Nfd", Pop_of_interest), "Memory_normto_Naive", 
                                         "Memory_normto_SP4"))) %>%
  ggplot(aes(x = time.post.BMT, y = Nfd)) +
  geom_point(aes(col=Tissue_location), size=2) +
  geom_hline(yintercept = 1.00, linetype = 2, size =1, col=1) +
  #scale_color_viridis_d(name = "Host age at BMT") + ylim(0,1.1)+
  scale_color_manual(name=NULL, values=c(2,4,7, 3)) + ylim(0,1.1)+
  #scale_x_log10(breaks=c(75, 150, 300)) + 
  labs(x = "Days post BMT", y = NULL, title = "Normalized donor fractions in Tregs") +
  facet_wrap(.~ pop_nfd) + #+ guides(col='none')+
  theme_bw() + myTheme + theme(legend.position = c(0.95, 0.15), legend.background = element_blank())
  

Treg_join %>%
  select(-contains('fd')) %>%
  rename(Periphery = total_naiveTregs_periph,
         Thymus = total_naiveTregs_thy) %>%
  gather(-c("mouse.ID", "time.post.BMT", "age.at.S1K", "age.at.BMT"), key = "Pop_of_interest", value = "Counts") %>%
  filter(Pop_of_interest !=  "total_naiveTregs_spl",
         Pop_of_interest !=  "total_naiveTregs_ln") %>%
  ggplot(aes(x = age.at.S1K, y = Counts)) +
  geom_point(aes(col=age.at.BMT), size=2) +
  scale_color_viridis_c(name = "Host age at BMT") +
  scale_y_log10(limits = c(5e3, 1e7), breaks = c(1e4, 1e5, 1e6,1e7)) +
  scale_x_continuous(limits = c(60, 600), trans = "log10", breaks = c(75, 150, 300, 600)) + 
  labs(x = "Host age (days)", y = NULL, title = "Total cell counts") + 
  facet_wrap(.~ Pop_of_interest) + #+ guides(col = 'none') +
  theme_bw() + myTheme + theme(legend.position = c(0.9, 0.77), legend.background = element_blank())
  

source_join %>%
  rename(DP1 = fd_DP1,
         FoxP3_Pos_SP4 = fd_FoxP3_pos_SP4,
         FoxP3_Neg_SP4 = fd_FoxP3_neg_SP4) %>%
  select(-contains("total"), -"FoxP3_Pos_SP4") %>%
  gather(-c("mouse.ID", "time.post.BMT", "age.at.S1K", "age.at.BMT"), key = "Pop_of_interest", value = "Chi") %>%
  ggplot(aes(x = time.post.BMT, y = Chi)) +
  geom_point(aes(col=age.at.BMT), size=2) +
  geom_hline(yintercept = 1.00, linetype = 2, size =1, col=1) +
  scale_color_viridis_c(name = "Host age at BMT") +
  #scale_color_manual(values = c(7, 2, 4), name = NULL, labels=c("DP1", "FoxP3- SP4", "FoxP3+ SP4")) +
  labs(x = "Days post BMT", y = NULL, title = "Donor fractions in Thymic precursors") +
  facet_wrap(.~ Pop_of_interest) +
  theme_bw() + myTheme #+ theme(legend.position = c(0.9, 0.77), legend.background = element_blank())

source_join %>%
  rename(DP1 = total_DP1,
         FoxP3_Pos_SP4 = total_FoxP3_pos_SP4,
         FoxP3_Neg_SP4 = total_FoxP3_neg_SP4) %>%
  select(-contains("fd"), -"FoxP3_Pos_SP4") %>%
  gather(-c("mouse.ID", "time.post.BMT", "age.at.S1K", "age.at.BMT"), key = "Pop_of_interest", value = "Chi") %>%
  ggplot(aes(x = age.at.S1K, y = Chi)) +
  geom_point(aes(col=age.at.BMT), size=2) +
  scale_color_viridis_c(name = "Host age at BMT") +
  #scale_color_manual(values = c(7, 2, 4), name = NULL, labels=c("DP1", "FoxP3- SP4", "FoxP3+ SP4")) +
  scale_y_log10(limits = c(1e5, 1e8), breaks = c(1e4, 1e5, 1e6,1e7)) +
  scale_x_continuous(limits = c(60, 600), trans = "log10", breaks = c(75, 150, 300, 600)) + 
  labs(x = "Host age (days)", y = NULL, title = "Total cell counts") +
  facet_wrap(.~ Pop_of_interest) +
  theme_bw() + myTheme #+ theme(legend.position = c(0.12, 0.25), legend.background = element_blank())

Treg_donorKi67 %>%
  rename(Periphery = Ki67_naiveTregs_periph,
         Thymus = Ki67_naiveTregs_thy) %>%
  gather(-c("mouse.ID", "time.post.BMT", "age.at.S1K", "age.at.BMT"), key = "Pop_of_interest", value = "Donor_Ki67") %>%
  ggplot(aes(x = time.post.BMT, y = Donor_Ki67)) +
  geom_point(aes(col=age.at.BMT), size=2) +
  scale_color_viridis_c(name = "Host age at BMT") + ylim(0, 0.4) +
  #scale_color_manual(values = c(7, 2, 4), name = NULL, labels=c("DP1", "FoxP3- SP4", "FoxP3+ SP4")) +
  labs(x = "Days post BMT", y = NULL, title = "Ki67 fractions in donor naive Tregs") +
  facet_wrap(.~ Pop_of_interest) +
  theme_bw() + myTheme

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
  mutate(Periphery = round(total_naiveTregs_periph, 0),
         Thymus = round(total_naiveTregs_thy, 0)) %>% 
  select(-contains("total"), -contains('fd')) %>%
  write.csv(file = "data/Counts_naiTreg.csv", row.names = FALSE)

Treg_fd_Norm %>%
  mutate(Periphery = round(Nfd_naiveTregs_periph, 4),
         Thymus = round(Nfd_naiveTregs_thy, 4)) %>%
  select(-contains("Nfd")) %>%
  write.csv(file = "data/Nfd_naiTreg.csv", row.names = FALSE)

Treg_hostKi67 %>%
  mutate(Periphery = round(Ki67_naiveTregs_periph, 4),
         Thymus = round(Ki67_naiveTregs_thy, 4)) %>%
  write.csv(file = "data/hostKi67_naiTreg.csv", row.names = FALSE)

Treg_donorKi67 %>%
  mutate(Periphery = round(Ki67_naiveTregs_periph, 4),
         Thymus = round(Ki67_naiveTregs_thy, 4)) %>%
  write.csv(file = "data/donorKi67_naiTreg.csv", row.names = FALSE)
  
source_join %>%
  mutate(FoxP3_Neg_SP4 = round(total_FoxP3_neg_SP4, 0))%>%
  select(-contains("total"), -contains("fd")) %>%
  write.csv(file = "data/Counts_thymicSource.csv", row.names = FALSE)

source_join %>% 
  mutate(FoxP3_Neg_SP4 = round(fd_FoxP3_neg_SP4, 4)) %>%
  select(-contains("total"), -contains("fd")) %>%
  write.csv(file = "data/Chimerism_thymicSource.csv", row.names = FALSE)


source_donorKi67 %>%
  write.csv(file = "data/donorKi67_thymicSource.csv", row.names = FALSE)


source_hostKi67 %>%
  write.csv(file = "data/hostKi67_thymicSource.csv", row.names = FALSE)


