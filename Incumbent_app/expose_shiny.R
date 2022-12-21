


## importing data  
counts_file <- file.path("data", "Counts_naiTreg.csv")
counts_data <- read.csv(counts_file) %>% 
  arrange(age.at.S1K) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 56, 'agebin1',
                             ifelse(age.at.BMT <= 70, 'agebin2',
                                    ifelse(age.at.BMT <= 84, 'agebin3', 'agebin4')))) %>%
  gather(c(Thymus, Periphery), key='location', value = 'total_counts')

Nfd_file <- file.path("data", "Nfd_naiTreg.csv")
Nfd_data <- read.csv(Nfd_file) %>% 
  arrange(age.at.S1K)%>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 56, 'agebin1',
                             ifelse(age.at.BMT <= 70, 'agebin2',
                                    ifelse(age.at.BMT <= 84, 'agebin3', 'agebin4')))) %>%
  gather(c(Thymus, Periphery), key='location', value = 'Nfd')

hostki_file <- file.path("data", "hostKi67_naiTreg.csv")
hostki_data <- read.csv(hostki_file) %>% 
  arrange(age.at.S1K) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 56, 'agebin1',
                             ifelse(age.at.BMT <= 70, 'agebin2',
                                    ifelse(age.at.BMT <= 84, 'agebin3', 'agebin4'))),
         subcomp='Host') %>%
  gather(c(Thymus, Periphery), key='location', value = 'prop_ki')

donorki_file <- file.path("data", "donorKi67_naiTreg.csv")
donorki_data <- read.csv(donorki_file) %>% 
  arrange(age.at.S1K) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 56, 'agebin1',
                             ifelse(age.at.BMT <= 70, 'agebin2',
                                    ifelse(age.at.BMT <= 84, 'agebin3', 'agebin4'))),
         subcomp='Donor')%>%
  gather(c(Thymus, Periphery), key='location', value = 'prop_ki')

ki_data <- rbind(donorki_data, hostki_data)

### R ode solver predictions
Theta_spline <- function(Time, psi){psi * (10^6.407133491 * exp(-0.002387866 * (Time - 49)))}
chi_spline <- function(Time){ifelse(Time-10 >0, 0.81548689 * (1 - exp(-0.06286984 * (Time - 10))), 0)}
Donor_eps_spline <- function(Time){exp(- 0.06799028 * (Time - 49)) + 0.37848471}

## ODE System
shm_chi <- function (Time, ageatBMT, y, parms) {
  eps_host = 0.326611
  kloss = 1/3.5
  with(as.list(c(y, parms)),{
    dy1 <- Theta_spline(Time, psi) * (1- chi_spline(Time - ageatBMT)) * eps_host + rho_D * (2 * y2 + y1) + beta * y3 - (kloss + alpha + delta_D) * y1
    dy2 <- Theta_spline(Time, psi) * (1- chi_spline(Time - ageatBMT)) * (1 - eps_host) + kloss * y1 + beta * y4  - (rho_D + alpha + delta_D) * y2
    dy3 <- alpha * y1 + rho_D * (2 * y4 + y3) - (kloss + beta + delta_D) * y3
    dy4 <- alpha * y2 + kloss * y3 - (rho_D + beta + delta_D) * y4
    dy5 <- alpha * y7 + rho_I * (2 * y6 + y5) - (kloss + beta + rho_I) * y5
    dy6 <- alpha * y8 + kloss * y5 - (rho_I + beta + rho_I) * y6
    dy7 <- beta * y5 + rho_I * (2 * y8 + y7) - (kloss + alpha + rho_I) * y7
    dy8 <- beta * y6 + kloss * y7 - (rho_I + alpha + rho_I) * y8
    
    dy9 <- Theta_spline(Time, psi) * chi_spline(Time - ageatBMT) * Donor_eps_spline(Time)+ rho_D * (2 * y10 + y9) + beta * y11  - (kloss + alpha + delta_D) * y9
    dy10 <- Theta_spline(Time, psi) * chi_spline(Time - ageatBMT) * (1 - Donor_eps_spline(Time)) + kloss * y9 + beta * y12  - (rho_D + alpha + delta_D) * y10
    dy11 <- alpha * y9 + rho_D * (2 * y12 + y11) - (kloss + beta + delta_D) * y11
    dy12 <- alpha * y10 + kloss * y11 - (rho_D + beta + delta_D) * y12
    
    list(c(dy1, dy2, dy3, dy4, dy5, dy6, dy7, dy8,
           dy9, dy10, dy11, dy12))
  })
}

# time sequence for predictions specific to age bins within the data
ts_pred1 <- 10^seq(log10(66), log10(450), length.out = 300)
ts_pred2 <- 10^seq(log10(91), log10(450), length.out = 300)
ts_pred3 <- 10^seq(log10(90), log10(450), length.out = 300)
ts_pred4 <- 10^seq(log10(174), log10(450), length.out = 300)
tb_pred1 <- rep(45, 300)
tb_pred2 <- rep(66, 300)
tb_pred3 <- rep(76, 300)
tb_pred4 <- rep(118, 300)


chivec1 <- sapply(ts_pred1-tb_pred1, chi_spline)
chivec2 <- sapply(ts_pred2-tb_pred2, chi_spline)
chivec3 <- sapply(ts_pred3-tb_pred3, chi_spline)
chivec4 <- sapply(ts_pred4-tb_pred4, chi_spline)

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
fac_labels <- c(`agebin1`= '6-8 weeks', `agebin2`= '8-10 weeks', `agebin3`= '10-12 weeks', `agebin4`= '12-25 weeks')



