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
thynaiveTregs_counts <- read_csv("data/Counts_naiTreg.csv") 
thynaiveTregs_chimerism <- read_csv("data/Nfd_naiTreg.csv") %>% arrange(time.post.BMT)
thynaiveTregs_donorki <- read_csv("data/donorKi67_naiTreg.csv") 
thynaiveTregs_hostki <- read_csv("data/hostKi67_naiTreg.csv") 


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

(chi_spline(thynaiveTregs_chimerism$time.post.BMT, 0.8, 0.01))
             
## LL function
logit_transf <- function(x){asin(sqrt(x))}
Thy_chi_nlm <- nls((Thymus) ~ (chi_spline(time.post.BMT, chiEst, qEst)),
                  data =  thynaiveTregs_chimerism,
                  start = list(chiEst=0.84, qEst=0.01))
Thy_chi_pars <- coef(Thy_chi_nlm)

# prediction
Thy_chi_fit <- data.frame(dpt_seq, "y_sp" = chi_spline(dpt_seq, Thy_chi_pars[1], Thy_chi_pars[2]))
Thy_chi_fit_m <- data.frame(dpt_seq, "y_sp" = chi_spline(dpt_seq, 0.8, 0.065))
ggplot() + 
  geom_point(data= thynaiveTregs_chimerism, aes(x=age.at.S1K-age.at.BMT, y=Thymus), size =2) +
  geom_line(data = Thy_chi_fit, aes(x = dpt_seq , y = y_sp), col=4, size =1) + 
  #geom_line(data = Thy_chi_fit_m, aes(x = dpt_seq , y = y_sp), col=2, size =1) + 
  ylim(0,1) + #scale_x_log10(limits=c(1, 350)) +
  labs(title = 'Donor fraction in FoxP3 positive naive SP4 cells',  y=NULL,  x = 'Days post BMT') 

## Counts fit
## phenomenological function
counts_spline <- function(age.at.S1K, b1, b0, nu){
  Time = age.at.S1K   ###
  t0 = 40
  return(10^b1 + (10^b0/(1+ ((Time-t0)/nu)^2)))
}

countsvec <- sapply(timeseq, counts_spline, b1=4, b0=5, nu=120)
ggplot()+
  geom_line(aes(x=timeseq, y= countsvec))

thy_counts_nlm <- nls(log(Thymus) ~ log(counts_spline(age.at.S1K, b0, nu)),
                   data =  thynaiveTregs_counts,
                   start = list(b0=5, nu=10))
thy_counts_pars <- coef(thy_counts_nlm)

thycounts_fit <- data.frame(timeseq, "y_sp" = counts_spline(timeseq, thy_counts_pars[1], thy_counts_pars[2]))
thycounts_fit_m <- data.frame(timeseq, "y_sp" = counts_spline(timeseq,4.3, 5.1, 40))

ggplot() + 
  geom_point(data=thynaiveTregs_counts, aes(x=age.at.S1K, y=Thymus), col=4, size =2) +
  #geom_line(data = thycounts_fit, aes(x = timeseq , y = y_sp), col=4, size =1) + 
  geom_line(data = thycounts_fit_m, aes(x = timeseq , y = y_sp), col=4, size =1) + 
  scale_y_log10(limits=c(1e3, 2e5)) + xlim(40, 450)+
  labs(title = 'Total counts of thymic FoxP3 positive naive SP4 cells',  y=NULL,  x = 'Host age (days)') 


thynaiveTregs_hostki$Periphery %>% mean()
thynaiveTregs_donorki$Periphery %>% mean()


ki_init <- function(ki){
  r_ki_init=exp(-1)
  if(ki >= 0.0 && ki <= 1.0){
    value = exp(-ki * r_ki_init)/((1 - exp(-r_ki_init))/r_ki_init);
  }  else {
    value = 0.0;
  }
}

kseq <- seq(0, 1, 0.001)
kivec <- sapply(kseq, ki_init)
plot(kivec ~ kseq)

# exporting plot as PDF
pdf(file = file.path("out_fit", paste(colnames(cnts)[4], "IndPlots%03d.pdf", sep = "")),
    width = 12, height = 4.5, onefile = FALSE, useDingbats = FALSE )

cowplot::plot_grid(p2, p1, nrow = 1, labels = c("A", "B"))

dev.off()

