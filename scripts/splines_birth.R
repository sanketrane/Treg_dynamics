### data wrangling for T regulatory cells from Busulfan chimera experiments

## claning environment and cache
rm(list = ls()); gc()

#Loading required libraries
library(tidyverse)

## Precursor population dynamics

#### Putative precursors for naive Tregs: Thymic DP1, Thymic FoxP3negative SP4
#### Putative precursors for memory Tregs: Thymic DP1/Thymic FoxP3negative SP4 (representing RTE), naive Tregs

### Ki67 Proportions
source_WTKi67 <- readxl::read_excel(path = "data/ontogeny_data.xlsx", sheet = 3) %>%
  select(contains("mouse.id"), contains("age.at.S1K.days"), contains("DP1"), contains("Fox25"))%>% 
  na.omit() %>% unique() 

### Total counts and donor fractions for the source population
source_data<- readxl::read_excel(path = "data/ontogeny_data.xlsx", sheet = 2) %>%
  select(contains("mouse.id"), contains("time"), contains("age.at.S1K.days"), contains("DP1"), contains("Fox25")) %>% 
  left_join(source_WTKi67, by = c("mouse.id", "age.at.S1K.days"), suffix= c("_counts", "_ki")) %>%
  na.omit() %>% unique() %>%
  mutate(### 1st and 4th quadrants of FOXP3 (x-axis) and CD25 (y-axis) Boolean gate, FoxP3- SP4 = Treg free SP4
    FoxP3_neg_SP4_counts = Fox25.Q1_counts + Fox25.Q4_counts,
    ### 2nd and 3rd quadrants of FOXP3 (x-axis) and CD25 (y-axis) Boolean gate, FoxP3+ SP4 = SP4 Tregs
    FoxP3_pos_SP4_counts = Fox25.Q2_counts + Fox25.Q3_counts,
    FoxP3_neg_SP4_ki =  (((Fox25.Q1_ki/100) * Fox25.Q1_counts) + ((Fox25.Q4_ki/100) * Fox25.Q4_counts))/(Fox25.Q1_counts + Fox25.Q4_counts),
    FoxP3_pos_SP4_ki =  (((Fox25.Q2_ki/100) * Fox25.Q2_counts) + ((Fox25.Q3_ki/100) * Fox25.Q3_counts))/(Fox25.Q2_counts + Fox25.Q3_counts))%>%   
  select(-contains("Fox25")) 

ggplot(source_data)+
  geom_point(aes(x=age.at.S1K.days, y= FoxP3_pos_SP4_counts))+
  scale_y_log10() + scale_x_log10()

ggplot(source_data)+
  geom_point(aes(x=age.at.S1K.days, y= FoxP3_pos_SP4_ki)) + scale_x_log10()

df_merge <- source_data %>%
  select(mouse.id, age.at.S1K.days, contains("FoxP3_pos")) %>%
  rename(age.at.S1K = age.at.S1K.days,
         total_counts = FoxP3_pos_SP4_counts,
         ki_prop = FoxP3_pos_SP4_ki) %>% arrange(age.at.S1K)

## Counts fit
## phenomenological function
count_spline <- function(Time, basl, theta, n, X, q){
  t = Time - 5
  return(exp(basl) + (theta * t^n) * (1 - ((t^q)/((X^q) + (t^q)))))
}

ki_spline <- function(Time, eps_0, eps_f, A){
  
  return(exp(- eps_f * (Time + A)) + eps_0)
}

## LL function for counts
count_obj <- function(params) {
  basl <- params[["basl"]]
  theta  <- params[["theta"]]
  n  <- params[["n"]]
  X  <- params[["X"]]
  q   <- params[["q"]]
  
  ## SSR
  df_merge %>%
    mutate(count_fit = count_spline(age.at.S1K,
                                    basl = basl,
                                    theta = theta,
                                    n = n,
                                    X = X,
                                    q = q),
           log_sr = (log(total_counts) - log(count_fit)) ^ 2) %>%
    summarize(log_ssr = sum(log_sr)) %>%
    unlist()
}

# initial guess for LL calculation
params <- list(basl = 12.5, theta = 50, n = 2 , X = 25, q = 2)

## model fit
counts_fit <- optim(params, count_obj, control = list(trace = TRUE, maxit = 1000))
counts_pars_est <- counts_fit$par

## LL function for ki proportions
ki_obj <- function(params) {
  eps_0 <- params[["eps_0"]]
  eps_f  <- params[["eps_f"]]
  A  <- params[["A"]]
  
  ## SSR
  df_merge %>%
    mutate(ki_fit = ki_spline(age.at.S1K,
                              eps_0 = eps_0,
                              eps_f = eps_f,
                              A = A),
           SqRes = ((ki_prop) - (ki_fit)) ^ 2) %>%
    summarize(SSR = sum(SqRes)) %>%
    unlist()
}

# initial guess for LL calculation
params <- list(eps_0 = 0.1, eps_f = 0.03, A = 3)

## model fit
ki_fit <- optim(params, ki_obj, control = list(trace = TRUE, maxit = 1000))
ki_pars_est <- ki_fit$par

## time course for prediction
ts_p <- seq(5, 350, 0.1)

myTheme <- theme(text = element_text(size = 12), axis.text = element_text(size = 12),
                 axis.title =  element_text(size = 12, face = "bold"),
                 plot.title = element_text(size=12,  hjust = 0.5, face = "bold"),
                 legend.background = element_blank(), legend.key = element_blank())

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

log10minorbreaks = as.numeric(1:10 %o% 10^(4:8))


## prediction plots
ki_fit_spline <- data.frame(ts_p, "y_sp" = ki_spline(ts_p, ki_pars_est[1], ki_pars_est[2], ki_pars_est[3]))

ggplot(df_merge) + 
  geom_point(aes(age.at.S1K, ki_prop), col = "#0096b4", size =2) + 
  geom_line(data = ki_fit_spline, aes(x = ts_p, y = y_sp), col = "#0096b4", size =1) + 
  scale_x_log10(limits= c(3, 350), breaks = c(3,10,30,100, 300)) + scale_y_continuous(limits = c(0, 1)) +
  labs(title = 'Ki67+ proprtions in SP4 T cells',  y=NULL,  x = 'Host age (days)') +
  myTheme + guides(col="none")


ggsave(filename = file.path("out_fit",  "Ki67_sp4.pdf"), last_plot(), device = "pdf", width = 4.8, height = 3.45)

## prediction plots
counts_fit_spline2 <- data.frame(ts_p, "y_sp" = count_spline(ts_p, counts_pars_est[1], counts_pars_est[2], counts_pars_est[3], counts_pars_est[4], counts_pars_est[5]))
counts_fit_spline <- data.frame(ts_p, "y_sp" = count_spline(ts_p, 10.5, 150, 3, 12, 4))

ggplot(df_merge) + geom_point(aes(age.at.S1K, total_counts), col = "#0096b4", size =2) + 
  geom_line(data = counts_fit_spline, aes(x = ts_p, y = y_sp), col = 2, size =1) + 
  geom_line(data = counts_fit_spline2, aes(x = ts_p, y = y_sp), col = "#0096b4",  size =1) + 
  scale_x_log10(limits= c(3, 350), breaks = c(3,10,30,100, 300)) + #scale_y_log10()+
  labs(y =NULL, x = 'Host age (days)', title =  'Total counts in SP4 T cells')+
  scale_y_continuous(trans="log10",minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  myTheme + guides(col="none")


ggsave(filename = file.path("out_fit", "counts_sp4.pdf"), last_plot(), device = "pdf", width = 4.8, height = 3.45)






