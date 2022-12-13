
## Predictions
init_pred1 <- ode(y=init_cond, times=c(40, 45), func=shm_chi, parms=params, ageatBMT=40)[2,2:13]
init_cond1 <- c(init_pred1[1] + init_pred1[9], init_pred1[2] + init_pred1[10], init_pred1[3] + init_pred1[11],
                init_pred1[4] + init_pred1[12], init_pred1[5], init_pred1[6], init_pred1[7], init_pred1[8],
                y9=0,y10=0,y11=0,y12=0)
R_ode_pred1 <- data.frame(ode(y=init_cond1, times=c(45, ts_pred1), func=shm_chi, parms=params, ageatBMT=45)) %>%
  filter(time != 45) %>%
  mutate(time_seq = time,
         counts_thy = y1 + y2 + y7 + y8 + y9 + y10,
         counts_per = y3 + y4 + y5 + y6 + y11 + y12,
         total_counts = counts_thy + counts_per,
         Nfd_thy = (y9 + y10)/(counts_thy * chivec4),
         Nfd_per = (y11 + y12)/(counts_per * chivec4),
         donor_ki_thy = (y9)/(y9 + y10),
         donor_ki_per = (y11)/(y11 + y12),
         host_ki_thy = (y1 + y7)/(y1 + y2 + y7 + y8),
         host_ki_per = (y3 + y5)/(y3 + y4 + y5 + y6),
         ageBMT_bin = 'agebin1') %>%
  select(time_seq, ageBMT_bin, contains('thy'), contains('per'))

init_pred2 <- ode(y=init_cond, times=c(40, 66), func=shm_chi, parms=params, ageatBMT=40)[2,2:13]
init_cond2 <- c(init_pred2[1] + init_pred2[9], init_pred2[2] + init_pred2[10], init_pred2[3] + init_pred2[11],
                init_pred2[4] + init_pred2[12], init_pred2[5], init_pred2[6], init_pred2[7], init_pred2[8],
                y9=0,y10=0,y11=0,y12=0)

R_ode_pred2 <- data.frame(ode(y=init_cond2,  times=c(66, ts_pred2), func=shm_chi, parms=params, ageatBMT=66)) %>%
  filter(time != 66) %>%
  mutate(time_seq = time,
         counts_thy = y1 + y2 + y7 + y8 + y9 + y10,
         counts_per = y3 + y4 + y5 + y6 + y11 + y12,
         total_counts = counts_thy + counts_per,
         Nfd_thy = (y9 + y10)/(counts_thy * chivec4),
         Nfd_per = (y11 + y12)/(counts_per * chivec4),
         donor_ki_thy = (y9)/(y9 + y10),
         donor_ki_per = (y11)/(y11 + y12),
         host_ki_thy = (y1 + y7)/(y1 + y2 + y7 + y8),
         host_ki_per = (y3 + y5)/(y3 + y4 + y5 + y6),
         ageBMT_bin = 'agebin2') %>%
  select(time_seq, ageBMT_bin, contains('thy'), contains('per'))


init_pred3 <-  ode(y=init_cond, times=c(40, 76), func=shm_chi, parms=params, ageatBMT=40)[2,2:13]
init_cond3 <- c(init_pred3[1] + init_pred3[9], init_pred3[2] + init_pred3[10], init_pred3[3] + init_pred3[11],
                init_pred3[4] + init_pred3[12], init_pred3[5], init_pred3[6], init_pred3[7], init_pred3[8],
                y9=0,y10=0,y11=0,y12=0)

R_ode_pred3 <-data.frame(ode(y=init_cond3,  times=c(76, ts_pred3), func=shm_chi, parms=params, ageatBMT=76)) %>%
  filter(time != 76) %>%
  mutate(time_seq = time,
         counts_thy = y1 + y2 + y7 + y8 + y9 + y10,
         counts_per = y3 + y4 + y5 + y6 + y11 + y12,
         total_counts = counts_thy + counts_per,
         Nfd_thy = (y9 + y10)/(counts_thy * chivec4),
         Nfd_per = (y11 + y12)/(counts_per * chivec4),
         donor_ki_thy = (y9)/(y9 + y10),
         donor_ki_per = (y11)/(y11 + y12),
         host_ki_thy = (y1 + y7)/(y1 + y2 + y7 + y8),
         host_ki_per = (y3 + y5)/(y3 + y4 + y5 + y6),
         ageBMT_bin = 'agebin3') %>%
  select(time_seq, ageBMT_bin, contains('thy'), contains('per'))


init_pred4 <-  ode(y=init_cond, times=c(40, 118), func=shm_chi, parms=params, ageatBMT=40)[2,2:13]
init_cond4 <- c(init_pred4[1] + init_pred4[9], init_pred4[2] + init_pred4[10], init_pred4[3] + init_pred4[11],
                init_pred4[4] + init_pred4[12], init_pred4[5], init_pred4[6], init_pred4[7], init_pred4[8],
                y9=0,y10=0,y11=0,y12=0)

R_ode_pred4 <- data.frame(ode(y=init_cond4,  times=c(118, ts_pred4), func=shm_chi, parms=params, ageatBMT=118)) %>%
  filter(time != 118) %>%
  mutate(time_seq = time,
         counts_thy = y1 + y2 + y7 + y8 + y9 + y10,
         counts_per = y3 + y4 + y5 + y6 + y11 + y12,
         total_counts = counts_thy + counts_per,
         Nfd_thy = (y9 + y10)/(counts_thy * chivec4),
         Nfd_per = (y11 + y12)/(counts_per * chivec4),
         donor_ki_thy = (y9)/(y9 + y10),
         donor_ki_per = (y11)/(y11 + y12),
         host_ki_thy = (y1 + y7)/(y1 + y2 + y7 + y8),
         host_ki_per = (y3 + y5)/(y3 + y4 + y5 + y6),
         ageBMT_bin = 'agebin4') %>%
  select(time_seq, ageBMT_bin, contains('thy'), contains('per'))


rbind(R_ode_pred1, R_ode_pred2, R_ode_pred3, R_ode_pred4)


