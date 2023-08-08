functions{
 //spline 1
 // Timecourse of thymic CD4 SP population -- changes with time
 real sp_numbers(real time) {
   real t0 = 49.0;     // mean mouse age at BMT for the first ageBMT bin
   real value;
   // spline fitted separately to the counts of thymic FoxP3 negative SP4 T cells to estimate the parameters
   real basl  = 6.407133491;  real nu = 0.002387866 ;
   //best fitting spline
   value = 10^basl * exp(-nu * (time - t0));
   return value;
 }

 // Total influx into the naive TReg cell compartment from Foxp negative Sp4 T cells (cells/day)
 real theta_spline(real time, real psi){
   //psi is the proportionality constant -- per capita rate of influx
   real value = psi * sp_numbers(time);
   return value;
 }

 // spline2 --
 real Chi_spline(real time) {
   // chiEst is the level of stabilised chimerism in the source (FoxP3 negative SP4) compartment
   // qEst is the rate with which cimerism changes in the source (FoxP3 negative SP4) compartment
   // spline fitted separately to the donor chimerism in the thymic FoxP3 negative SP4 T cells to estimate the parameters
   real value;  real chiEst = 0.81548689;   real qEst = 0.06286984;
   if (time - 10 < 0){              // t0 = 14 -- assumption: no donor cells seen in FoxP3neg SP4 compartment for 2 weeks
     value = 0;
   } else {
     value = chiEst * (1 - exp(-qEst * (time - 10)));
   }
   return value;
 }

 // spline3 --
 // proportions of ki67hi cells in the donor-derived FoxP3 negative SP4 T cells -- varies with time
 real donor_eps_spline(real time){
   real t0 = 49.0;     // mean mouse age at BMT for the first ageBMT bin
  //parameters estimated from spline fit to the timecourse of ki67 fraction in the donor-derived FoxP3 negative SP4 T cells
  real eps_0 = 0.37848471; real eps_f = 0.06799028;
  return 0.0; //exp(- eps_f * (time - t0)) + eps_0;
}

real[] shm_chi(real time, real[] y, real[] parms, real[] rdata,  int[] idata) {
  real psi = parms[1];
  real rho_D = parms[2];
  real alpha = parms[3];
  real delta_D = parms[4];
  real rho_I = parms[5];
  real Beta = parms[6];

  real dydt[12];
  real kloss  = 1/3.5;            //rate of loss of ki67
  real eps_host = 0.0; //0.326611;      // Mean Ki67 hi fraction in host-BM-derived FoxP3 negative Sp4 T cells

  // age of BMT in each recipient
  real ageAtBMT = parms[7];

  // model that assumes that tranistionals divide and die at different rates than mature naive T cells
  // Host naive Tregs
  // Thymic ki  hi displaceable
  //dydt[1] = theta_spline(time, psi) * (1- Chi_spline(time - ageAtBMT)) * eps_host + rho_D * (2 * y[2] + y[1]) + Beta * y[3] - (kloss + alpha + delta_D) * y[1];
  dydt[1] = theta_spline(time, psi) * (1- Chi_spline(time - ageAtBMT)) * eps_host + rho_D * (2 * y[2] + y[1]) + Beta * y[3] - (kloss + alpha + delta_D) * y[1];
  // Thymic ki lo displaceable
  dydt[2] = theta_spline(time, psi) * (1- Chi_spline(time - ageAtBMT)) * (1 - eps_host) + kloss * y[1] + Beta * y[4]  - (rho_D + alpha + delta_D) * y[2];
  // Peripheral ki hi displaceable
  dydt[3] = alpha * y[1] + rho_D * (2 * y[4] + y[3]) - (kloss + Beta + delta_D) * y[3];
  // Peripheral ki lo displaceable
  dydt[4] = alpha * y[2] + kloss * y[3] - (rho_D + Beta + delta_D) * y[4];
  // Peripheral ki hi Incumbent
  dydt[5] = alpha * y[7] + rho_I * (2 * y[6] + y[5]) - (kloss + Beta + rho_I) * y[5];
  // Peripheral ki lo Incumbent
  dydt[6] = alpha * y[8] + kloss * y[5] - (rho_I + Beta + rho_I) * y[6];
  // Thymic ki hi Incumbent
  dydt[7] = Beta * y[5] + rho_I * (2 * y[8] + y[7]) - (kloss + alpha + rho_I) * y[7];
  // Thymic ki lo Incumbent
  dydt[8] = Beta * y[6] + kloss * y[7] - (rho_I + alpha + rho_I) * y[8];

  // Donor naive Tregs
  // Thymic ki  hi displaceable
  //dydt[9] = theta_spline(time, psi) * Chi_spline(time - ageAtBMT) * donor_eps_spline(time) + rho_D * (2 * y[10] + y[9]) + Beta * y[11]  - (kloss + alpha + delta_D) * y[9];
  dydt[9] = theta_spline(time, psi) * Chi_spline(time - ageAtBMT) * donor_eps_spline(time) + rho_D * (2 * y[10] + y[9]) + Beta * y[11]  - (kloss + alpha + delta_D) * y[9];
  // Thymic ki lo displaceable
  dydt[10] = theta_spline(time, psi) * Chi_spline(time - ageAtBMT) * (1 - donor_eps_spline(time)) + kloss * y[9] + Beta * y[12]  - (rho_D + alpha + delta_D) * y[10];
  // Peripheral ki hi displaceable
  dydt[11] = alpha * y[9] + rho_D * (2 * y[12] + y[11]) - (kloss + Beta + delta_D) * y[11];
  // Peripheral ki lo displaceable
  dydt[12] = alpha * y[10] + kloss * y[11] - (rho_D + Beta + delta_D) * y[12];
  return dydt;
}

// solving for total counts of thymic and peripheral naive Tregs at time of BMT assuming the youngest animal as the t0
// these counts form the initial conditions for other recipients N(0) = N_host(0)
real[] solve_init(real ageAtBMT,
  real[] init_cond,                          // initial conditions at BMT in the youngest mouse
  real[] parms){

    real ta = 40;                            // age at BMT for the youngest host
    real y_init[2, 12];
    real params_init[7];

    params_init[1:6] = parms[1:6];
    params_init[7] = ta;

    y_init[1] = init_cond;                 // init conditions at the earliest BMT (i.e. in younegst animal)
    y_init[1] = init_cond;                 // init conditions at the earliest BMT (i.e. in younegst animal)
    if (ageAtBMT==40) {
      y_init[2] = init_cond;
    } else {
      y_init[2] = to_array_1d(integrate_ode_rk45(shm_chi, init_cond, ta, rep_array(ageAtBMT, 1), params_init, {0.0}, {0}));
    }
    return y_init[2];
}

real[] solve_chi(real solve_time, real ageAtBMT, real[] init_cond, real[] parms){
     real y_solve[12];
    real params[7];

    real y0[12];
    real init_tb[12];                         // init conditions at the mean age of BMT for the group

    //solution for the initial conditions at the mean age of BMT for the group
    y0 = solve_init(ageAtBMT, init_cond, parms);

    // init conditions at the BMT
    init_tb[1] = y0[1] + y0[9];
    init_tb[2] = y0[2] + y0[10];
    init_tb[3] = y0[3] + y0[11];                               //at tbmt - all cells are host
    init_tb[4] = y0[4] + y0[12];
    init_tb[5] = y0[5];
    init_tb[6] = y0[6];
    init_tb[7] = y0[7];
    init_tb[8] = y0[8];
    init_tb[9] = 0.0;                               //at tbmt - donor population = 0
    init_tb[10] = 0.0;
    init_tb[11] = 0.0;
    init_tb[12] = 0.0;

    params[1:6] = parms[1:6];
    params[7] = ageAtBMT;                                           // age at BMT

    y_solve = to_array_1d(integrate_ode_rk45(shm_chi, init_tb, ageAtBMT, rep_array(solve_time, 1), params, {0.0}, {0}));

    return y_solve;
  }

 real[,] solve_ode_chi(real[] solve_time, real[] ageAtBMT, real[] init_cond, real[] parms){
    int numdim = size(solve_time);
    real y_solve[numdim, 12];
    for (i in 1:numdim) {
      y_solve[i] = solve_chi(solve_time[i], ageAtBMT[i], init_cond, parms);
    }
    return y_solve;
  }

  vector math_reduce(vector global_params, vector local_params, real[] x_r, int[] x_i){
    // data for each shard
    int n = size(x_i);
    real solve_time[n] = x_r[1:n];
    int ageAtBMT[n] = x_i[1:n];                          // time zero -- for chimeras age at BMT
    real tb_time[n];

    //params
    real y1_0 = global_params[7]; real y2_0 = global_params[8];  real y3_0 = global_params[9];
    real y4_0 = global_params[10]; real y5_0 = global_params[11]; real y6_0 = global_params[12];
    real y7_0 = global_params[13]; real y8_0 = global_params[14];

    real init_cond[12];

    // ODE solution -- predictions for the observed timecourse
    real chi_solve[n, 12];

    real counts_thy[n]; real counts_per[n]; real donor_counts_thy[n]; real donor_counts_per[n];
    real host_counts_thy[n]; real host_counts_per[n]; real donor_ki_thy[n]; real donor_ki_per[n];
    real host_ki_thy[n]; real host_ki_per[n];

    vector[8*n] y_mean_stacked;

    // ODE solution -- predictions for the observed timecourse
    init_cond[1] = y1_0; init_cond[2] = y2_0; init_cond[3] = y3_0; init_cond[4] = y4_0;
    init_cond[5] = y5_0; init_cond[6] = y6_0; init_cond[7] = y7_0; init_cond[8] = y8_0;
    init_cond[9] =  0; init_cond[10] = 0; init_cond[11] = 0; init_cond[12] = 0;

    for (i in 1:n){
      tb_time[i] = ageAtBMT[i]/1.0;
    }

    // each shard has a single datpoint so its unique ****
    // PDE solution for chimera dataset -- x_r = data time and x_i = time at BMT
      chi_solve = solve_ode_chi(solve_time, tb_time, init_cond, to_array_1d(global_params));

      for (i in 1:n){
        counts_thy[i] = chi_solve[i, 1] + chi_solve[i, 2] + chi_solve[i, 7] + chi_solve[i, 8] + chi_solve[i, 9] + chi_solve[i, 10];
        counts_per[i] = chi_solve[i, 3] + chi_solve[i, 4] + chi_solve[i, 5] + chi_solve[i, 6] + chi_solve[i, 11] + chi_solve[i, 12];
        donor_counts_thy[i] = chi_solve[i, 9] + chi_solve[i, 10];
        donor_counts_per[i] = chi_solve[i, 11] + chi_solve[i, 12];
        host_counts_thy[i] = chi_solve[i, 1] + chi_solve[i, 2] + chi_solve[i, 7] + chi_solve[i, 8];
        host_counts_per[i] = chi_solve[i, 3] + chi_solve[i, 4] + chi_solve[i, 5] + chi_solve[i, 6];
        donor_ki_thy[i] = (chi_solve[i, 9])/donor_counts_thy[i];
        donor_ki_per[i] = (chi_solve[i, 11])/donor_counts_thy[i];
        host_ki_thy[i] = (chi_solve[i, 1] + chi_solve[i, 7])/host_counts_thy[i];
        host_ki_per[i] = (chi_solve[i, 3] + chi_solve[i, 5])/host_counts_per[i];

        y_mean_stacked[8*i - 7] = counts_thy[i];
        y_mean_stacked[8*i - 6] = counts_per[i];
        y_mean_stacked[8*i - 5] = donor_counts_thy[i]/(counts_thy[i] * Chi_spline(solve_time[i] - tb_time[i]));
        y_mean_stacked[8*i - 4] = donor_counts_per[i]/(counts_per[i] * Chi_spline(solve_time[i] - tb_time[i]));
        y_mean_stacked[8*i - 3] = host_ki_thy[i];
        y_mean_stacked[8*i - 2] = host_ki_per[i];
        y_mean_stacked[8*i - 1] = donor_ki_thy[i];
        y_mean_stacked[8*i - 0] = donor_ki_per[i];
      }

      return y_mean_stacked;
    }

    // functions for transformation of fractions in (0,a), where a >=1
    real logit_inverse(real x){
       real ans;
         ans = exp(x)/(1+exp(x));
         return ans;
    }

    // functions for transformation of fractions in (0,a), where a >=1
    real[] asinsqrt_array(real[] x){
      int ndims = size(x);
      real answer[ndims];
      real a = 1.2;

      for (i in 1: ndims){
        answer[i] = asin(sqrt(x[i])/sqrt(a));
      }
      return answer;
    }

    real asinsqrt_real(real x){
      real a = 1.2;

      real answer = asin(sqrt(x)/sqrt(a));
      return answer;
    }

    real asinsqrt_inv(real x){
      real a = 1.2;

      real answer = a * (sin(x))^2;
      return answer;
    }
  }


data{
   int numObs1;
   int numObs2;
   int numObs3;
   int numObs4;
   int n_solve;
   int n_shards;
   int numPred;
   int solve_time[n_solve];
   int ageAtBMT[n_solve];
   int dpBMT[n_solve];
   int time_index_counts[numObs1];
   int time_index_chi[numObs2];
   int time_index_donorki[numObs3];
   int time_index_hostki[numObs4];
   real counts_per[numObs1];
   real counts_thy[numObs1];
   real Nfd_thy[numObs2];
   real Nfd_per[numObs2];
   real ki_donor_per[numObs3];
   real ki_donor_thy[numObs3];
   real ki_host_per[numObs4];
   real ki_host_thy[numObs4];
   real ts_pred1[numPred];
   real ts_pred2[numPred];
   real ts_pred3[numPred];
   real ts_pred4[numPred];
   real tb_pred1[numPred];
   real tb_pred2[numPred];
   real tb_pred3[numPred];
   real tb_pred4[numPred];
}

transformed data{
   int M = n_solve/n_shards;     // pershard numobs
   int x_i[n_shards, M];         // each shard gets a single data point
   real x_r[n_shards, M];        // each shard gets a single data point

   // empty set of per shard params
   vector[0] local_params[n_shards];  // shard specific params --  useful for hierarchical modelling

   // data split into shards
   for (s in 1: n_shards){
     int i = 1 + (s-1) * M;                       // start index for ith shard
     int j = s * M;                               // end index for ith shard
     x_i[s, 1:M] = ageAtBMT[i:j];                       // age at BMT split
     x_r[s, 1:M] = solve_time[i:j];                     // solve time split
   }
}

parameters {
  real<lower= 0, upper= 1> psi;
  real<lower= 0, upper= 1> rho_D;
  real<lower= 0, upper= 1> alpha;
  real<lower= 0, upper= 1> delta_D;
  real<lower= 0, upper= 1> rho_I;
  real<lower= 0, upper= 1> Beta;
  real<lower= 0> y1_0;
  real<lower= 0> y2_0;
  real<lower= 0> y3_0;
  real<lower= 0> y4_0;
  real<lower= 0> y5_0;
  real<lower= 0> y6_0;
  real<lower= 0> y7_0;
  real<lower= 0> y8_0;


  real<lower=0> sigma_counts_per;
  real<lower=0> sigma_counts_thy;
  real<lower=0> sigma_Nfd_per;
  real<lower=0> sigma_Nfd_thy;
  real<lower=0> sigma_donor_ki_per;
  real<lower=0> sigma_donor_ki_thy;
  real<lower=0> sigma_host_ki_per;
  real<lower=0> sigma_host_ki_thy;
}

transformed parameters{
  vector[14] global_params;
  vector[n_solve] counts_thy_solve;               // ODE predictions for naive Treg counts in thymus
  vector[n_solve] counts_per_solve;               // ODE predictions for naive Treg counts in Periphery
  vector[n_solve] Nfd_thy_solve;                  // ODE predictions for naive Treg Nfd in thymus
  vector[n_solve] Nfd_per_solve;                  // ODE predictions for naive Treg Nfd in periphery
  vector[n_solve] ki_host_thy_solve;              // ODE predictions for naive Treg host ki67 proportions in thymus
  vector[n_solve] ki_host_per_solve;              // ODE predictions for naive Treg host ki67 proportions in periphery
  vector[n_solve] ki_donor_thy_solve;             // ODE predictions for naive Treg donor ki67 proportions in thymus
  vector[n_solve] ki_donor_per_solve;             // ODE predictions for naive Treg donor ki67 proportions in periphery
  vector[8*n_solve] y_mean_stacked;               // compliled output across all nodes

  vector[numObs1] counts_thy_mean;               // ODE predictions for naive Treg counts in thymus
  vector[numObs1] counts_per_mean;               // ODE predictions for naive Treg counts in Periphery
  vector[numObs2] Nfd_thy_mean;                  // ODE predictions for naive Treg Nfd in thymus
  vector[numObs2] Nfd_per_mean;                  // ODE predictions for naive Treg Nfd in periphery
  vector[numObs4] ki_host_thy_mean;              // ODE predictions for naive Treg host ki67 proportions in thymus
  vector[numObs4] ki_host_per_mean;              // ODE predictions for naive Treg host ki67 proportions in periphery
  vector[numObs3] ki_donor_thy_mean;             // ODE predictions for naive Treg donor ki67 proportions in thymus
  vector[numObs3] ki_donor_per_mean;             // ODE predictions for naive Treg donor ki67 proportions in periphery

  global_params[1] = psi;
  global_params[2] = rho_D;
  global_params[3] = alpha;
  global_params[4] = delta_D;
  global_params[5] = rho_I;
  global_params[6] = Beta;
  global_params[7] = exp(y1_0);
  global_params[8] = exp(y2_0);
  global_params[9] = exp(y3_0);
  global_params[10] = exp(y4_0);
  global_params[11] = exp(y5_0);
  global_params[12] = exp(y6_0);
  global_params[13] = exp(y7_0);
  global_params[14] = exp(y8_0);

  // combining the output from all the shards
  y_mean_stacked = map_rect(math_reduce, global_params, local_params, x_r, x_i);

  for (i in 1:n_solve){
    counts_thy_solve[i]   = y_mean_stacked[8*i - 7];
    counts_per_solve[i]   = y_mean_stacked[8*i - 6];
    Nfd_thy_solve[i]      = y_mean_stacked[8*i - 5];
    Nfd_per_solve[i]      = y_mean_stacked[8*i - 4];
    ki_host_thy_solve[i]  = y_mean_stacked[8*i - 3];
    ki_host_per_solve[i]  = y_mean_stacked[8*i - 2];
    ki_donor_thy_solve[i] = y_mean_stacked[8*i - 1];
    ki_donor_per_solve[i] = y_mean_stacked[8*i];
  }

  for (i in 1:numObs1){
    counts_thy_mean[i]   = counts_thy_solve[time_index_counts[i]];
    counts_per_mean[i]   = counts_per_solve[time_index_counts[i]];
  }
  for (i in 1:numObs2){
    Nfd_thy_mean[i]   = Nfd_thy_solve[time_index_chi[i]];
    Nfd_per_mean[i]   = Nfd_per_solve[time_index_chi[i]];
  }
  for (i in 1:numObs3){
    ki_donor_thy_mean[i]   = ki_donor_thy_solve[time_index_donorki[i]];
    ki_donor_per_mean[i]   = ki_donor_per_solve[time_index_donorki[i]];
  }
  for (i in 1:numObs4){
    ki_host_thy_mean[i]   = ki_host_thy_solve[time_index_hostki[i]];
    ki_host_per_mean[i]   = ki_host_per_solve[time_index_hostki[i]];
  }
}

model{
  psi ~ normal(0.3, 0.2);
  rho_D ~ normal(0.005, 0.25);
  alpha ~ normal(0.1, 0.25);
  delta_D ~ normal(0.01, 0.25);
  rho_I ~ normal(0.01, 0.25);
  Beta ~ normal(0.01, 0.25);
  y1_0 ~ normal(9, 2.5);
  y2_0 ~ normal(11, 2.5);
  y3_0 ~ normal(9, 2.5);
  y4_0 ~ normal(11, 2.5);
  y5_0 ~ normal(12, 2.5);
  y6_0 ~ normal(13, 2.5);
  y7_0 ~ normal(11, 2.5);
  y8_0 ~ normal(12, 2.5);

  sigma_counts_per ~ normal(0.5, 0.1);
  sigma_counts_thy ~ normal(0.7, 0.1);
  sigma_Nfd_per ~ normal(0.4, 0.05);
  sigma_Nfd_thy ~ normal(0.7, 0.05);
  sigma_donor_ki_per ~ normal(0.3, 0.5);
  sigma_donor_ki_thy ~ normal(0.3, 0.5);
  sigma_host_ki_per ~ normal(0.2, 2);
  sigma_host_ki_thy ~ normal(0.3, 2);

  log(counts_thy) ~ normal(log(counts_thy_mean), sigma_counts_thy);
  log(counts_per) ~ normal(log(counts_per_mean), sigma_counts_per);
  asinsqrt_array(Nfd_thy) ~ normal(asinsqrt_array(to_array_1d(Nfd_thy_mean)), sigma_Nfd_thy);
  asinsqrt_array(Nfd_per) ~ normal(asinsqrt_array(to_array_1d(Nfd_per_mean)), sigma_Nfd_per);
  asinsqrt_array(ki_donor_thy) ~ normal(asinsqrt_array(to_array_1d(ki_donor_thy_mean)), sigma_donor_ki_thy);
  asinsqrt_array(ki_donor_per) ~ normal(asinsqrt_array(to_array_1d(ki_donor_per_mean)), sigma_donor_ki_per);
  asinsqrt_array(ki_host_thy) ~ normal(asinsqrt_array(to_array_1d(ki_host_thy_mean)), sigma_host_ki_thy);
  asinsqrt_array(ki_host_per) ~ normal(asinsqrt_array(to_array_1d(ki_host_per_mean)), sigma_host_ki_per);
}

generated quantities{
  real y_chi_pred1[numPred, 12];
  real y_chi_pred2[numPred, 12];
  real y_chi_pred3[numPred, 12];
  real y_chi_pred4[numPred, 12];

  real counts_thy_mean_pred1[numPred];  real counts_per_mean_pred1[numPred];
  real counts_thy_mean_pred2[numPred];  real counts_per_mean_pred2[numPred];
  real counts_thy_mean_pred3[numPred];  real counts_per_mean_pred3[numPred];
  real counts_thy_mean_pred4[numPred];  real counts_per_mean_pred4[numPred];

  real counts_thy_pred1[numPred];  real counts_per_pred1[numPred];
  real counts_thy_pred2[numPred];  real counts_per_pred2[numPred];
  real counts_thy_pred3[numPred];  real counts_per_pred3[numPred];
  real counts_thy_pred4[numPred];  real counts_per_pred4[numPred];

  real Nfd_thy_mean_pred1[numPred];  real Nfd_per_mean_pred1[numPred];
  real Nfd_thy_mean_pred2[numPred];  real Nfd_per_mean_pred2[numPred];
  real Nfd_thy_mean_pred3[numPred];  real Nfd_per_mean_pred3[numPred];
  real Nfd_thy_mean_pred4[numPred];  real Nfd_per_mean_pred4[numPred];

  real Nfd_thy_pred1[numPred];  real Nfd_per_pred1[numPred];
  real Nfd_thy_pred2[numPred];  real Nfd_per_pred2[numPred];
  real Nfd_thy_pred3[numPred];  real Nfd_per_pred3[numPred];
  real Nfd_thy_pred4[numPred];  real Nfd_per_pred4[numPred];

  real ki_donor_per_mean_pred1[numPred]; real ki_donor_thy_mean_pred1[numPred];
  real ki_donor_per_mean_pred2[numPred]; real ki_donor_thy_mean_pred2[numPred];
  real ki_donor_per_mean_pred3[numPred]; real ki_donor_thy_mean_pred3[numPred];
  real ki_donor_per_mean_pred4[numPred]; real ki_donor_thy_mean_pred4[numPred];

  real ki_donor_per_pred1[numPred]; real ki_donor_thy_pred1[numPred];
  real ki_donor_per_pred2[numPred]; real ki_donor_thy_pred2[numPred];
  real ki_donor_per_pred3[numPred]; real ki_donor_thy_pred3[numPred];
  real ki_donor_per_pred4[numPred]; real ki_donor_thy_pred4[numPred];

  real ki_host_per_mean_pred1[numPred]; real ki_host_thy_mean_pred1[numPred];
  real ki_host_per_mean_pred2[numPred]; real ki_host_thy_mean_pred2[numPred];
  real ki_host_per_mean_pred3[numPred]; real ki_host_thy_mean_pred3[numPred];
  real ki_host_per_mean_pred4[numPred]; real ki_host_thy_mean_pred4[numPred];

  real ki_host_per_pred1[numPred]; real ki_host_thy_pred1[numPred];
  real ki_host_per_pred2[numPred]; real ki_host_thy_pred2[numPred];
  real ki_host_per_pred3[numPred]; real ki_host_thy_pred3[numPred];
  real ki_host_per_pred4[numPred]; real ki_host_thy_pred4[numPred];

  // log likelihoods
  vector[numObs1] log_lik_counts_thy;
  vector[numObs1] log_lik_counts_per;
  vector[numObs2] log_lik_Nfd_thy;
  vector[numObs2] log_lik_Nfd_per;
  vector[numObs3] log_lik_ki_donor_thy;
  vector[numObs3] log_lik_ki_donor_per;
  vector[numObs4] log_lik_ki_host_thy;
  vector[numObs4] log_lik_ki_host_per;

  // initial conditions
  real init_cond[12];
  init_cond[1] = exp(y1_0); init_cond[2] = exp(y2_0); init_cond[3] = exp(y3_0); init_cond[4] = exp(y4_0);
  init_cond[5] = exp(y5_0); init_cond[6] = exp(y6_0); init_cond[7] = exp(y7_0); init_cond[8] = exp(y8_0);
  init_cond[9] =  0; init_cond[10] = 0; init_cond[11] = 0; init_cond[12] = 0;

  // predictions for the whole timecourse
  y_chi_pred1 = solve_ode_chi(ts_pred1, tb_pred1, init_cond, to_array_1d(global_params));
  y_chi_pred2 = solve_ode_chi(ts_pred2, tb_pred2, init_cond, to_array_1d(global_params));
  y_chi_pred3 = solve_ode_chi(ts_pred3, tb_pred3, init_cond, to_array_1d(global_params));
  y_chi_pred4 = solve_ode_chi(ts_pred4, tb_pred4, init_cond, to_array_1d(global_params));

  for (i in 1:numPred){
    counts_thy_mean_pred1[i] = y_chi_pred1[i, 1] + y_chi_pred1[i, 2] + y_chi_pred1[i, 7] + y_chi_pred1[i, 8] + y_chi_pred1[i, 9] + y_chi_pred1[i, 10];
    counts_thy_mean_pred2[i] = y_chi_pred2[i, 1] + y_chi_pred2[i, 2] + y_chi_pred2[i, 7] + y_chi_pred2[i, 8] + y_chi_pred2[i, 9] + y_chi_pred2[i, 10];
    counts_thy_mean_pred3[i] = y_chi_pred3[i, 1] + y_chi_pred3[i, 2] + y_chi_pred3[i, 7] + y_chi_pred3[i, 8] + y_chi_pred3[i, 9] + y_chi_pred3[i, 10];
    counts_thy_mean_pred4[i] = y_chi_pred4[i, 1] + y_chi_pred4[i, 2] + y_chi_pred4[i, 7] + y_chi_pred4[i, 8] + y_chi_pred4[i, 9] + y_chi_pred4[i, 10];

    counts_thy_pred1[i] = exp(normal_rng(log(counts_thy_mean_pred1[i]), sigma_counts_thy));
    counts_thy_pred2[i] = exp(normal_rng(log(counts_thy_mean_pred2[i]), sigma_counts_thy));
    counts_thy_pred3[i] = exp(normal_rng(log(counts_thy_mean_pred3[i]), sigma_counts_thy));
    counts_thy_pred4[i] = exp(normal_rng(log(counts_thy_mean_pred4[i]), sigma_counts_thy));

    counts_per_mean_pred1[i] = y_chi_pred1[i, 3] + y_chi_pred1[i, 4] + y_chi_pred1[i, 5] + y_chi_pred1[i, 6] + y_chi_pred1[i, 11] + y_chi_pred1[i, 12];
    counts_per_mean_pred2[i] = y_chi_pred2[i, 3] + y_chi_pred2[i, 4] + y_chi_pred2[i, 5] + y_chi_pred2[i, 6] + y_chi_pred2[i, 11] + y_chi_pred2[i, 12];
    counts_per_mean_pred3[i] = y_chi_pred3[i, 3] + y_chi_pred3[i, 4] + y_chi_pred3[i, 5] + y_chi_pred3[i, 6] + y_chi_pred3[i, 11] + y_chi_pred3[i, 12];
    counts_per_mean_pred4[i] = y_chi_pred4[i, 3] + y_chi_pred4[i, 4] + y_chi_pred4[i, 5] + y_chi_pred4[i, 6] + y_chi_pred4[i, 11] + y_chi_pred4[i, 12];

    counts_per_pred1[i] = exp(normal_rng(log(counts_per_mean_pred1[i]), sigma_counts_per));
    counts_per_pred2[i] = exp(normal_rng(log(counts_per_mean_pred2[i]), sigma_counts_per));
    counts_per_pred3[i] = exp(normal_rng(log(counts_per_mean_pred3[i]), sigma_counts_per));
    counts_per_pred4[i] = exp(normal_rng(log(counts_per_mean_pred4[i]), sigma_counts_per));

    Nfd_thy_mean_pred1[i] = (y_chi_pred1[i, 9] + y_chi_pred1[i, 10])/(counts_thy_mean_pred1[i] * Chi_spline(ts_pred1[i] - tb_pred1[i]));
    Nfd_thy_mean_pred2[i] = (y_chi_pred2[i, 9] + y_chi_pred2[i, 10])/(counts_thy_mean_pred2[i] * Chi_spline(ts_pred2[i] - tb_pred2[i]));
    Nfd_thy_mean_pred3[i] = (y_chi_pred3[i, 9] + y_chi_pred3[i, 10])/(counts_thy_mean_pred3[i] * Chi_spline(ts_pred3[i] - tb_pred3[i]));
    Nfd_thy_mean_pred4[i] = (y_chi_pred4[i, 9] + y_chi_pred4[i, 10])/(counts_thy_mean_pred4[i] * Chi_spline(ts_pred4[i] - tb_pred4[i]));

    Nfd_thy_pred1[i] = asinsqrt_inv(normal_rng(asinsqrt_real(Nfd_thy_mean_pred1[i]), sigma_Nfd_thy));
    Nfd_thy_pred2[i] = asinsqrt_inv(normal_rng(asinsqrt_real(Nfd_thy_mean_pred2[i]), sigma_Nfd_thy));
    Nfd_thy_pred3[i] = asinsqrt_inv(normal_rng(asinsqrt_real(Nfd_thy_mean_pred3[i]), sigma_Nfd_thy));
    Nfd_thy_pred4[i] = asinsqrt_inv(normal_rng(asinsqrt_real(Nfd_thy_mean_pred4[i]), sigma_Nfd_thy));

    Nfd_per_mean_pred1[i] = (y_chi_pred1[i, 11] + y_chi_pred1[i, 12])/(counts_per_mean_pred1[i] * Chi_spline(ts_pred1[i] - tb_pred1[i]));
    Nfd_per_mean_pred2[i] = (y_chi_pred2[i, 11] + y_chi_pred2[i, 12])/(counts_per_mean_pred2[i] * Chi_spline(ts_pred2[i] - tb_pred2[i]));
    Nfd_per_mean_pred3[i] = (y_chi_pred3[i, 11] + y_chi_pred3[i, 12])/(counts_per_mean_pred3[i] * Chi_spline(ts_pred3[i] - tb_pred3[i]));
    Nfd_per_mean_pred4[i] = (y_chi_pred4[i, 11] + y_chi_pred4[i, 12])/(counts_per_mean_pred4[i] * Chi_spline(ts_pred4[i] - tb_pred4[i]));

    Nfd_per_pred1[i] = asinsqrt_inv(normal_rng(asinsqrt_real(Nfd_per_mean_pred1[i]), sigma_Nfd_per));
    Nfd_per_pred2[i] = asinsqrt_inv(normal_rng(asinsqrt_real(Nfd_per_mean_pred2[i]), sigma_Nfd_per));
    Nfd_per_pred3[i] = asinsqrt_inv(normal_rng(asinsqrt_real(Nfd_per_mean_pred3[i]), sigma_Nfd_per));
    Nfd_per_pred4[i] = asinsqrt_inv(normal_rng(asinsqrt_real(Nfd_per_mean_pred4[i]), sigma_Nfd_per));

    ki_donor_thy_mean_pred1[i] = (y_chi_pred1[i, 9])/(y_chi_pred1[i, 9] + y_chi_pred1[i, 10]);
    ki_donor_thy_mean_pred2[i] = (y_chi_pred2[i, 9])/(y_chi_pred2[i, 9] + y_chi_pred2[i, 10]);
    ki_donor_thy_mean_pred3[i] = (y_chi_pred3[i, 9])/(y_chi_pred3[i, 9] + y_chi_pred3[i, 10]);
    ki_donor_thy_mean_pred4[i] = (y_chi_pred4[i, 9])/(y_chi_pred4[i, 9] + y_chi_pred4[i, 10]);

    ki_donor_thy_pred1[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_donor_thy_mean_pred1[i]), sigma_donor_ki_thy));
    ki_donor_thy_pred2[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_donor_thy_mean_pred2[i]), sigma_donor_ki_thy));
    ki_donor_thy_pred3[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_donor_thy_mean_pred3[i]), sigma_donor_ki_thy));
    ki_donor_thy_pred4[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_donor_thy_mean_pred4[i]), sigma_donor_ki_thy));

    ki_donor_per_mean_pred1[i] = (y_chi_pred1[i, 11])/(y_chi_pred1[i, 11] + y_chi_pred1[i, 12]);
    ki_donor_per_mean_pred2[i] = (y_chi_pred2[i, 11])/(y_chi_pred2[i, 11] + y_chi_pred2[i, 12]);
    ki_donor_per_mean_pred3[i] = (y_chi_pred3[i, 11])/(y_chi_pred3[i, 11] + y_chi_pred3[i, 12]);
    ki_donor_per_mean_pred4[i] = (y_chi_pred4[i, 11])/(y_chi_pred4[i, 11] + y_chi_pred4[i, 12]);

    ki_donor_per_pred1[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_donor_thy_mean_pred1[i]), sigma_donor_ki_per));
    ki_donor_per_pred2[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_donor_thy_mean_pred2[i]), sigma_donor_ki_per));
    ki_donor_per_pred3[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_donor_thy_mean_pred3[i]), sigma_donor_ki_per));
    ki_donor_per_pred4[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_donor_thy_mean_pred4[i]), sigma_donor_ki_per));

    ki_host_thy_mean_pred1[i] = (y_chi_pred1[i, 1] + y_chi_pred1[i, 7])/(y_chi_pred1[i, 1] + y_chi_pred1[i, 2] + y_chi_pred1[i, 7] + y_chi_pred1[i, 8]);
    ki_host_thy_mean_pred2[i] = (y_chi_pred2[i, 1] + y_chi_pred2[i, 7])/(y_chi_pred2[i, 1] + y_chi_pred2[i, 2] + y_chi_pred2[i, 7] + y_chi_pred2[i, 8]);
    ki_host_thy_mean_pred3[i] = (y_chi_pred3[i, 1] + y_chi_pred3[i, 7])/(y_chi_pred3[i, 1] + y_chi_pred3[i, 2] + y_chi_pred3[i, 7] + y_chi_pred3[i, 8]);
    ki_host_thy_mean_pred4[i] = (y_chi_pred4[i, 1] + y_chi_pred4[i, 7])/(y_chi_pred4[i, 1] + y_chi_pred4[i, 2] + y_chi_pred4[i, 7] + y_chi_pred4[i, 8]);

    ki_host_thy_pred1[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_host_thy_mean_pred1[i]), sigma_host_ki_thy));
    ki_host_thy_pred2[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_host_thy_mean_pred2[i]), sigma_host_ki_thy));
    ki_host_thy_pred3[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_host_thy_mean_pred3[i]), sigma_host_ki_thy));
    ki_host_thy_pred4[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_host_thy_mean_pred4[i]), sigma_host_ki_thy));

    ki_host_per_mean_pred1[i] = (y_chi_pred1[i, 3] + y_chi_pred1[i, 5])/(y_chi_pred1[i, 3] + y_chi_pred1[i, 4] + y_chi_pred1[i, 5] + y_chi_pred1[i, 6]);
    ki_host_per_mean_pred2[i] = (y_chi_pred2[i, 3] + y_chi_pred2[i, 5])/(y_chi_pred2[i, 3] + y_chi_pred2[i, 4] + y_chi_pred2[i, 5] + y_chi_pred2[i, 6]);
    ki_host_per_mean_pred3[i] = (y_chi_pred3[i, 3] + y_chi_pred3[i, 5])/(y_chi_pred3[i, 3] + y_chi_pred3[i, 4] + y_chi_pred3[i, 5] + y_chi_pred3[i, 6]);
    ki_host_per_mean_pred4[i] = (y_chi_pred4[i, 3] + y_chi_pred4[i, 5])/(y_chi_pred4[i, 3] + y_chi_pred4[i, 4] + y_chi_pred4[i, 5] + y_chi_pred4[i, 6]);

    ki_host_per_pred1[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_host_per_mean_pred1[i]), sigma_host_ki_per));
    ki_host_per_pred2[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_host_per_mean_pred2[i]), sigma_host_ki_per));
    ki_host_per_pred3[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_host_per_mean_pred3[i]), sigma_host_ki_per));
    ki_host_per_pred4[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_host_per_mean_pred4[i]), sigma_host_ki_per));

  }

  // calculating log likelihoods
  for (i in 1:numObs1) {
    log_lik_counts_thy[i] = normal_lpdf(log(counts_thy[i]) | log(counts_thy_mean[i]), sigma_counts_thy);
    log_lik_counts_per[i] = normal_lpdf(log(counts_per[i]) | log(counts_per_mean[i]), sigma_counts_per);
  }
  for (i in 1:numObs2) {
    log_lik_Nfd_thy[i] = normal_lpdf(asinsqrt_real(Nfd_thy[i]) | asinsqrt_real(Nfd_thy_mean[i]), sigma_Nfd_thy);
    log_lik_Nfd_per[i] = normal_lpdf(asinsqrt_real(Nfd_per[i]) | asinsqrt_real(Nfd_per_mean[i]), sigma_Nfd_per);
  }
  for (i in 1:numObs3) {
    log_lik_ki_donor_thy[i] = normal_lpdf(asinsqrt_real(ki_donor_thy[i]) | asinsqrt_real(ki_donor_thy_mean[i]), sigma_donor_ki_thy);
    log_lik_ki_donor_per[i] = normal_lpdf(asinsqrt_real(ki_donor_per[i]) | asinsqrt_real(ki_donor_per_mean[i]), sigma_donor_ki_per);
  }
  for (i in 1:numObs4) {
    log_lik_ki_host_thy[i] = normal_lpdf(asinsqrt_real(ki_host_thy[i]) | asinsqrt_real(ki_host_thy_mean[i]), sigma_host_ki_thy);
    log_lik_ki_host_per[i] = normal_lpdf(asinsqrt_real(ki_host_per[i]) | asinsqrt_real(ki_host_per_mean[i]), sigma_host_ki_per);
  }
}
