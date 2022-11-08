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
   real value;  real chiEst = 0.79495993;   real qEst = 0.09859144;
   if (time - 14 < 0){              // t0 = 14 -- assumption: no donor cells seen in FoxP3neg SP4 compartment for 2 weeks
     value = 0;
   } else {
     value = chiEst * (1 - exp(-qEst * (time - 14)));
   }
   return value;
 }

 // spline3 --
 // proportions of ki67hi cells in the donor-derived FoxP3 negative SP4 T cells -- varies with time
 real donor_eps_spline(real time){
   real t0 = 49.0;     // mean mouse age at BMT for the first ageBMT bin
  //parameters estimated from spline fit to the timecourse of ki67 fraction in the donor-derived FoxP3 negative SP4 T cells
  real eps_0 = 0.37848471; real eps_f = 0.06799028;
  return exp(- eps_f * (time - t0)) + eps_0;
}

real[] shm_chi(real time, real[] y, real[] parms, real[] rdata,  int[] idata) {
  real psi = parms[1];
  real rho_0 = parms[2];
  real alpha = parms[3];
  real delta_0 = parms[4];
  real mu = parms[5];
  real rho = parms[6];
  real beta = parms[7];
  real delta= parms[8];

  real dydt[16];
  real kloss  = 1/3.5;            //rate of loss of ki67
  real eps_host = 0.326611;      // Mean Ki67 hi fraction in host-BM-derived FoxP3 negative Sp4 T cells

  // age of BMT in each recipient
  real ageAtBMT = parms[9];

  // model that assumes that tranistionals divide and die at different rates than mature naive T cells
  // Host naive Tregs
  // Thymic ki  hi tranistionals
  dydt[1] = theta_spline(time, psi) * (1- Chi_spline(time - ageAtBMT)) * eps_host + rho_0 * (2 * y[2] + y[1]) - (kloss + alpha + delta_0) * y[1];
  // Thymic ki lo tranistionals
  dydt[2] = theta_spline(time, psi) * (1- Chi_spline(time - ageAtBMT)) * (1 - eps_host) + kloss * y[1] - (rho_0 + alpha + delta_0) * y[2];
  // Peripheral ki hi tranistionals
  dydt[3] = alpha * y[1] + rho_0 * (2 * y[4] + y[3]) - (kloss + mu + delta_0) * y[3];
  // Peripheral ki lo tranistionals
  dydt[4] = alpha * y[2] + kloss * y[3] - (rho_0 + mu + delta_0) * y[4];
  // Peripheral ki hi mature
  dydt[5] = mu * y[3] + alpha * y[7] + rho * (2 * y[6] + y[5]) - (kloss + beta + delta) * y[5];
  // Peripheral ki lo mature
  dydt[6] = mu * y[4] + alpha * y[8] + kloss * y[5] - (rho + beta + delta) * y[6];
  // Thymic ki hi mature
  dydt[7] = beta * y[5] + rho * (2 * y[8] + y[7]) - (kloss + alpha + delta) * y[7];
  // Thymic ki lo mature
  dydt[8] = beta * y[6] + kloss * y[7] - (rho + alpha + delta) * y[8];

  // Donor naive Tregs
  // Thymic ki  hi tranistionals
  dydt[9] = theta_spline(time, psi) * Chi_spline(time - ageAtBMT) * donor_eps_spline(time) + rho_0 * (2 * y[10] + y[9]) - (kloss + alpha + delta_0) * y[9];
  // Thymic ki lo tranistionals
  dydt[10] = theta_spline(time, psi) * Chi_spline(time - ageAtBMT) * (1 - donor_eps_spline(time)) + kloss * y[9] - (rho_0 + alpha + delta_0) * y[10];
  // Peripheral ki hi tranistionals
  dydt[11] = alpha * y[9] + rho_0 * (2 * y[12] + y[11]) - (kloss + mu + delta_0) * y[11];
  // Peripheral ki lo tranistionals
  dydt[12] = alpha * y[10] + kloss * y[11] - (rho_0 + mu + delta_0) * y[12];
  // Peripheral ki hi mature
  dydt[13] = mu * y[11] + alpha * y[15] + rho * (2 * y[14] + y[13]) - (kloss + beta + delta) * y[13];
  // Peripheral ki lo mature
  dydt[14] = mu * y[12] + alpha * y[16] + kloss * y[13] - (rho + beta + delta) * y[14];
  // Thymic ki hi mature
  dydt[15] = beta * y[13] + rho * (2 * y[16] + y[15]) - (kloss + alpha + delta) * y[15];
  // Thymic ki lo mature
  dydt[16] = beta * y[14] + kloss * y[15] - (rho + alpha + delta) * y[16];

  return dydt;
}

 real[] solve_chi(real solve_time, real ageAtBMT, real[] init_cond, real[] parms){
    real y_solve[16];
    real params[9];
    params[1:8] = parms[1:8];
    params[9] = ageAtBMT;                      // age at BMT
    y_solve = to_array_1d(integrate_ode_rk45(shm_chi, init_cond, ageAtBMT, rep_array(solve_time, 1), params, {0.0}, {0}));
    return y_solve;
  }

 real[,] solve_ode_chi(real[] solve_time, real[] ageAtBMT, real[] init_cond, real[] parms){
    int numdim = size(solve_time);
    real y_solve[numdim, 16];
    for (i in 1:numdim) {
      y_solve[i] = solve_chi(solve_time[i], ageAtBMT[i], init_cond, parms);
    }
    return y_solve;
  }

  vector math_reduce(vector global_params, vector local_params, real[] x_r, int[] x_i){
    // data for each shard
    int n = size(x_i); // n = 1
    real solve_time = x_r[1];
    int ageAtBMT = x_i[1];                          // time zero -- for chimeras age at BMT
    real tb_time = dat_t0/1.0;

    //params
    real N0 = global_params[7];
    real kappa0 = global_params[8];
    real init_cond[4];

    // ODE solution -- predictions for the observed timecourse
    real chi_solve[8];

    real chi_counts_mean;
    real host_counts_mean;
    real donor_counts_mean;
    real host_ki_mean;
    real donor_ki_mean;

    vector[4*n] y_mean_stacked;

    // ODE solution -- predictions for the observed timecourse
    init_cond[1] = kappa0 * N0;
    init_cond[2] = (1 - kappa0) * N0;
    init_cond[3] = 0.0;
    init_cond[4] = 0.0;

    // each shard has a single datpoint so its unique ****
    // PDE solution for chimera dataset -- x_r = data time and x_i = time at BMT
      chi_solve = solve_chi(dat_time, tb_time, init_cond, to_array_1d(global_params));
      chi_counts_mean = chi_solve[1] + chi_solve[2] + chi_solve[3] + chi_solve[4] + chi_solve[5] + chi_solve[6] + chi_solve[7] + chi_solve[8];
      donor_counts_mean = chi_solve[1] + chi_solve[2] + chi_solve[3] + chi_solve[4];
      host_counts_mean = chi_solve[5] + chi_solve[6] + chi_solve[7] + chi_solve[8];

      donor_ki_mean = (chi_solve[1] + chi_solve[3])/donor_counts_mean;
      host_ki_mean = (chi_solve[5] + chi_solve[7])/host_counts_mean;
      y_mean_stacked[1] = chi_counts_mean;
      y_mean_stacked[2] = donor_counts_mean/(chi_counts_mean * Chi_spline(dat_time - tb_time));
      y_mean_stacked[3] = host_ki_mean;
      y_mean_stacked[4] = donor_ki_mean;

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
   real Nfd_per[numObs2];
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
   int x_i[n_shards, 1];         // each shard gets a single data point
   real x_r[n_shards, 1];        // each shard gets a single data point
   // empty set of per shard params
   vector[0] local_params[n_shards];  // shard specific params --  useful for hierarchical modelling
   // data split into shards
   for (s in 1:n_shards){
    x_i[s, 1] = ageAtBMT[s];                       // age at BMT split
    x_r[s, 1] = solve_time[s];                     // time split
   }
}

parameters {
  real<lower= 0, upper= 1> psi;
  real<lower= 0, upper= 1> rho_0;
  real<lower= 0, upper= 1> alpha;
  real<lower= 0, upper= 1> delta_0;
  real<lower= 0, upper= 1> mu;
  real<lower= 0, upper= 1> rho;
  real<lower= 0, upper= 1> beta;
  real<lower= 0, upper= 1> delta;
  real<lower= 11, upper= 16> N0_per;
  real<lower= 0, upper= 1> k0_per;
  real<lower= 9, upper= 13> N0_thy;
  real<lower= 0, upper= 1> k0_thy;

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
  vector[8] global_params;
  vector[numObs1] counts_thy_mean;               // ODE predictions for naive Treg counts in thymus
  vector[numObs1] counts_per_mean;               // ODE predictions for naive Treg counts in Periphery
  vector[numObs2] Nfd_thy_mean;                  // ODE predictions for naive Treg Nfd in thymus
  vector[numObs2] Nfd_per_mean;                  // ODE predictions for naive Treg Nfd in periphery
  vector[numObs3] ki_host_thy_mean;              // ODE predictions for naive Treg host ki67 proportions in thymus
  vector[numObs3] ki_host_per_mean;              // ODE predictions for naive Treg host ki67 proportions in periphery
  vector[numObs4] ki_donor_thy_mean;             // ODE predictions for naive Treg donor ki67 proportions in thymus
  vector[numObs4] ki_donor_per_mean;             // ODE predictions for naive Treg donor ki67 proportions in periphery
  vector[2* (numObs1 + numObs2 + numObs3 + numObs4)] y_mean_stacked;               // compliled output across all nodes

  global_params[1] = psi;
  global_params[2] = rho_0;
  global_params[3] = alpha;
  global_params[4] = delta_0;
  global_params[5] = mu;
  global_params[6] = rho;
  global_params[7] = beta;
  global_params[8] = delta;

  // combining the output from all the shards
  y_mean_stacked = map_rect(math_reduce, global_params, local_params, x_r, x_i);

  for (i in 1:numChi){
    counts_thy_mean[i] = y_mean_stacked[4*i - 3];
    y4_mean[i] = y_mean_stacked[4*i - 2];
    y5_mean[i] = y_mean_stacked[4*i - 1];
    y6_mean[i] = y_mean_stacked[4*i];
  }
}

model{
  psi ~ normal(0.3, 0.2);
  rho_nai ~ normal(0.005, 0.25);
  rho_rte ~ normal(0.1, 0.25);
  delta_nai ~ normal(0.01, 0.25);
  delta_rte ~ normal(0.01, 0.25);
  mu ~ normal(0.1, 0.025);

  N0 ~ normal(9E5, 1E5);
  kappa0 ~ normal(0.8, 0.1);

  sigma_chi_counts ~ normal(0, 2);
  sigma_Nfd ~ normal(0, 2);
  sigma_donor_ki ~ normal(0, 2);
  sigma_host_ki ~ normal(0, 2);

  log(chi_counts) ~ normal(log(y3_mean), sigma_chi_counts);
  logit(N_donor_fraction) ~ normal(logit(to_array_1d(y4_mean)), sigma_Nfd);
  asinsqrt_array(host_ki) ~ normal(asinsqrt_array(to_array_1d(y5_mean)), sigma_host_ki);
  asinsqrt_array(donor_ki) ~ normal(asinsqrt_array(to_array_1d(y6_mean)), sigma_donor_ki);
}

generated quantities{
  real y_chi_pred1[numPred, 8];
  real y_chi_pred2[numPred, 8];
  real y_chi_pred3[numPred, 8];

  real y3_mean_pred1[numPred];  real y4_mean_pred1[numPred]; real y5_mean_pred1[numPred];  real y6_mean_pred1[numPred];
  real y3_mean_pred2[numPred];  real y4_mean_pred2[numPred]; real y5_mean_pred2[numPred];  real y6_mean_pred2[numPred];
  real y3_mean_pred3[numPred];  real y4_mean_pred3[numPred]; real y5_mean_pred3[numPred];  real y6_mean_pred3[numPred];

  real chicounts_pred1[numPred]; real Nfd_pred1[numPred];
  real chicounts_pred2[numPred]; real Nfd_pred2[numPred];
  real chicounts_pred3[numPred]; real Nfd_pred3[numPred];

  real donorki_pred1[numPred]; real hostki_pred1[numPred];
  real donorki_pred2[numPred]; real hostki_pred2[numPred];
  real donorki_pred3[numPred]; real hostki_pred3[numPred];

  // log likelihoods
  vector[numChi] log_lik_chi_counts;
  vector[numChi] log_lik_Nfd;
  vector[numChi] log_lik_donor_ki;
  vector[numChi] log_lik_host_ki;

  // initial conditions
  real init_cond[4];
  init_cond[1] = kappa0 * N0;
  init_cond[2] = (1 - kappa0) * N0;
  init_cond[3] = 0.0;
  init_cond[4] = 0.0;

  // predictions for the whole timecourse
  y_chi_pred1 = solve_ode_chi(ts_pred_chi1, tb_pred1, init_cond, to_array_1d(global_params));
  y_chi_pred2 = solve_ode_chi(ts_pred_chi2, tb_pred2, init_cond, to_array_1d(global_params));
  y_chi_pred3 = solve_ode_chi(ts_pred_chi3, tb_pred3, init_cond, to_array_1d(global_params));

  for (i in 1:numPred){
    y3_mean_pred1[i] = y_chi_pred1[i, 1] + y_chi_pred1[i, 2] + y_chi_pred1[i, 3] + y_chi_pred1[i, 4] + y_chi_pred1[i, 5] + y_chi_pred1[i, 6] + y_chi_pred1[i, 7] + y_chi_pred1[i, 8];
    y3_mean_pred2[i] = y_chi_pred2[i, 1] + y_chi_pred2[i, 2] + y_chi_pred2[i, 3] + y_chi_pred2[i, 4] + y_chi_pred2[i, 5] + y_chi_pred2[i, 6] + y_chi_pred2[i, 7] + y_chi_pred2[i, 8];
    y3_mean_pred3[i] = y_chi_pred3[i, 1] + y_chi_pred3[i, 2] + y_chi_pred3[i, 3] + y_chi_pred3[i, 4] + y_chi_pred3[i, 5] + y_chi_pred3[i, 6] + y_chi_pred3[i, 7] + y_chi_pred3[i, 8];

    y4_mean_pred1[i] = (y_chi_pred1[i, 1] + y_chi_pred1[i, 2] + y_chi_pred1[i, 3] + y_chi_pred1[i, 4])/(y3_mean_pred1[i] * Chi_spline(ts_pred_chi1[i] - 54));
    y4_mean_pred2[i] = (y_chi_pred2[i, 1] + y_chi_pred2[i, 2] + y_chi_pred2[i, 3] + y_chi_pred2[i, 4])/(y3_mean_pred2[i] * Chi_spline(ts_pred_chi2[i] - 71));
    y4_mean_pred3[i] = (y_chi_pred3[i, 1] + y_chi_pred3[i, 2] + y_chi_pred3[i, 3] + y_chi_pred3[i, 4])/(y3_mean_pred3[i] * Chi_spline(ts_pred_chi3[i] - 97));

    y5_mean_pred1[i] = (y_chi_pred1[i, 5] + y_chi_pred1[i, 7])/(y_chi_pred1[i, 5] + y_chi_pred1[i, 6] + y_chi_pred1[i, 7] + y_chi_pred1[i, 8]);
    y5_mean_pred2[i] = (y_chi_pred2[i, 5] + y_chi_pred2[i, 7])/(y_chi_pred2[i, 5] + y_chi_pred2[i, 6] + y_chi_pred2[i, 7] + y_chi_pred2[i, 8]);
    y5_mean_pred3[i] = (y_chi_pred3[i, 5] + y_chi_pred3[i, 7])/(y_chi_pred3[i, 5] + y_chi_pred3[i, 6] + y_chi_pred3[i, 7] + y_chi_pred3[i, 8]);

    y6_mean_pred1[i] = (y_chi_pred1[i, 1] + y_chi_pred1[i, 3])/(y_chi_pred1[i, 1] + y_chi_pred1[i, 2] + y_chi_pred1[i, 3] + y_chi_pred1[i, 4]);
    y6_mean_pred2[i] = (y_chi_pred2[i, 1] + y_chi_pred2[i, 3])/(y_chi_pred2[i, 1] + y_chi_pred2[i, 2] + y_chi_pred2[i, 3] + y_chi_pred2[i, 4]);
    y6_mean_pred3[i] = (y_chi_pred3[i, 1] + y_chi_pred3[i, 3])/(y_chi_pred3[i, 1] + y_chi_pred3[i, 2] + y_chi_pred3[i, 3] + y_chi_pred3[i, 4]);


    chicounts_pred1[i] = exp(normal_rng(log(y3_mean_pred1[i]), sigma_chi_counts));
    chicounts_pred2[i] = exp(normal_rng(log(y3_mean_pred2[i]), sigma_chi_counts));
    chicounts_pred3[i] = exp(normal_rng(log(y3_mean_pred3[i]), sigma_chi_counts));

    Nfd_pred1[i] = logit_inverse(normal_rng(logit(y4_mean_pred1[i]), sigma_Nfd));
    Nfd_pred2[i] = logit_inverse(normal_rng(logit(y4_mean_pred2[i]), sigma_Nfd));
    Nfd_pred3[i] = logit_inverse(normal_rng(logit(y4_mean_pred3[i]), sigma_Nfd));

    hostki_pred1[i] = asinsqrt_inv(normal_rng(asinsqrt_real(y5_mean_pred1[i]), sigma_host_ki));
    hostki_pred2[i] = asinsqrt_inv(normal_rng(asinsqrt_real(y5_mean_pred2[i]), sigma_host_ki));
    hostki_pred3[i] = asinsqrt_inv(normal_rng(asinsqrt_real(y5_mean_pred3[i]), sigma_host_ki));

    donorki_pred1[i] = asinsqrt_inv(normal_rng(asinsqrt_real(y6_mean_pred1[i]), sigma_donor_ki));
    donorki_pred2[i] = asinsqrt_inv(normal_rng(asinsqrt_real(y6_mean_pred2[i]), sigma_donor_ki));
    donorki_pred3[i] = asinsqrt_inv(normal_rng(asinsqrt_real(y6_mean_pred3[i]), sigma_donor_ki));
  }

  // calculating log likelihoods
  for (i in 1:numChi) {
    log_lik_chi_counts[i] = normal_lpdf(log(chi_counts[i]) | log(y3_mean[i]), sigma_chi_counts);
    log_lik_Nfd[i]        = normal_lpdf(logit(N_donor_fraction[i]) | logit(y4_mean[i]), sigma_Nfd);
    log_lik_host_ki[i]    = normal_lpdf(asinsqrt_real(host_ki[i]) | asinsqrt_real(y5_mean[i]), sigma_host_ki);
    log_lik_donor_ki[i]   = normal_lpdf(asinsqrt_real(donor_ki[i]) | asinsqrt_real(y6_mean[i]), sigma_donor_ki);
  }
}
