functions{
 //spline 1
 // Timecourse of thymic CD4 SP population -- changes with time
 real sp_numbers(real time) {
   real t0 = 40.0;     // mean mouse age at BMT for the first ageBMT bin
   real value;
   // spline fitted separately to the counts of thymic FoxP3 negative SP4 T cells to estimate the parameters
   real b0  = 4.3; real b1 = 5.1;  real nu = 40 ;
   //best fitting spline
   value = 10^b0 + (10^b1/(1 + ((time - t0)/nu)^2));
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
    real value;  real chiEst = 0.8;   real qEst = 0.05;
    if (time - 10 < 0){              // t0 = 14 -- assumption: no donor cells seen in FoxP3neg SP4 compartment for 2 weeks
      value = 0;
    } else {
      value = chiEst * (1 - exp(-qEst * (time - 10)));
    }
    return value;
  }

real[] shm_chi(real time, real[] y, real[] parms, real[] rdata,  int[] idata) {
  real psi = parms[1];
  real rho_D = parms[2];
  real delta_D = parms[3];
  real rho_I = parms[4];

  real dydt[6];
  real kloss  = 1/3.5;            //rate of loss of ki67
  real eps_host = 0.15;     // Mean Ki67 hi fraction in host-BM-derived FoxP3 positive Sp4 T cells
  real eps_donor = 0.09;     // Mean Ki67 hi fraction in donor-BM-derived FoxP3 positive Sp4 T cells

  // age of BMT in each recipient
  real ageAtBMT = parms[5];

  // Host naive Tregs
  // ki hi displaceable
  dydt[1] = theta_spline(time, psi) * (1- Chi_spline(time - ageAtBMT)) * eps_host + rho_D * (2 * y[2] + y[1]) - (kloss + delta_D) * y[1];
  // ki lo displaceable
  dydt[2] = theta_spline(time, psi) * (1- Chi_spline(time - ageAtBMT)) * (1 - eps_host) + kloss * y[1]  - (rho_D + delta_D) * y[2];
  // Peripheral ki hi Incumbent
  dydt[3] = rho_I * (2 * y[4] + y[3]) - (kloss + rho_I) * y[3];
  // Peripheral ki lo Incumbent
  dydt[4] = kloss * y[3] - (rho_I + rho_I) * y[4];

  // Donor naive Tregs
  // ki hi displaceable
  dydt[5] = theta_spline(time, psi) * Chi_spline(time - ageAtBMT) * eps_host + rho_D * (2 * y[6] + y[5]) - (kloss + delta_D) * y[5];
  // ki lo displaceable
  dydt[6] = theta_spline(time, psi) * Chi_spline(time - ageAtBMT) * (1 - eps_host) + kloss * y[5]  - (rho_D + delta_D) * y[6];
  
  return dydt;
}

// solving for total counts of thymic and peripheral naive Tregs at time of BMT assuming the youngest animal as the t0
// these counts form the initial conditions for other recipients N(0) = N_host(0)
real[] solve_init(real ageAtBMT,
  real[] init_cond,                          // initial conditions at BMT in the youngest mouse
  real[] parms){

    real ta = 40;                            // age at BMT for the youngest host
    real y_init[2, 6];
    real params_init[5];

    params_init[1:4] = parms[1:4];
    params_init[5] = ta;

    y_init[1] = init_cond;                   // init conditions at the earliest BMT (i.e. in younegst animal)
    if (ageAtBMT==40) {
      y_init[2] = init_cond;
    } else {
      y_init[2] = to_array_1d(integrate_ode_rk45(shm_chi, init_cond, ta, rep_array(ageAtBMT, 1), params_init, {0.0}, {0}));
    }
    return y_init[2];
}

real[] solve_chi(real solve_time, real ageAtBMT, real[] init_cond, real[] parms){
     real y_solve[6];
    real params[5];

    real y0[6];
    real init_tb[6];                         // init conditions at the mean age of BMT for the group

    //solution for the initial conditions at the mean age of BMT for the group
    y0 = solve_init(ageAtBMT, init_cond, parms);

    // init conditions at the BMT
    init_tb[1] = y0[1] + y0[5];
    init_tb[2] = y0[2] + y0[6];
    init_tb[3] = y0[3];                               //at tbmt - all cells are host
    init_tb[4] = y0[4];
    init_tb[5] = 0.0;
    init_tb[6] = 0.0;

    params[1:4] = parms[1:4];
    params[5] = ageAtBMT;                                           // age at BMT

    y_solve = to_array_1d(integrate_ode_rk45(shm_chi, init_tb, ageAtBMT, rep_array(solve_time, 1), params, {0.0}, {0}));

    return y_solve;
  }

 real[,] solve_ode_chi(real[] solve_time, real[] ageAtBMT, real[] init_cond, real[] parms){
    int numdim = size(solve_time);
    real y_solve[numdim, 6];
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
    real y1_0 = global_params[5]; real y2_0 = global_params[6];  real y3_0 = global_params[7];
    real y4_0 = global_params[8]; 

    real init_cond[6];

    // ODE solution -- predictions for the observed timecourse
    real chi_solve[n, 6];

    real counts_naive[n]; real donor_counts_naive[n]; real host_counts_naive[n]; real donor_ki_naive[n]; real host_ki_naive[n];

    vector[4*n] y_mean_stacked;

    // ODE solution -- predictions for the observed timecourse
    init_cond[1] = y1_0; init_cond[2] = y2_0; init_cond[3] = y3_0; init_cond[4] = y4_0;
    init_cond[5] = 0.0; init_cond[6] = 0.0; 

    for (i in 1:n){
      tb_time[i] = ageAtBMT[i]/1.0;
    }

    // each shard has a single datpoint so its unique ****
    // PDE solution for chimera dataset -- x_r = data time and x_i = time at BMT
      chi_solve = solve_ode_chi(solve_time, tb_time, init_cond, to_array_1d(global_params));

      for (i in 1:n){
        counts_naive[i] = chi_solve[i, 1] + chi_solve[i, 2] + chi_solve[i, 3] + chi_solve[i, 4] + chi_solve[i, 5] + chi_solve[i, 6];
        donor_counts_naive[i] = chi_solve[i, 5] + chi_solve[i, 6];
        host_counts_naive[i] = chi_solve[i, 1] + chi_solve[i, 2] + chi_solve[i, 3] + chi_solve[i, 4];
        donor_ki_naive[i] = (chi_solve[i, 5])/donor_counts_naive[i];
        host_ki_naive[i] = (chi_solve[i, 1] + chi_solve[i, 3])/host_counts_naive[i];

        y_mean_stacked[4*i - 3] = counts_naive[i];
        y_mean_stacked[4*i - 2] = donor_counts_naive[i]/(counts_naive[i] * Chi_spline(solve_time[i] - tb_time[i]));
        y_mean_stacked[4*i - 1] = donor_ki_naive[i];
        y_mean_stacked[4*i - 0] = host_ki_naive[i];
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
   real<lower = 0> counts_naive[numObs1];
   real<lower = 0> Nfd_naive[numObs2];
   real<lower = 0> ki_donor_naive[numObs3];
   real<lower = 0> ki_host_naive[numObs4];
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
  real<lower= 0, upper= 1> delta_D;
  real<lower= 0, upper= 1> rho_I;
  real<lower= 0> y1_0;
  real<lower= 0> y2_0;
  real<lower= 0> y3_0;
  real<lower= 0> y4_0;


  real<lower=0> sigma_counts_naive;
  real<lower=0> sigma_Nfd_naive;
  real<lower=0> sigma_donor_ki_naive;
  real<lower=0> sigma_host_ki_naive;
}

transformed parameters{
  vector[8] global_params;
  vector[n_solve] counts_naive_solve;               // ODE predictions for naive Treg counts in thymus
  vector[n_solve] Nfd_naive_solve;                  // ODE predictions for naive Treg Nfd in thymus
  vector[n_solve] ki_host_naive_solve;              // ODE predictions for naive Treg host ki67 proportions in thymus
  vector[n_solve] ki_donor_naive_solve;             // ODE predictions for naive Treg donor ki67 proportions in periphery
  vector[4*n_solve] y_mean_stacked;               // compliled output across all nodes

  vector[numObs1] counts_naive_mean;               // ODE predictions for naive Treg counts in thymus
  vector[numObs2] Nfd_naive_mean;                  // ODE predictions for naive Treg Nfd in thymus
  vector[numObs3] ki_donor_naive_mean;             // ODE predictions for naive Treg donor ki67 proportions in thymus
  vector[numObs4] ki_host_naive_mean;              // ODE predictions for naive Treg host ki67 proportions in thymus

  global_params[1] = psi;
  global_params[2] = rho_D;
  global_params[3] = delta_D;
  global_params[4] = rho_I;
  global_params[5] = exp(y1_0);
  global_params[6] = exp(y2_0);
  global_params[7] = exp(y3_0);
  global_params[8] = exp(y4_0);

  // combining the output from all the shards
  y_mean_stacked = map_rect(math_reduce, global_params, local_params, x_r, x_i);

  for (i in 1:n_solve){
    counts_naive_solve[i] = y_mean_stacked[4*i - 3];
    Nfd_naive_solve[i] = y_mean_stacked[4*i - 2];
    ki_donor_naive_solve[i] = y_mean_stacked[4*i - 1];
    ki_host_naive_solve[i] = y_mean_stacked[4*i - 0];
  }

  for (i in 1:numObs1){
    counts_naive_mean[i]   = counts_naive_solve[time_index_counts[i]];
  }
  for (i in 1:numObs2){
    Nfd_naive_mean[i]   = Nfd_naive_solve[time_index_chi[i]];
  }
  for (i in 1:numObs3){
    ki_donor_naive_mean[i]   = ki_donor_naive_solve[time_index_donorki[i]];
  }
  for (i in 1:numObs4){
    ki_host_naive_mean[i]   = ki_host_naive_solve[time_index_hostki[i]];
  }
}

model{
  psi ~ normal(0.3, 0.2);
  rho_D ~ normal(0.005, 0.25);
  delta_D ~ normal(0.01, 0.25);
  rho_I ~ normal(0.01, 0.25);
  y1_0 ~ normal(9, 2.5);
  y2_0 ~ normal(11, 2.5);
  y3_0 ~ normal(9, 2.5);
  y4_0 ~ normal(11, 2.5);

  sigma_counts_naive ~ normal(0.5, 0.1);
  sigma_Nfd_naive ~ normal(0.4, 0.05);
  sigma_donor_ki_naive ~ normal(0.3, 0.5);
  sigma_host_ki_naive ~ normal(0.2, 2);

  log(counts_naive) ~ normal(log(counts_naive_mean), sigma_counts_naive);
  asinsqrt_array(Nfd_naive) ~ normal(asinsqrt_array(to_array_1d(Nfd_naive_mean)), sigma_Nfd_naive);
  asinsqrt_array(ki_donor_naive) ~ normal(asinsqrt_array(to_array_1d(ki_donor_naive_mean)), sigma_donor_ki_naive);
  asinsqrt_array(ki_host_naive) ~ normal(asinsqrt_array(to_array_1d(ki_host_naive_mean)), sigma_host_ki_naive);
}

generated quantities{
  real y_chi_pred1[numPred, 6];
  real y_chi_pred2[numPred, 6];
  real y_chi_pred3[numPred, 6];
  real y_chi_pred4[numPred, 6];

  real counts_naive_mean_pred1[numPred];    real counts_naive_pred1[numPred]; 
  real counts_naive_mean_pred2[numPred];    real counts_naive_pred2[numPred]; 
  real counts_naive_mean_pred3[numPred];    real counts_naive_pred3[numPred]; 
  real counts_naive_mean_pred4[numPred];    real counts_naive_pred4[numPred]; 

  real Nfd_naive_mean_pred1[numPred];     real Nfd_naive_pred1[numPred];       
  real Nfd_naive_mean_pred2[numPred];     real Nfd_naive_pred2[numPred];    
  real Nfd_naive_mean_pred3[numPred];     real Nfd_naive_pred3[numPred];    
  real Nfd_naive_mean_pred4[numPred];     real Nfd_naive_pred4[numPred];  

  real ki_donor_naive_mean_pred1[numPred];     real ki_donor_naive_pred1[numPred];  
  real ki_donor_naive_mean_pred2[numPred];     real ki_donor_naive_pred2[numPred];  
  real ki_donor_naive_mean_pred3[numPred];     real ki_donor_naive_pred3[numPred];  
  real ki_donor_naive_mean_pred4[numPred];     real ki_donor_naive_pred4[numPred]; 
   

  real ki_host_naive_mean_pred1[numPred];   real ki_host_naive_pred1[numPred]; 
  real ki_host_naive_mean_pred2[numPred];   real ki_host_naive_pred2[numPred]; 
  real ki_host_naive_mean_pred3[numPred];   real ki_host_naive_pred3[numPred]; 
  real ki_host_naive_mean_pred4[numPred];   real ki_host_naive_pred4[numPred]; 

  // log likelihoods
  vector[numObs1] log_lik_counts_naive;
  vector[numObs2] log_lik_Nfd_naive;
  vector[numObs3] log_lik_ki_donor_naive;
  vector[numObs4] log_lik_ki_host_naive;

  // initial conditions
  real init_cond[6];
  init_cond[1] = exp(y1_0); init_cond[2] = exp(y2_0); init_cond[3] = exp(y3_0); init_cond[4] = exp(y4_0);
  init_cond[5] = 0.0; init_cond[6] = 0.0; 

  // predictions for the whole timecourse
  y_chi_pred1 = solve_ode_chi(ts_pred1, tb_pred1, init_cond, to_array_1d(global_params));
  y_chi_pred2 = solve_ode_chi(ts_pred2, tb_pred2, init_cond, to_array_1d(global_params));
  y_chi_pred3 = solve_ode_chi(ts_pred3, tb_pred3, init_cond, to_array_1d(global_params));
  y_chi_pred4 = solve_ode_chi(ts_pred4, tb_pred4, init_cond, to_array_1d(global_params));

  for (i in 1:numPred){
    counts_naive_mean_pred1[i] = y_chi_pred1[i, 1] + y_chi_pred1[i, 2] + y_chi_pred1[i, 3] + y_chi_pred1[i, 4] + y_chi_pred1[i, 5] + y_chi_pred1[i, 6];
    counts_naive_mean_pred2[i] = y_chi_pred2[i, 1] + y_chi_pred2[i, 2] + y_chi_pred2[i, 3] + y_chi_pred2[i, 4] + y_chi_pred2[i, 5] + y_chi_pred2[i, 6];
    counts_naive_mean_pred3[i] = y_chi_pred3[i, 1] + y_chi_pred3[i, 2] + y_chi_pred3[i, 3] + y_chi_pred3[i, 4] + y_chi_pred3[i, 5] + y_chi_pred3[i, 6];
    counts_naive_mean_pred4[i] = y_chi_pred4[i, 1] + y_chi_pred4[i, 2] + y_chi_pred4[i, 3] + y_chi_pred4[i, 4] + y_chi_pred4[i, 5] + y_chi_pred4[i, 6];

    counts_naive_pred1[i] = exp(normal_rng(log(counts_naive_mean_pred1[i]), sigma_counts_naive));
    counts_naive_pred2[i] = exp(normal_rng(log(counts_naive_mean_pred2[i]), sigma_counts_naive));
    counts_naive_pred3[i] = exp(normal_rng(log(counts_naive_mean_pred3[i]), sigma_counts_naive));
    counts_naive_pred4[i] = exp(normal_rng(log(counts_naive_mean_pred4[i]), sigma_counts_naive));

    Nfd_naive_mean_pred1[i] = (y_chi_pred1[i, 5] + y_chi_pred1[i, 6])/(counts_naive_mean_pred1[i] * Chi_spline(ts_pred1[i] - tb_pred1[i]));
    Nfd_naive_mean_pred2[i] = (y_chi_pred2[i, 5] + y_chi_pred2[i, 6])/(counts_naive_mean_pred2[i] * Chi_spline(ts_pred2[i] - tb_pred2[i]));
    Nfd_naive_mean_pred3[i] = (y_chi_pred3[i, 5] + y_chi_pred3[i, 6])/(counts_naive_mean_pred3[i] * Chi_spline(ts_pred3[i] - tb_pred3[i]));
    Nfd_naive_mean_pred4[i] = (y_chi_pred4[i, 5] + y_chi_pred4[i, 6])/(counts_naive_mean_pred4[i] * Chi_spline(ts_pred4[i] - tb_pred4[i]));
   
    Nfd_naive_pred1[i] = asinsqrt_inv(normal_rng(asinsqrt_real(Nfd_naive_mean_pred1[i]), sigma_Nfd_naive));
    Nfd_naive_pred2[i] = asinsqrt_inv(normal_rng(asinsqrt_real(Nfd_naive_mean_pred2[i]), sigma_Nfd_naive));
    Nfd_naive_pred3[i] = asinsqrt_inv(normal_rng(asinsqrt_real(Nfd_naive_mean_pred3[i]), sigma_Nfd_naive));
    Nfd_naive_pred4[i] = asinsqrt_inv(normal_rng(asinsqrt_real(Nfd_naive_mean_pred4[i]), sigma_Nfd_naive));

    ki_donor_naive_mean_pred1[i] = (y_chi_pred1[i, 5])/(y_chi_pred1[i, 5] + y_chi_pred1[i, 6]);
    ki_donor_naive_mean_pred2[i] = (y_chi_pred2[i, 5])/(y_chi_pred2[i, 5] + y_chi_pred2[i, 6]);
    ki_donor_naive_mean_pred3[i] = (y_chi_pred3[i, 5])/(y_chi_pred3[i, 5] + y_chi_pred3[i, 6]);
    ki_donor_naive_mean_pred4[i] = (y_chi_pred4[i, 5])/(y_chi_pred4[i, 5] + y_chi_pred4[i, 6]);
    
    ki_donor_naive_pred1[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_donor_naive_mean_pred1[i]), sigma_donor_ki_naive));
    ki_donor_naive_pred2[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_donor_naive_mean_pred2[i]), sigma_donor_ki_naive));
    ki_donor_naive_pred3[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_donor_naive_mean_pred3[i]), sigma_donor_ki_naive));
    ki_donor_naive_pred4[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_donor_naive_mean_pred4[i]), sigma_donor_ki_naive));

    ki_host_naive_mean_pred1[i] = (y_chi_pred1[i, 1] + y_chi_pred1[i, 3])/(y_chi_pred1[i, 1] + y_chi_pred1[i, 2] + y_chi_pred1[i, 3] + y_chi_pred1[i, 4]);
    ki_host_naive_mean_pred2[i] = (y_chi_pred2[i, 1] + y_chi_pred2[i, 3])/(y_chi_pred2[i, 1] + y_chi_pred2[i, 2] + y_chi_pred2[i, 3] + y_chi_pred2[i, 4]);
    ki_host_naive_mean_pred3[i] = (y_chi_pred3[i, 1] + y_chi_pred3[i, 3])/(y_chi_pred3[i, 1] + y_chi_pred3[i, 2] + y_chi_pred3[i, 3] + y_chi_pred3[i, 4]);
    ki_host_naive_mean_pred4[i] = (y_chi_pred4[i, 1] + y_chi_pred4[i, 3])/(y_chi_pred4[i, 1] + y_chi_pred4[i, 2] + y_chi_pred4[i, 3] + y_chi_pred4[i, 4]);
    
    ki_host_naive_pred1[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_host_naive_mean_pred1[i]), sigma_host_ki_naive));
    ki_host_naive_pred2[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_host_naive_mean_pred2[i]), sigma_host_ki_naive));
    ki_host_naive_pred3[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_host_naive_mean_pred3[i]), sigma_host_ki_naive));
    ki_host_naive_pred4[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_host_naive_mean_pred4[i]), sigma_host_ki_naive));
  }

  // calculating log likelihoods
  for (i in 1:numObs1) {
    log_lik_counts_naive[i] = normal_lpdf(log(counts_naive[i]) | log(counts_naive_mean[i]), sigma_counts_naive);
  }
  for (i in 1:numObs2) {
    log_lik_Nfd_naive[i] = normal_lpdf(asinsqrt_real(Nfd_naive[i]) | asinsqrt_real(Nfd_naive_mean[i]), sigma_Nfd_naive);
  }
  for (i in 1:numObs3) {
    log_lik_ki_donor_naive[i] = normal_lpdf(asinsqrt_real(ki_donor_naive[i]) | asinsqrt_real(ki_donor_naive_mean[i]), sigma_donor_ki_naive);
  }
  for (i in 1:numObs4) {
    log_lik_ki_host_naive[i] = normal_lpdf(asinsqrt_real(ki_host_naive[i]) | asinsqrt_real(ki_host_naive_mean[i]), sigma_host_ki_naive);
  }
}
