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

  // initial age distribution of cells exisiting at t0
  real g_age(real age, real[] parms) {
    real t0 = 1.0;
    real N0 = parms[1];
    real value;

    // Flat, normalised age-distribution of initial cells
    if(age >= 0 && age <= t0) {
      value = N0/t0;
    } else {
      value = 0;
    }
    return value;
  }

  // Total influx into the naive T cell compartment from the thymus (cells/day)
  real theta_spline(real time, real[] parms){
    real N0 = parms[1];
    real t0 = 1.0;
    real psi;  real value;

    psi = g_age(0.0, parms)/sp_numbers(t0);
    value = psi * sp_numbers(time);
    return value;
  }

  real Chi_spline( real time) {
    // chiEst is the level if stabilised chimerism in the source compartment
    // qEst is the rate with which cimerism chnages in the source compartment
    real chi;
    real chiEst = 0.83155297;
    real qEst = 0.047555868;

    if (time < 0){
      chi = 0;                       // conditioning the function to adapt to the timepoints before BMT
    } else {
      chi = chiEst * (1 - exp(-qEst * time));
    }
    return chi;
  }


  // influx of donor cells into the naive donor T cell compartment from the thymus (cells/day)
  real theta_donor(real time, real[] parms){
    real value;
    real tBMT = parms[4];
    real tC = 42;
    real m = 2;

    //value = theta_spline(time, parms) * (time/tC)^m;
    if (time <= tC) {
      value = theta_spline(time, parms) * (time/tC)^m;
    } else {
      value = theta_spline(time, parms);
    }
    return value;
  }

  // influx of host derived cells into the naive host T cell compartment from the thymus (cells/day)
  real theta_host(real time, real[] parms){

    return theta_spline(time, parms) - theta_donor(time, parms);
  }

  // Rate of loss -- varies with cell age
  // integrand function
  real[] lambda_age(real age, real[] y, real[] parms, real[] rdata, int[] idata) {
    real lambda0 = parms[2];
    real r_lambda = parms[3];

    real value = lambda0 * exp(- r_lambda * age);
    return {value};
  }

  // function that calculates intgral of net loss rate --  solved analytically to speed up the numerical integration
  real lambda_integ(real lo_lim, real up_lim, real[] parms){
    real lambda0 = parms[2];
    real r_lambda = parms[3];

    real value = (lambda0/r_lambda) * (exp(-r_lambda * lo_lim) - exp(-r_lambda * up_lim));
    return value;
  }

  // Cell age distribution of the initial cohort
  real Asm_init_age(real age, real time, real[] parms) {
    real t0 = 1.0;

    real value = g_age(t0, parms) * exp(- lambda_integ(age - time + t0, age, parms));
    return value;
   }

   // Cell age distribution of the total theta cohort
   real Asm_theta_age(real age, real time, real[] parms) {

     real value =  theta_spline(time - age, parms) * exp(- lambda_integ(0, age, parms));
     return value;
    }

  // Cell age distribution of the whole population
  real Asm_total_age(real age, real time, real[] parms){
    real value;
    real t0 = 1.0;

    if(age < (time - t0)) {
      value =  theta_spline(time - age, parms) * exp(- lambda_integ(0, age, parms));
      } else {
        value = g_age(t0, parms) * exp(- lambda_integ(age - time + t0, age, parms));
    }

    return value;
  }

  // Cell age distribution of the initial cohort
  real Asm_Host_init_age(real age, real time, real[] parms) {
    real tBMT = parms[4];
    real value;

    if (age >= time - tBMT){
      value = Asm_total_age(age - time + tBMT, tBMT, parms) * exp(- lambda_integ(age - time + tBMT, age, parms));
    } else {
      value = 0.0;
    }
    return value;
   }

   // Cell age distribution of the host theta cohort
   real Asm_Host_theta_age(real age, real time, real[] parms) {
     real tBMT = parms[4];
     real value;

     if (age < time - tBMT){
       value = theta_host(time - age, parms) * exp(- lambda_integ(0, age, parms));
     } else {
       value = 0.0;
     }
     return value;
   }

   // Cell age distribution of the donor theta cohort
   real Asm_Donor_theta_age(real age, real time, real[] parms) {
     real tBMT = parms[4];
     real value;

     if (age < time - tBMT){
       value = theta_donor(time - age, parms) * exp(- lambda_integ(0, age, parms));
     } else {
       value = 0.0;
     }
     return value;
   }

   // Cell age distribution of the total theta cohort
   real Asm_pooled_age(real age, real time, real[] parms) {
     real tBMT = parms[4];
     real value;

     if (age < time  - tBMT){
       value = Asm_Host_theta_age(age, time,  parms) + Asm_Donor_theta_age(age, time, parms);
     } else {
       value = Asm_Host_init_age(age, time, parms);
     }

     return value;
   }

   // Cell age distribution of the host cohort
   real Asm_host_age(real age, real time, real[] parms) {
     real tBMT = parms[4];
     real value;

     if (age < time  - tBMT){
       value = Asm_Host_theta_age(age, time,  parms);
     } else {
       value = Asm_Host_init_age(age, time, parms);
     }

     return value;
   }

   // Cell age distribution of the donor cohort
   real Asm_donor_age(real age, real time, real[] parms) {
     real tBMT = parms[4];
     real value;

     if (age < time - tBMT){
       value = Asm_Donor_theta_age(age, time, parms);
     } else {
       value = 0.0;
     }

     return value;
   }

  real[] Asm_total_ode(real age,  real[] y, real[] parms, real[] x_r,  int[] x_i) {
    real value;
    real time = parms[4];   // time (host age) input  as a param

    value = Asm_total_age(age, time, parms);

    return {value};
  }

  real[] Asm_pooled_ode(real age,  real[] y, real[] parms, real[] x_r,  int[] x_i) {
    real value;
    real time = parms[5];   // time (host age) input  as a param

    value = Asm_pooled_age(age, time, parms);

    return {value};
  }

  real[] Asm_host_ode(real age,  real[] y, real[] parms, real[] x_r,  int[] x_i) {
    real value;
    real time = parms[5];   // time (host age) input  as a param

    value = Asm_host_age(age, time, parms);

    return {value};
  }

  real[] Asm_donor_ode(real age,  real[] y, real[] parms, real[] x_r,  int[] x_i) {
    real value;
    real time = parms[5];   // time (host age) input  as a param

    value = Asm_donor_age(age, time, parms);

    return {value};
  }

  real solve_total_counts(real[] parms) {
    int x_i[0];
    real value;
    real time = parms[4];   // time (host age) input  as a param

    // integrate_ode_rk45(function, y0, t0, t, theta, x_r, x_i);
    value = integrate_ode_rk45(Asm_total_ode, {0.0}, 0.0, rep_array(time, 1), parms, {0.0}, x_i)[1, 1];
    return value;
  }

  real solve_pooled_counts(real[] parms) {
    int x_i[0];
    real value;
    real time = parms[5];   // time (host age) input  as a param

    value = integrate_ode_rk45(Asm_pooled_ode, {0.0}, 0.0, rep_array(time, 1), parms, {0.0}, x_i)[1, 1];
    return value;
  }

  real solve_host_counts(real[] parms) {
    int x_i[0];
    real value;
    real time = parms[5];   // time (host age) input  as a param

    value = integrate_ode_rk45(Asm_host_ode, {0.0}, 0.0, rep_array(time, 1), parms, {0.0}, x_i)[1, 1];
    return value;
  }

  real solve_donor_counts(real[] parms) {
    int x_i[0];
    real value;
    real time = parms[5];   // time (host age) input  as a param

    value = integrate_ode_rk45(Asm_donor_ode, {0.0}, 0.0, rep_array(time, 1), parms, {0.0}, x_i)[1, 1];
    return value;
  }

  // Vectorised function for total pool size
  real[] N_total_time(real[] time, real[] parms){
   int ndim = size(time);
   real y_solve[ndim];
   real params[5];
   params[1:3] = parms[1:3];

   for (i in 1:ndim){
     params[4] = time[i];
     y_solve[i] = solve_total_counts(params);
   }
   return y_solve;
  }

  // Vectorised function for total pool size
  real[] N_pooled_time(real[] time, real[] tBMT, real[] parms){
    int ndim = size(time);
    real y_solve[ndim];
    real params[5];
    params[1:3] = parms[1:3];

    for (i in 1:ndim){
      params[4] = tBMT[i];
      params[5] = time[i];
      y_solve[i] = solve_pooled_counts(params);
    }
    return y_solve;
  }

  // Vectorised function for total pool size
  real[] N_host_time(real[] time, real[] tBMT, real[] parms){
    int ndim = size(time);
    real y_solve[ndim];
    real params[5];
    params[1:3] = parms[1:3];

    for (i in 1:ndim){
      params[4] = tBMT[i];
      params[5] = time[i];
      y_solve[i] = solve_host_counts(params);
    }
    return y_solve;
  }

  // Vectorised function for total pool size
  real[] N_donor_time(real[] time, real[] tBMT, real[] parms){
    int ndim = size(time);
    real y_solve[ndim];
    real params[5];
    params[1:3] = parms[1:3];

    for (i in 1:ndim){
      params[4] = tBMT[i];
      params[5] = time[i];
      y_solve[i] = solve_donor_counts(params);
    }
    return y_solve;
  }


  vector math_reduce(vector global_params, vector local_params, real[] x_r, int[] x_i){
    // data for each shard
    int n = size(x_i); // n = 1
    int dat_t0 = x_i[1];                          // time zero -- for chimeras age at BM

    real chi_counts_mean[n];
    real host_counts_mean[n];
    real donor_counts_mean[n];

    vector[2*n] y_mean_stacked;

    // each shard has a single datpoint so its unique ****
    // PDE solution for chimera dataset -- x_r = data time and x_i = time at BMT
    chi_counts_mean = N_pooled_time(x_r,  to_array_1d(to_vector(x_i)/1.0), to_array_1d(global_params));
    host_counts_mean = N_host_time(x_r, to_array_1d(to_vector(x_i)/1.0), to_array_1d(global_params));
    donor_counts_mean = N_donor_time(x_r, to_array_1d(to_vector(x_i)/1.0), to_array_1d(global_params));

    y_mean_stacked[1] = chi_counts_mean[1];
    y_mean_stacked[2] = donor_counts_mean[1]/chi_counts_mean[1];

    return y_mean_stacked;
  }

  // functions for transformation of fractions in (0,a), where a >=1
  real logit_inverse(real x){
     real ans;

     ans = exp(x)/(1+exp(x));

     return ans;
   }
}

data{
  int<lower = 0> numOnt;
  int<lower = 1> numChi;
  real<lower = 0> chi_counts[numChi];
  real<lower = 0> N_donor_fraction[numChi];
  int<lower  = 1> numPred;
  real<lower = 0> ts_pred_ont[numPred];
  real<lower = 0> ts_pred_chi1[numPred];
  real<lower = 0> ts_pred_chi2[numPred];
  real<lower = 0> ts_pred_chi3[numPred];
  real<lower = 0> tb_pred1[numPred];
  real<lower = 0> tb_pred2[numPred];
  real<lower = 0> tb_pred3[numPred];
  int n_shards;
  int<lower = 0> dat_t0[n_shards];       // nshards = numOnt + numChi
  int<lower = 0> dat_time[n_shards];     // nshards = numOnt + numChi
}

transformed data{
  int x_i[n_shards, 1];         // each shard gets a single data point
  real x_r[n_shards, 1];        // each shard gets a single data point

  // empty set of per shard params
  vector[0] local_params[n_shards];  // shard specific params --  useful for hierarchical modelling

  // data split into shards
  for (s in 1:n_shards){
   x_i[s, 1] = dat_t0[s];                       // age at BMT split
   x_r[s, 1] = dat_time[s];                     // time split
  }
}

parameters{
  real<lower=1E4, upper=2E6> N0;                  // total cells counts at t0
  real<lower=0.001, upper=0.5> lambda0;
  real r_lambda;
  real<lower=0> sigma_chi_counts;
  real<lower=0> sigma_Nfd;
}

transformed parameters{
  vector[3] global_params;
  vector[numChi] y3_mean;               // PDE prediction for counts from chimera data
  vector[numChi] y4_mean;               // PDE prediction for Nfd from chimera data
  vector[(2*numChi)] y_mean_stacked;        // compliled output across all nodes

  global_params[1] = N0;
  global_params[2] = lambda0;
  global_params[3] = r_lambda;

  // combining the output from all the shards
  y_mean_stacked = map_rect(math_reduce, global_params, local_params, x_r, x_i);

  for (i in 1:numChi){
    y3_mean[i] = y_mean_stacked[(4*numOnt)+ 2*i - 1];
    y4_mean[i] = y_mean_stacked[(4*numOnt)+ 2*i - 0];
  }
}

model{
  N0 ~ normal(5E5, 1.5E5);
  lambda0 ~ normal(0.05, 0.2);
  r_lambda ~ normal(0.0, 0.2);

  sigma_chi_counts ~ normal(0, 2);
  sigma_Nfd ~ normal(0, 2);

  log(chi_counts) ~ normal(log(y3_mean), sigma_chi_counts);
  logit(N_donor_fraction) ~ normal(logit(y4_mean), sigma_Nfd);
}

generated quantities{
  real y_chi_pred1[numPred, 2];
  real y_chi_pred2[numPred, 2];
  real y_chi_pred3[numPred, 2];

  real y3_mean_pred1[numPred];  real y4_mean_pred1[numPred];
  real y3_mean_pred2[numPred];  real y4_mean_pred2[numPred];
  real y3_mean_pred3[numPred];  real y4_mean_pred3[numPred];

  real chicounts_pred1[numPred]; real Nfd_pred1[numPred];
  real chicounts_pred2[numPred]; real Nfd_pred2[numPred];
  real chicounts_pred3[numPred]; real Nfd_pred3[numPred];

  real host_counts_pred1[numPred]; real host_counts_pred2[numPred];  real host_counts_pred3[numPred];
  real donor_counts_pred1[numPred]; real donor_counts_pred2[numPred];  real donor_counts_pred3[numPred];

  // log likelihoods
  vector[numChi] log_lik_chi_counts;
  vector[numChi] log_lik_Nfd;

  // PDE solution -- predictions for total counts, Nfd
  y3_mean_pred1 = N_pooled_time(ts_pred_chi1,  tb_pred1, to_array_1d(global_params));
  y3_mean_pred2 = N_pooled_time(ts_pred_chi2,  tb_pred2, to_array_1d(global_params));
  y3_mean_pred3 = N_pooled_time(ts_pred_chi3,  tb_pred3, to_array_1d(global_params));

  host_counts_pred1 = N_host_time(ts_pred_chi1, tb_pred1, to_array_1d(global_params));
  host_counts_pred2 = N_host_time(ts_pred_chi1, tb_pred2, to_array_1d(global_params));
  host_counts_pred3 = N_host_time(ts_pred_chi1, tb_pred3, to_array_1d(global_params));

  donor_counts_pred1 = N_donor_time(ts_pred_chi1, tb_pred1, to_array_1d(global_params));
  donor_counts_pred2 = N_donor_time(ts_pred_chi2, tb_pred2, to_array_1d(global_params));
  donor_counts_pred3 = N_donor_time(ts_pred_chi3, tb_pred3, to_array_1d(global_params));

  for (i in 1:numPred){
    y4_mean_pred1[i] = donor_counts_pred1[i]/y3_mean_pred1[i];
    y4_mean_pred2[i] = donor_counts_pred2[i]/y3_mean_pred2[i];
    y4_mean_pred3[i] = donor_counts_pred3[i]/y3_mean_pred3[i];

    chicounts_pred1[i] = exp(normal_rng(log(y3_mean_pred1[i]), sigma_chi_counts));
    chicounts_pred2[i] = exp(normal_rng(log(y3_mean_pred2[i]), sigma_chi_counts));
    chicounts_pred3[i] = exp(normal_rng(log(y3_mean_pred3[i]), sigma_chi_counts));

    Nfd_pred1[i] = logit_inverse(normal_rng(logit(y4_mean_pred1[i]), sigma_Nfd));
    Nfd_pred2[i] = logit_inverse(normal_rng(logit(y4_mean_pred2[i]), sigma_Nfd));
    Nfd_pred3[i] = logit_inverse(normal_rng(logit(y4_mean_pred3[i]), sigma_Nfd));
}

  // calculating log likelihoods
  for (i in 1:numChi) {
    log_lik_chi_counts[i] = normal_lpdf(log(chi_counts[i]) | log(y3_mean[i]), sigma_chi_counts);
    log_lik_Nfd[i]        = normal_lpdf(logit(N_donor_fraction[i]) | logit(y4_mean[i]), sigma_Nfd);
  }
}
