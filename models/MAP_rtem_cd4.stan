 functions{
 //spline 1
  // Timecourse of thymic CD4 SP population -- changes with time
  real sp_numbers(real time) {
    real t0 = 1.0;
    real dpt0 = time - t0;     // days post t0
    real value; real fit1;
    // spline fitted separately to the counts of thymic SP4 cells
    // parameters estimated from spline fit to the timecourse of counts of source compartment -- SP CD4
    real theta0  =  4.3E5;  real theta_f = 1.8E3;  real n = 2.1;   real X = 30.0;   real q = 3.7;
    //best fitting spline
    fit1 = theta0 + (theta_f * dpt0^n) * (1 - ((dpt0^q)/((X^q) + (dpt0^q))));

    if(time < t0){
      value = 0.0;
    } else {
      value = fit1;
    }
    return value;
  }

  // Total influx into the naive T cell compartment from the thymus (cells/day)
  real theta_spline(real time, real psi){
    real value;

    value = psi * sp_numbers(time);
    return value;
  }

  // spline2 --
   // proportions of ki67 hi cells in source -- varies with time
   real eps_spline(real time){
     real value;
     //parameters estimated from spline fit to the timecourse of ki67 proportions of source compartment -- SP CD4
     real eps_0 = 0.14965320; real eps_f = 0.03470231; real A = 3.43078629;
     //best fitting spline
     real fit = exp(-eps_f * (time + A)) + eps_0;

     return fit;
  }

  real Chi_spline( real time) {
    // chiEst is the level if stabilised chimerism in the source compartment
    // qEst is the rate with which cimerism chnages in the source compartment
    real chi;
    real chiEst = 0.847543332;
    real qEst = 0.050944623;

    if (time < 0){
      chi = 0;                       // conditioning the function to adapt to the timepoints before BMT
    } else {
      chi = chiEst * (1 - exp(-qEst * time));
    }
    return chi;
  }

  real[] shm_ont(real t, real[] y, real[] parms, real[] rdata,  int[] idata) {
    real psi       = parms[1];
    real delta_nai = parms[2];
    real delta_rte = parms[3];
    real rho_nai   = parms[4];
    real rho_rte   = parms[5];
    real mu        = parms[6];

    real dydt[4];
    real beta  = 1/3.5;            //rate of loss of ki67 -- mature naive cells

    // model that assumes that RTEs divide and die at different rates than mature naive T cells
    // ki hi RTE
    dydt[1] = theta_spline(t, psi) * eps_spline(t) + rho_rte * (2 * y[2] + y[1]) - (beta + delta_rte + mu) * y[1];
    // ki lo RTE
    dydt[2] = theta_spline(t, psi) * (1 - eps_spline(t)) + beta * y[1] - (rho_rte + delta_rte + mu) * y[2];
    // ki hi mN
    dydt[3] = mu * y[1] + rho_nai * (2 * y[4] + y[3]) - (beta + delta_nai) * y[3];
    // ki lo mN
    dydt[4] = mu * y[2] + beta * y[3] - (rho_nai + delta_nai) * y[4];

    return dydt;
  }

  real[] solve_ont(real solve_time, real[] init_cond, real[] parms) {
    // solves the ode for each timepoint from t0
    return to_array_1d(integrate_ode_rk45(shm_ont, init_cond, 1.0, rep_array(solve_time, 1), parms, {0.0}, {0}));
   }

  real[,] solve_ode_ont(real[] solve_time, real[] init_cond, real[] parms){
    int num_solve = size(solve_time);
    real y_hat[num_solve, 4];
    // ode solution for the whole timecourse
    y_hat[1] = init_cond;
    for (i in 2:num_solve){
      y_hat[i] = solve_ont(solve_time[i], init_cond, parms);
    }
  return y_hat;
  }

  real[] shm_chi(real time, real[] y, real[] parms, real[] rdata,  int[] idata) {
    real psi       = parms[1];
    real delta_nai = parms[2];
    real delta_rte = parms[3];
    real rho_nai   = parms[4];
    real rho_rte   = parms[5];
    real mu        = parms[6];

    real dydt[8];
    real beta  = 1/3.5;            //rate of loss of ki67 -- mature naive cells

    // age of BMT in each recipient
    real ageAtBMT = parms[7];

    // ki hi donor RTE
    dydt[1] = theta_spline(time, psi) * Chi_spline(time - ageAtBMT) * eps_spline(time) + rho_rte * (2 * y[2] + y[1]) - (beta + delta_rte + mu) * y[1];
    // ki lo donor RTE
    dydt[2] = theta_spline(time, psi) * Chi_spline(time - ageAtBMT) * (1 - eps_spline(time)) + beta * y[1] - (rho_rte + delta_rte + mu) * y[2];

    // ki hi donor mN
    dydt[3] = mu * y[1] + rho_nai * (2 * y[4] + y[3]) - (beta + delta_nai) * y[3];
    // ki lo donor mN
    dydt[4] = mu * y[2] + beta * y[3] - (rho_nai + delta_nai) * y[4];


    // ki hi host RTE
    dydt[5] = theta_spline(time, psi) * (1-Chi_spline(time - ageAtBMT)) * eps_spline(time) + rho_rte * (2 * y[6] + y[5]) - (beta + delta_rte + mu) * y[5];
    // ki lo host RTE
    dydt[6] = theta_spline(time, psi) * (1-Chi_spline(time - ageAtBMT)) * (1 - eps_spline(time)) +  beta * y[5] - (rho_rte + delta_rte + mu) * y[6];

    // ki hi mN
    dydt[7] = mu * y[5] + rho_nai * (2 * y[6] + y[7]) - (beta + delta_nai) * y[7];
    // ki lo mN
    dydt[8] = mu * y[6] + beta * y[7] - (rho_nai + delta_nai) * y[8];

    return dydt;
  }

  real[] solve_chi(real solve_time,            // time point of observation
    real ageAtBMT,
    real[] init_cond,
    real[] parms){

      real y_solve[8];
      real params[7];

      real y0[4];
      real init_tb[8];                         // init conditions at the mean age of BMT for the group

      //solution for the initial conditions at the mean age of BMT for the group
      y0 = solve_ont(ageAtBMT, init_cond, parms);

      // init conditions at the BMT
      init_tb[1] = 0;                                           //at tbmt - # donor is zero
      init_tb[2] = 0;                                           //at tbmt - # donor is zero
      init_tb[3] = 0;                                           //at tbmt - # donor is zero
      init_tb[4] = 0;                                           //at tbmt - # donor is zero
      init_tb[5] = y0[1];                                       //at tbmt - all ki67Hi rte cells are host
      init_tb[6] = y0[2];                                       //at tbmt - all ki67Lo rte cells are host
      init_tb[7] = y0[3];                                       //at tbmt - all ki67Hi mN cells are host
      init_tb[8] = y0[4];                                       //at tbmt - all ki67Lo mN cells are host

      params[1:6] = parms[1:6];
      params[7] = ageAtBMT;                                           // age at BMT

      y_solve = to_array_1d(integrate_ode_rk45(shm_chi, init_tb, ageAtBMT, rep_array(solve_time, 1), params, {0.0}, {0}));

      return y_solve;
    }

    real[,] solve_ode_chi(real[] solve_time,
      real[] ageAtBMT,
      real[] init_cond,
      real[] parms){

        int numdim = size(solve_time);
        real y_solve[numdim, 8];

        for (i in 1:numdim) {
          y_solve[i] = solve_chi(solve_time[i], ageAtBMT[i], init_cond, parms);
        }

        return y_solve;
      }

  vector math_reduce(vector global_params, vector local_params, real[] x_r, int[] x_i){
    // data for each shard
    int n = size(x_i); // n = 1
    real dat_time = x_r[1];
    int dat_t0 = x_i[1];                          // time zero -- for chimeras age at BMT
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
   int<lower = 1> numChi;
   real<lower = 0> chi_counts[numChi];
   real<lower = 0> N_donor_fraction[numChi];
   real<lower = 0> donor_ki[numChi];
   real<lower = 0> host_ki[numChi];
   int<lower  = 1> numPred;
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

parameters {
  real<lower= 0, upper= 1> psi;
  real<lower= 0, upper= 1> delta_nai;
  real<lower= 0, upper= 1> delta_rte;
  real<lower= 0, upper= 1> rho_nai;
  real<lower= 0, upper= 1> rho_rte;
  real<lower= 0, upper= 1> mu;
  real<lower=1E4, upper=5E6> N0;                  // total cells counts at t0
  real<lower= 0, upper= 1> kappa0;

  real<lower=0> sigma_chi_counts;
  real<lower=0> sigma_Nfd;
  real<lower=0> sigma_donor_ki;
  real<lower=0> sigma_host_ki;
}

transformed parameters{
  vector[8] global_params;
  vector[numChi] y3_mean;               // PDE prediction for counts from chimera data
  vector[numChi] y4_mean;               // PDE prediction for Nfd from chimera data
  vector[numChi] y5_mean;               // PDE prediction for ki proportions in donor compartment from chimera data
  vector[numChi] y6_mean;               // PDE prediction for ki proportions in host compartment from chimera data
  vector[4*numChi] y_mean_stacked;        // compliled output across all nodes

  global_params[1] = psi;
  global_params[2] = delta_nai;
  global_params[3] = delta_rte;
  global_params[4] = rho_nai;
  global_params[5] = rho_rte;
  global_params[6] = mu;
  global_params[7] = N0;
  global_params[8] = kappa0;

  // combining the output from all the shards
  y_mean_stacked = map_rect(math_reduce, global_params, local_params, x_r, x_i);

  for (i in 1:numChi){
    y3_mean[i] = y_mean_stacked[4*i - 3];
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
