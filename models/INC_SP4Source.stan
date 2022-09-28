functions{
  real theta_spline(real time) {
      real theta0 = exp(15.29);   // theta0 gives the initial counts of the source compartment
      real nu = 0.00396;          //rate of decline of source influx assuming that its proportional to the counts of source compartment

      real theta =  theta0 * exp(-nu * time);
      return theta;
   }

   real Chi_spline(real time) {
       real chi;

       real chiEst = 0.82;       // chiEst is the level if stabilised chimerism in the source compartment
       real qEst = 0.047;        // qEst is the rate with which cimerism chnages in the source compartment

       if (time < 0){
         chi = 0;                       // conditioning the function to adapt to the timepoints before BMT
       } else {
         chi = chiEst * (1 - exp(-qEst * time));
       }
       return chi;
    }

  real[] inc(real time, real[] y, real[] parms, real[] rdata, int[] idata) {
     real dydt[5];             // system of ODEs

     real alpha1 = parms[1];                  // rate of influx into Thymic naive Tregs
     real alpha2 = parms[2];                  // rate of influx into Thymic naive Tregs
     real mu  = parms[3];                     // rate of migration from thymic naive Tregs to peripheral naive Tregs
     real beta  = parms[4];                   // rate of backcirculation into Thymic naive Tregs from peripheral naive Tregs
     real lambda = parms[5];                  // Net rate of loss of naive Tregs
     // age of BMT in each recipient
     real ageAtBMT = parms[6];

     real lambda_inc = 0.0;

     // thymic naive Tregs donor compartment
     dydt[1] = (alpha1 * theta_spline(time) * Chi_spline(time - ageAtBMT)) + beta * y[3] - (mu + lambda) * y[1];
     // thymic naive Tregs host compartment
     dydt[2] = (alpha1 * theta_spline(time) * (1 - Chi_spline(time - ageAtBMT))) + beta * y[4] - (mu + lambda) * y[2];

     // Peripheral naive Tregs donor compartment
     dydt[3] = (alpha2 * theta_spline(time) * Chi_spline(time - ageAtBMT)) + mu * y[1]  - (beta + lambda) * y[3];
     // Peripheral naive Tregs host compartment
     dydt[4] = (alpha2 * theta_spline(time) * (1 - Chi_spline(time - ageAtBMT))) + mu * y[2] - (beta + lambda) * y[4];

     // peripheral incumbent cells
     dydt[5] = - lambda_inc * y[5];

     return dydt;
 }

   // solving for total counts of thymic and peripheral naive Tregs at time of BMT assuming the youngest animal as the t0
   // these counts form the initial conditions for other recipients N(0) = N_host(0)
   real[,] solve_init(real[] tb_time,           // age at BMT
     real[] init_cond,                          // initial conditions at BMT in the youngest mouse
     real[] parms){

       int ndim = size(tb_time);
       int x_i[0];
       real y_init[ndim, 5];
       real params_init[6];

       params_init[1:5] = parms[1:5];
       params_init[6] = 40;                     // age at BMT for the youngest animal

       y_init[1] = init_cond;                 // init conditions at the earliest BMT (i.e. in younegst animal)
       y_init[2:ndim] = integrate_ode_rk45(inc, init_cond, tb_time[1], tb_time[2:ndim], params_init, {0.0}, x_i);

       return y_init;
   }

   real[] solve_ode(real time_point,            // time point of observation
     real ageAtBMT,
     real[] init_tBMT,                          // initial conditions at BMT in the youngest mouse
     real[] parms){

       real y_hat[2, 5];
       real init_tb[5];                         // init conditions at the time of BMT
       real params[6];

       params[1:5] = parms[1:5];
       params[6] = ageAtBMT;                                           // age at BMT

       // init conditions at the BMT
       init_tb[1] = 0;                                           //at tbmt - # donor thymic naive Treg is zero
       init_tb[2] = init_tBMT[1] + init_tBMT[2];                 //at tbmt - all thymic naive Treg cells are host
       init_tb[3] = 0;                                           //at tbmt - # donor Peripheral naive Treg is zero
       init_tb[4] = init_tBMT[3] + init_tBMT[4];                 //at tbmt - all Peripheral naive Treg cells are host
       init_tb[5] = init_tBMT[5];                 //at tbmt - all Peripheral naive Treg cells are host


       y_hat[1] = init_tb;
       y_hat[2] = to_array_1d(integrate_ode_rk45(inc, init_tb, ageAtBMT, rep_array(time_point, 1), params, {0.0}, {0}));

       return y_hat[2];
     }


  real[,] solve_timecourse(real[] data_time,            //timepoints in the data
      int[] ageAtBMT,                                  // age at BMT in data
      real[] tb_time,                                   // unique age at BMT to solve for using solve_init()
      int[] tb_index,                                  // indices of tb_time that match with age at BMT
      real[] init_cond,
      real[] parms){

        int num_obs = size(data_time);
        real y_solve[num_obs, 5];

        int num_tb = size(tb_time);
        real init_tBMT[num_tb, 5];

        //array of initial conditions at each unique age at BMT
        init_tBMT = solve_init(tb_time, init_cond, parms);

        // solving for counts in each mouse with its matching age at BMT and using initial conditions solved above
        for (i in 1:num_obs){
          y_solve[i] = solve_ode(data_time[i], ageAtBMT[i], init_tBMT[tb_index[i]], parms);
        }

        return y_solve;
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
  int<lower  = 1> num_obs;
  int<lower  = 1> num_tb;
  int<lower  = 1> num_index;
  real<lower = 0> data_time[num_obs];
  int<lower  = 0> ageAtBMT[num_obs];
  real<lower = 0> tb_time[num_tb];
  int<lower  = 0> tb_index[num_index];
  real<lower = 0> per_counts[num_obs];
  real<lower = 0> per_Nfd[num_obs];
  real<lower = 0> thy_counts[num_obs];
  real<lower = 0> thy_Nfd[num_obs];
  int<lower  = 1> num_pred;
  real ts_pred1[num_pred];
  real ts_pred2[num_pred];
  real ts_pred3[num_pred];
}

transformed data{
  real y1[num_obs];
  real y2[num_obs];
  real y3[num_obs];
  real y4[num_obs];

  y1 = log(thy_counts);                         // transforming cell counts of donor compartments
  y2 = asinsqrt_array(thy_Nfd);                 // untransfored cell counts of donor fractions normalised to source chimerism
  y3 = log(per_counts);                         // transforming cell counts of donor compartments
  y4 = asinsqrt_array(per_Nfd);                 // untransfored cell counts of donor fractions normalised to source chimerism
}

parameters{
  real<lower = 0> thy_N0Log;
  real<lower = 0> per_N0Log;
  real<lower = 0> per_I0Log;
  real<lower = 0, upper=1> alpha1;
  real<lower = 0, upper=1> alpha2;
  real<lower = 0, upper=1> mu;
  real<lower = 0, upper=1> beta;
  real lambda;

  real<lower = 0> sigma1;
  real<lower = 0> sigma2;
  real<lower = 0> sigma3;
  real<lower = 0> sigma4;
  }

transformed parameters{
  real y_hat[num_obs, 5];
  real y1_mean[num_obs];
  real y2_mean[num_obs];
  real y3_mean[num_obs];
  real y4_mean[num_obs];
  real parms[5];
  real init_cond[5];

  real thy_N0 = exp(thy_N0Log);
  real per_N0 = exp(per_N0Log);
  real per_I0 = exp(per_I0Log);

  init_cond[1] = 0;                                // count of donor cells in thymic naive T reg compartment at earliest ageAtBMT
  init_cond[2] = thy_N0 ;                          // count of host cells in thymic naive T reg compartment at earliest ageAtBMT
  init_cond[3] = 0;                                // count of donor cells in peripheral naive T reg compartment at earliest ageAtBMT
  init_cond[4] = per_N0 ;                          // count of host cells in thymic naive T reg compartment at earliest ageAtBMT
  init_cond[5] = per_I0 ;                          // count of host cells in thymic naive T reg compartment at earliest ageAtBMT
  parms[1] = alpha1;
  parms[2] = alpha2;
  parms[3] = mu;
  parms[4] = beta;
  parms[5] = lambda;

  // solution for ODEs for the whole timecourse
  y_hat = solve_timecourse(data_time, ageAtBMT, tb_time, tb_index, init_cond, parms);

  for (i in 1:num_obs){
    //Thymic Tregs
    // total counts
    y1_mean[i] = y_hat[i, 1] + y_hat[i, 2];

    // donor fractions normalised with chimerism in the source
    y2_mean[i] = y_hat[i, 1]/(y1_mean[i] * Chi_spline(data_time[i] - ageAtBMT[i]));

    //Peripheral Tregs
    // total counts
    y3_mean[i] = y_hat[i, 3] + y_hat[i, 4] + y_hat[i, 5];

    // donor fractions normalised with chimerism in the source
    y4_mean[i] = y_hat[i, 3]/(y3_mean[i] * Chi_spline(data_time[i] - ageAtBMT[i]));
  }
}

model{
  alpha1 ~ uniform(0.0, 1.0);
  alpha2 ~ uniform(0.0, 1.0);
  mu ~ normal(0.3, 0.2);
  beta ~ normal(0.1, 0.2);
  lambda ~ normal(0.05, 0.2);

  thy_N0Log ~ normal(10, 2);
  per_N0Log ~ normal(14, 2);
  per_I0Log ~ normal(10, 2);

  sigma1 ~ cauchy(0.1, 1);
  sigma2 ~ cauchy(0.1, 1);
  sigma3 ~ cauchy(0.1, 1);
  sigma4 ~ cauchy(0.1, 1);

  y1 ~ normal(log(y1_mean), sigma1);
  y2 ~ normal(asinsqrt_array(y2_mean), sigma2);
  y3 ~ normal(log(y3_mean), sigma3);
  y4 ~ normal(asinsqrt_array(y4_mean), sigma4);
}

generated quantities{
  real y_init_pred[3, 5];
  real y_solve_pred1[num_pred, 5]; real y_solve_pred2[num_pred, 5];  real y_solve_pred3[num_pred, 5];
  real y1_mean_pred1[num_pred];  real y1_mean_pred2[num_pred];  real y1_mean_pred3[num_pred];
  real y2_mean_pred1[num_pred];  real y2_mean_pred2[num_pred];  real y2_mean_pred3[num_pred];
  real y3_mean_pred1[num_pred];  real y3_mean_pred2[num_pred];  real y3_mean_pred3[num_pred];
  real y4_mean_pred1[num_pred];  real y4_mean_pred2[num_pred];  real y4_mean_pred3[num_pred];

  // log likelihoods for 4 peoperties i.e. counts of Thy and Per nai Tregs and theire respective NFDs.
  vector[num_obs] log_lik1;
  vector[num_obs] log_lik2;
  vector[num_obs] log_lik3;
  vector[num_obs] log_lik4;

  // init conditions for 3 prediction groups with mean age at BMT == 49, 72 and 128 viz.
  y_init_pred = solve_init({49.0, 72.0, 128.0}, init_cond, parms);

  // solutions of the ode model for different timecourses of 3 prediction groups
  for (i in 1:num_pred){
    y_solve_pred1[i] = solve_ode(ts_pred1[i], 49.0, y_init_pred[1], parms);
    y_solve_pred2[i] = solve_ode(ts_pred2[i], 72.0, y_init_pred[2], parms);
    y_solve_pred3[i] = solve_ode(ts_pred3[i], 128.0, y_init_pred[3], parms);
  }


  // Predictions for each property for the 3 prediction groups
  for (i in 1:num_pred){
    //Thymic Tregs
    // total counts
    y1_mean_pred1[i] = y_solve_pred1[i, 1] + y_solve_pred1[i, 2];
    y1_mean_pred2[i] = y_solve_pred2[i, 1] + y_solve_pred2[i, 2];
    y1_mean_pred3[i] = y_solve_pred3[i, 1] + y_solve_pred3[i, 2];

    // donor fractions normalised with chimerism in the source
    y2_mean_pred1[i] = y_solve_pred1[i, 1]/(y1_mean_pred1[i] * Chi_spline(ts_pred1[i] - 49.0));
    y2_mean_pred2[i] = y_solve_pred2[i, 1]/(y1_mean_pred2[i] * Chi_spline(ts_pred2[i] - 72.0));
    y2_mean_pred3[i] = y_solve_pred3[i, 1]/(y1_mean_pred3[i] * Chi_spline(ts_pred3[i] - 128.0));

    //Peripheral Tregs
    // total counts
    y3_mean_pred1[i] = y_solve_pred1[i, 3] + y_solve_pred1[i, 4] + y_solve_pred1[i, 5];
    y3_mean_pred2[i] = y_solve_pred2[i, 3] + y_solve_pred2[i, 4] + y_solve_pred1[i, 5];
    y3_mean_pred3[i] = y_solve_pred3[i, 3] + y_solve_pred3[i, 4] + y_solve_pred1[i, 5];

    // donor fractions normalised with chimerism in the source
    y4_mean_pred1[i] = y_solve_pred1[i, 3]/(y3_mean_pred1[i] * Chi_spline(ts_pred1[i] - 49.0));
    y4_mean_pred2[i] = y_solve_pred2[i, 3]/(y3_mean_pred2[i] * Chi_spline(ts_pred2[i] - 72.0));
    y4_mean_pred3[i] = y_solve_pred3[i, 3]/(y3_mean_pred3[i] * Chi_spline(ts_pred3[i] - 128.0));
  }

  // calculating the log predictive accuracy for each point
  for (n in 1:num_obs) {
    log_lik1[n] = normal_lpdf(y1[n] | log(y1_mean[n]), sigma1);
    log_lik2[n] = normal_lpdf(y2[n] | asinsqrt_real(y2_mean[n]), sigma2);
    log_lik3[n] = normal_lpdf(y3[n] | log(y3_mean[n]), sigma3);
    log_lik4[n] = normal_lpdf(y4[n] | asinsqrt_real(y4_mean[n]), sigma4);
  }
}
