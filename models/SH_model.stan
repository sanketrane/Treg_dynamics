functions{
  real theta_spline(real t,
    real nu,               // nu is the rate of decline of source influx assuming that its proportional to the counts of source compartment
    real theta0) {         // theta0 gives the initial counts of the source compartment

      real theta;

      theta = theta0 * exp(-nu * t);
      return theta;
   }

   real Chi_spline( real t,
     real chiEst,         // chiEst is the level if stabilised chimerism in the source compartment
     real qEst) {         // qEst is the rate with which cimerism chnages in the source compartment

       real chi;

       if (t < 0){
         chi = 0;                       // conditioning the function to adapt to the timepoints before BMT
       } else {
         chi = chiEst * (1 - exp(-qEst * t));
       }
       return chi;
    }

  real[] shm(real t, real[] y, real[] parms, real[] rdata, int[] idata) {
     real dydt[4];             // system of ODEs

     real psi        = parms[1];
     real alpha      = parms[2];
     real lambda_thy = parms[3];
     real lambda_per = parms[4];

     real theta0 = rdata[1];
     real nu     = rdata[2];
     real chiEst = rdata[3];
     real qEst   = rdata[4];

     dydt[1] = (psi * theta_spline(t, nu, theta0) * Chi_spline(t, chiEst, qEst)) - lambda_thy * y[1];
     dydt[2] = (psi * theta_spline(t, nu, theta0) * (1 - Chi_spline(t, chiEst, qEst))) - lambda_thy * y[2];

     dydt[3] = alpha * y[1]  - lambda_per * y[3];
     dydt[4] = alpha * y[2]  - lambda_per * y[4];

     return dydt;
 }

 real[] foreach_init(real tb, real ta, real[] init_cond, real[] parms, real[] rdata, int x_i){

   return to_array_1d(integrate_ode_rk45(shm, init_cond, ta, rep_array(tb, 1), parms, rdata, rep_array(x_i, 1)));
 }

 real[,] solve_init(real[] tb_time, real[] init_cond, real[] parms, real[] rdata, int x_i, int num_tb){
  real y_init[num_tb, 4];

  y_init[1] = init_cond;
  for (i in 2:num_tb){
    y_init[i] = foreach_init(tb_time[i], tb_time[1], init_cond, parms, rdata, x_i);
  }

  return y_init;
 }

 real[] foreach_ode(real ts, real t0, real[] init_cond, real[] parms, real[] rdata, int x_i) {

   return to_array_1d(integrate_ode_rk45(shm, init_cond, t0, rep_array(ts, 1), parms, rdata, rep_array(x_i, 1)));
  }

 real[,] solve_ode(real[] solve_time, real[] init_cond, real[] parms, real[] rdata, int[] tb, int num_index, real[] tb_time, int[] tb_index, int num_tb){
 real y_hat[num_index, 4];
 real y0[num_tb, 4];
 real init_tb[4];

 y0 = solve_init(tb_time, init_cond, parms, rdata, 40, num_tb);

 init_tb[1] = init_cond[1];                                           //at tbmt - donor ki67Hi subset size is zero
 init_tb[2] = init_cond[2];                                           //at tbmt - donor ki67Lo subset size is zero

 for (i in 1:num_index){
   init_tb[3] = y0[tb_index[i], 1] + y0[tb_index[i], 3];              //at tbmt - all ki67Hi cells would be host
   init_tb[4] = y0[tb_index[i], 2] + y0[tb_index[i], 4];              //at tbmt - all ki67Lo cells would be host
   y_hat[i] = foreach_ode(solve_time[i], tb[i], init_tb, parms, rdata, tb[i]);
 }

 return y_hat;
 }

 real[,] solve_ode_pred(real[] solve_time, real[] init_cond, real[] parms, real[] rdata, int[] tb, int num_index, real[] tb_time){
 real y_hat[num_index, 4];
 real y0[2, 4];
 real init_tb[4];

 y0 = solve_init(tb_time, init_cond, parms, rdata, 40, 2);

 init_tb[1] = init_cond[1];                                         //at tbmt - donor ki67Lo subset size is zero
 init_tb[2] = init_cond[1];                                         //at tbmt - donor ki67Lo subset size is zero

 init_tb[3] = y0[2, 1] + y0[2, 3];                                  //at tbmt - all ki67Hi cells would be host
 init_tb[4] = y0[2, 2] + y0[2, 4];                                  //at tbmt - all ki67Hi cells would be host

 y_hat[1] = init_tb;
 for (i in 2:num_index){
   y_hat[i] = foreach_ode(solve_time[i], tb[i], init_tb, parms, rdata, tb[i]);
 }

 return y_hat;
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
  int<lower  = 1> numObs;
  int<lower  = 1> num_index;
  real<lower = 0> solve_time[num_index];
  int<lower  = 0> time_index[numObs];
  int<lower  = 1> num_tb;
  real<lower = 0> tb_time[num_tb];
  int<lower  = 0> tb_index[num_index];
  int<lower  = 0> ageBMT[num_index];
  int dpBMT[numObs];
  real<lower = 0> per_counts[numObs];
  real<lower = 0> per_Nfd[numObs];
  real<lower = 0> thy_counts[numObs];
  real<lower = 0> thy_Nfd[numObs];
  int<lower  = 1> numPred1;
  int<lower  = 1> numPred2;
  int<lower  = 1> numPred3;
  real ts_pred1[numPred1];
  real ts_pred2[numPred2];
  real ts_pred3[numPred3];
  int tb_pred1[numPred1];
  int tb_pred2[numPred2];
  int tb_pred3[numPred3];
  real tb_time_pred1[2];
  real tb_time_pred2[2];
  real tb_time_pred3[2];
  real thy_Nd0;
  real per_Nd0;
  real theta0;
  real nu;
  real chiEst;
  real qEst;
  }

transformed data{
  real y1[numObs];
  real y2[numObs];
  real y3[numObs];
  real y4[numObs];
  real rdata[4];
  int tb[num_index];

  rdata[1] = theta0;
  rdata[2] = nu;
  rdata[3] = chiEst;
  rdata[4] = qEst;

  y1 = log(thy_counts);                 // transforming cell counts of donor compartments to feed in to ODEs
  y2 = asinsqrt_array(thy_Nfd);                       // untransfored cell counts of donor fractions normalised to source chimerism to feed in to ODEs
  y3 = log(per_counts);                 // transforming cell counts of donor compartments to feed in to ODEs
  y4 = asinsqrt_array(per_Nfd);                       // untransfored cell counts of donor fractions normalised to source chimerism to feed in to ODEs
}

parameters{
  real thy_y0Log;
  real per_y0Log;
  real<lower = 0, upper=1> psi;
  real<lower = 0, upper=1> alpha;
  real lambda_thy;
  real lambda_per;

  real<lower = 0> sigma1;
  real<lower = 0> sigma2;
  real<lower = 0> sigma3;
  real<lower = 0> sigma4;
  }

transformed parameters{
  real y_hat[num_index, 4];
  real y1_mean[numObs];
  real y2_mean[numObs];
  real y3_mean[numObs];
  real y4_mean[numObs];
  real parms[4];
  real init_cond[4];
  real thy_y0;
  real per_y0;

  thy_y0 = exp(thy_y0Log);
  per_y0 = exp(per_y0Log);

  init_cond[1] = thy_Nd0;
  init_cond[2] = thy_y0 ;
  init_cond[3] = per_Nd0;
  init_cond[4] = per_y0 ;
  parms[1] = psi;
  parms[2] = alpha;
  parms[3] = lambda_thy;
  parms[4] = lambda_per;

  // solution for ODEs
  y_hat = solve_ode(solve_time, init_cond, parms, rdata, tb, num_index, tb_time, tb_index, num_tb);

  for (i in 1:numObs){
    //Thymic Tregs
    // total counts
    y1_mean[i] = y_hat[time_index[i], 1] + y_hat[time_index[i], 2];

    // donor fractions normalised with chimerism in the source
    y2_mean[i] = (y_hat[time_index[i], 1])/(y1_mean[i] * Chi_spline(dpBMT[i], chiEst, qEst));

    //Peripheral Tregs
    // total counts
    y3_mean[i] = y_hat[time_index[i], 3] + y_hat[time_index[i], 4];

    // donor fractions normalised with chimerism in the source
    y4_mean[i] = (y_hat[time_index[i], 3])/(y3_mean[i] * Chi_spline(dpBMT[i], chiEst, qEst));
  }
}

model{
  psi ~ normal(0.01, 0.5);
  alpha ~ normal(0.01, 0.5);
  thy_y0Log ~ normal(13, 2);
  per_y0Log ~ normal(11, 2);
  lambda_thy ~ normal(0.01, 2);
  lambda_per ~ normal(0.01, 2);

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
  real y_hat_pred_age1[numPred1, 4];
  real y_hat_pred_age2[numPred2, 4];
  real y_hat_pred_age3[numPred3, 4];
  real y1_mean_pred_age1[numPred1];
  real y1_mean_pred_age2[numPred2];
  real y1_mean_pred_age3[numPred3];
  real y2_mean_pred_age1[numPred1];
  real y2_mean_pred_age2[numPred2];
  real y2_mean_pred_age3[numPred3];
  real y3_mean_pred_age1[numPred1];
  real y3_mean_pred_age2[numPred2];
  real y3_mean_pred_age3[numPred3];
  real y4_mean_pred_age1[numPred1];
  real y4_mean_pred_age2[numPred2];
  real y4_mean_pred_age3[numPred3];
  real thy_countspred_age1[numPred1];
  real thy_countspred_age2[numPred2];
  real thy_countspred_age3[numPred3];
  real thy_fdpred_age1[numPred1];
  real thy_fdpred_age2[numPred2];
  real thy_fdpred_age3[numPred3];
  real per_countspred_age1[numPred1];
  real per_countspred_age2[numPred2];
  real per_countspred_age3[numPred3];
  real per_fdpred_age1[numPred1];
  real per_fdpred_age2[numPred2];
  real per_fdpred_age3[numPred3];
  vector[numObs] log_lik;
  vector[numObs] log_lik1;
  vector[numObs] log_lik2;
  vector[numObs] log_lik3;
  vector[numObs] log_lik4;

  //parameters of interest
  real lambdaThy_inv = 1/lambda_thy;
  real lambdaPer_inv = 1/lambda_per;


  //ODE solution for different age bins
  y_hat_pred_age1 = solve_ode_pred(ts_pred1, init_cond, parms, rdata, tb_pred1, numPred1, tb_time_pred1);
  y_hat_pred_age2 = solve_ode_pred(ts_pred2, init_cond, parms, rdata, tb_pred2, numPred2, tb_time_pred2);
  y_hat_pred_age3 = solve_ode_pred(ts_pred3, init_cond, parms, rdata, tb_pred3, numPred3, tb_time_pred3);

  // Initial conditions for total counts anf fd of thymic and peripheral Tregs
  for (i in 1:numPred1){
    y1_mean_pred_age1[i] = y_hat_pred_age1[i, 1] + y_hat_pred_age1[i, 2];
    thy_countspred_age1[i] = exp(normal_rng(log(y1_mean_pred_age1[i]), sigma1));

    y1_mean_pred_age2[i] = y_hat_pred_age2[i, 1] + y_hat_pred_age2[i, 2];
    thy_countspred_age2[i] = exp(normal_rng(log(y2_mean_pred_age2[i]), sigma1));

    y1_mean_pred_age3[i] = y_hat_pred_age3[i, 1] + y_hat_pred_age3[i, 2];
    thy_countspred_age3[i] = exp(normal_rng(log(y3_mean_pred_age3[i]), sigma1));

    y3_mean_pred_age1[i] = y_hat_pred_age1[i, 3] + y_hat_pred_age1[i, 4];
    per_countspred_age1[i] = exp(normal_rng(log(y3_mean_pred_age1[i]), sigma3));

    y3_mean_pred_age2[i] = y_hat_pred_age2[i, 3] + y_hat_pred_age2[i, 4];
    per_countspred_age2[i] = exp(normal_rng(log(y3_mean_pred_age2[i]), sigma3));

    y3_mean_pred_age3[i] = y_hat_pred_age3[i, 3] + y_hat_pred_age3[i, 4];
    per_countspred_age3[i] = exp(normal_rng(log(y3_mean_pred_age3[i]), sigma3));
  }

  y2_mean_pred_age1[1] = thy_Nd0;
  thy_fdpred_age1[1]   = thy_Nd0;
  y2_mean_pred_age2[1] = thy_Nd0;
  thy_fdpred_age2[1]   = thy_Nd0;
  y2_mean_pred_age3[1] = thy_Nd0;
  thy_fdpred_age3[1]   = thy_Nd0;

  y4_mean_pred_age1[1] = per_Nd0;
  per_fdpred_age1[1]   = per_Nd0;
  y4_mean_pred_age2[1] = per_Nd0;
  per_fdpred_age2[1]   = per_Nd0;
  y4_mean_pred_age3[1] = per_Nd0;
  per_fdpred_age3[1]   = per_Nd0;


  for (i in 1:numPred1){
    y2_mean_pred_age1[i] = y_hat_pred_age1[i, 1] / (y1_mean_pred_age1[i] * Chi_spline(ts_pred1[i], chiEst, qEst));
    thy_fdpred_age1[i] = asinsqrt_inv(normal_rng(asinsqrt_real(y2_mean_pred_age1[i]), sigma2));

    y2_mean_pred_age2[i] = y_hat_pred_age2[i, 1] / (y1_mean_pred_age2[i] * Chi_spline(ts_pred2[i], chiEst, qEst));
    thy_fdpred_age2[i] = asinsqrt_inv(normal_rng(asinsqrt_real(y2_mean_pred_age2[i]), sigma2));

    y2_mean_pred_age3[i] = y_hat_pred_age3[i, 1] / (y1_mean_pred_age3[i] * Chi_spline(ts_pred3[i], chiEst, qEst));
    thy_fdpred_age3[i] = asinsqrt_inv(normal_rng(asinsqrt_real(y2_mean_pred_age3[i]), sigma2));


    y4_mean_pred_age1[i] = y_hat_pred_age1[i, 3] / (y3_mean_pred_age1[i] * Chi_spline(ts_pred1[i], chiEst, qEst));
    per_fdpred_age1[i] = asinsqrt_inv(normal_rng(asinsqrt_real(y2_mean_pred_age1[i]), sigma4));

    y4_mean_pred_age2[i] = y_hat_pred_age2[i, 3] / (y3_mean_pred_age2[i] * Chi_spline(ts_pred2[i], chiEst, qEst));
    per_fdpred_age2[i] = asinsqrt_inv(normal_rng(asinsqrt_real(y2_mean_pred_age2[i]), sigma4));

    y4_mean_pred_age3[i] = y_hat_pred_age3[i, 3] / (y3_mean_pred_age3[i] * Chi_spline(ts_pred3[i], chiEst, qEst));
    per_fdpred_age3[i] = asinsqrt_inv(normal_rng(asinsqrt_real(y2_mean_pred_age3[i]), sigma4));
    }

  // calculating the log predictive accuracy for each point
  for (n in 1:numObs) {
    log_lik1[n] = normal_lpdf(y1[n] | log(y1_mean[n]), sigma1);
    log_lik2[n] = normal_lpdf(y2[n] | asinsqrt_real(y2_mean[n]), sigma2);
    log_lik3[n] = normal_lpdf(y3[n] | log(y3_mean[n]), sigma3);
    log_lik4[n] = normal_lpdf(y4[n] | asinsqrt_real(y4_mean[n]), sigma4);
    log_lik[n] = log_lik1[n] + log_lik2[n];
  }
}
