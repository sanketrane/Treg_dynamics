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
    real solve_time[n] = x_r[1:n];
    int ageAtBMT[n] = x_i[1:n];                          // time zero -- for chimeras age at BMT
    real tb_time[n];

    //params
    real y1_0 = global_params[9]; real y2_0 = global_params[10];  real y3_0 = global_params[11];
    real y4_0 = global_params[12]; real y5_0 = global_params[13]; real y6_0 = global_params[14];
    real y7_0 = global_params[15]; real y8_0 = global_params[16];

    real init_cond[16];

    // ODE solution -- predictions for the observed timecourse
    real chi_solve[n, 16];

    real counts_thy[n]; real counts_per[n]; real donor_counts_thy[n]; real donor_counts_per[n];
    real host_counts_thy[n]; real host_counts_per[n]; real donor_ki_thy[n]; real donor_ki_per[n];
    real host_ki_thy[n]; real host_ki_per[n];

    vector[8*n] y_mean_stacked;

    // ODE solution -- predictions for the observed timecourse
    init_cond[1] = y1_0; init_cond[2] = y2_0; init_cond[3] = y3_0; init_cond[4] = y4_0;
    init_cond[5] = y5_0; init_cond[6] = y6_0; init_cond[7] = y7_0; init_cond[8] = y8_0;
    init_cond[9] =  0; init_cond[10] = 0; init_cond[11] = 0; init_cond[12] = 0;
    init_cond[13] = 0; init_cond[14] = 0; init_cond[15] = 0; init_cond[16] = 0;

    for (i in 1:n){
      tb_time[i] = ageAtBMT[i]/1.0;
    }

    // each shard has a single datpoint so its unique ****
    // PDE solution for chimera dataset -- x_r = data time and x_i = time at BMT
      chi_solve = solve_ode_chi(solve_time, tb_time, init_cond, to_array_1d(global_params));

      for (i in 1:n){
        counts_thy[i] = chi_solve[i, 1] + chi_solve[i, 2] + chi_solve[i, 7] + chi_solve[i, 8] + chi_solve[i, 9] + chi_solve[i, 10] + chi_solve[i, 15] + chi_solve[i, 16];
        counts_per[i] = chi_solve[i, 3] + chi_solve[i, 4] + chi_solve[i, 5] + chi_solve[i, 6] + chi_solve[i, 11] + chi_solve[i, 12] + chi_solve[i, 13] + chi_solve[i, 14];
        donor_counts_thy[i] = chi_solve[i, 9] + chi_solve[i, 10] + chi_solve[i, 15] + chi_solve[i, 16];
        donor_counts_per[i] = chi_solve[i, 11] + chi_solve[i, 12] + chi_solve[i, 13] + chi_solve[i, 14];
        host_counts_thy[i] = chi_solve[i, 1] + chi_solve[i, 2] + chi_solve[i, 7] + chi_solve[i, 8];
        host_counts_per[i] = chi_solve[i, 3] + chi_solve[i, 4] + chi_solve[i, 5] + chi_solve[i, 6];
        donor_ki_thy[i] = (chi_solve[i, 9] + chi_solve[i, 15])/donor_counts_thy[i];
        donor_ki_per[i] = (chi_solve[i, 11] + chi_solve[i, 13])/donor_counts_thy[i];
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
