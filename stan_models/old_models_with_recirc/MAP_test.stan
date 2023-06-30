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
  return exp(- eps_f * (time - t0)) + eps_0;
}

real[] shm_chi(real time, real[] y, real[] parms, real[] rdata,  int[] idata) {
  real psi = parms[1];
  real rho_D = parms[2];
  real alpha = parms[3];
  real delta_D = parms[4];
  real rho_I = parms[5];
  real beta = parms[6];
  real delta_I = parms[7];

  real dydt[12];
  real kloss  = 1/3.5;            //rate of loss of ki67
  real eps_host = 0.326611;      // Mean Ki67 hi fraction in host-BM-derived FoxP3 negative Sp4 T cells

  // age of BMT in each recipient
  real ageAtBMT = parms[8];

  // model that assumes that tranistionals divide and die at different rates than mature naive T cells
  // Host naive Tregs
  // Thymic ki  hi displaceable
  dydt[1] = theta_spline(time, psi) * (1- Chi_spline(time - ageAtBMT)) * eps_host + rho_D * (2 * y[2] + y[1]) - (kloss + alpha + delta_D) * y[1];
  // Thymic ki lo displaceable
  dydt[2] = theta_spline(time, psi) * (1- Chi_spline(time - ageAtBMT)) * (1 - eps_host) + kloss * y[1] - (rho_D + alpha + delta_D) * y[2];
  // Peripheral ki hi displaceable
  dydt[3] = alpha * y[1] + rho_D * (2 * y[4] + y[3]) - (kloss + beta + delta_D) * y[3];
  // Peripheral ki lo displaceable
  dydt[4] = alpha * y[2] + kloss * y[3] - (rho_D + beta + delta_D) * y[4];
  // Peripheral ki hi Incumbent
  dydt[5] = alpha * y[7] + rho_I * (2 * y[6] + y[5]) - (kloss + beta + delta_I) * y[5];
  // Peripheral ki lo Incumbent
  dydt[6] = alpha * y[8] + kloss * y[5] - (rho_I + beta + delta_I) * y[6];
  // Thymic ki hi Incumbent
  dydt[7] = beta * y[5] + rho_I * (2 * y[8] + y[7]) - (kloss + alpha + delta_I) * y[7];
  // Thymic ki lo Incumbent
  dydt[8] = beta * y[6] + kloss * y[7] - (rho_I + alpha + delta_I) * y[8];

  // Donor naive Tregs
  // Thymic ki  hi displaceable
  dydt[9] = theta_spline(time, psi) * Chi_spline(time - ageAtBMT) * donor_eps_spline(time) + rho_D * (2 * y[10] + y[9]) - (kloss + alpha + delta_D) * y[9];
  // Thymic ki lo displaceable
  dydt[10] = theta_spline(time, psi) * Chi_spline(time - ageAtBMT) * (1 - donor_eps_spline(time)) + kloss * y[9] - (rho_D + alpha + delta_D) * y[10];
  // Peripheral ki hi displaceable
  dydt[11] = alpha * y[9] + rho_D * (2 * y[12] + y[11]) - (kloss + beta + delta_D) * y[11];
  // Peripheral ki lo displaceable
  dydt[12] = alpha * y[10] + kloss * y[11] - (rho_D + beta + delta_D) * y[12];
  return dydt;
}

// solving for total counts of thymic and peripheral naive Tregs at time of BMT assuming the youngest animal as the t0
// these counts form the initial conditions for other recipients N(0) = N_host(0)
real[] solve_init(real ageAtBMT,
  real[] init_cond,                          // initial conditions at BMT in the youngest mouse
  real[] parms){

    real ta = 40;                            // age at BMT for the youngest host
    real y_init[2, 12];
    real params_init[8];

    params_init[1:7] = parms[1:7];
    params_init[8] = ta;

    y_init[1] = init_cond;                 // init conditions at the earliest BMT (i.e. in younegst animal)
    y_init[2] = to_array_1d(integrate_ode_rk45(shm_chi, init_cond, ta, rep_array(ageAtBMT, 1), params_init, {0.0}, {0}));

    return y_init[2];
}

real[] solve_chi(real solve_time, real ageAtBMT, real[] init_cond, real[] parms){
     real y_solve[12];
    real params[8];

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

    params[1:7] = parms[1:7];
    params[8] = ageAtBMT;                                           // age at BMT

    y_solve = to_array_1d(integrate_ode_rk45(shm_chi, init_tb, ageAtBMT, rep_array(solve_time, 1), params, {0.0}, {0}));

    return y_solve;
  }

 real[,] solve_ode_chi(real[] solve_time, real[] ageAtBMT, real[] init_cond, real[] parms){
    int numdim = size(solve_time);
    real y_solve[numdim, 12];
    for (i in 1:numdim) {
      y_solve[i] = solve_chi2(solve_time[i], ageAtBMT[i], init_cond, parms);
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
    real y1_0 = global_params[8]; real y2_0 = global_params[9];  real y3_0 = global_params[10];
    real y4_0 = global_params[11]; real y5_0 = global_params[12]; real y6_0 = global_params[13];
    real y7_0 = global_params[14]; real y8_0 = global_params[15];

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
    // ODE solution for chimera dataset -- x_r = data time and x_i = time at BMT
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

  }
