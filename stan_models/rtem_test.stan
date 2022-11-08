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
}
