functions{
   //spline 1
 // Timecourse of thymic CD4 SP population -- changes with time
 real sp_numbers(real time) {
   real t0 = 40.0;     // mean mouse age at BMT for the first ageBMT bin
   real value;
   // spline fitted separately to the counts of thymic total SP4 T cells to estimate the parameters
   real b0  = 5.94; real b1 = 6.52;  real nu = 91 ;
   //best fitting spline
   value = 10^b0 + (10^b1/(1 + ((time - t0)/nu)^2));
   return value;
 }
 
  // initial age distribution of cells exisiting at t0
  real g_age(real age, real[] parms) {
    real t0 = 40.0;
    real N0 = parms[1];                         // cell counts at t0
    real p_age = parms[2];                      // parameter controling age distribution at t0
    real psi = (N0 * p_age)/(sp_numbers(t0) * (1-exp(-p_age * t0)));
    real theta_0 = psi * sp_numbers(t0);       // thymic output at t0 i.e. number of most recent emigrants (cells of age 0)
    real value;

    // Flat, normalised age-distribution of initial cells
    if(age >= 0 && age <= t0) {
        value = theta_0 * exp(-p_age * age);
        } else {
            value = 0;
            }
    return value;
  }

  // Total influx into the naive T cell compartment from the thymus (cells/day)
  real theta_spline(real time, real[] parms){
    real t0 = 40.0;
    real N0 = parms[1];                         // cell counts at t0
    real p_age = parms[2];                      // parameter controling age distribution at t0
    real psi = (N0 * p_age)/(sp_numbers(t0) * (1-exp(-p_age * t0)));
    
    return psi * sp_numbers(time);
  }
  
  // spline2 --
  real Chi_spline(real time) {
    // chiEst is the level of stabilised chimerism in the source (FoxP3 negative SP4) compartment
    // qEst is the rate with which cimerism changes in the source (FoxP3 negative SP4) compartment
    // spline fitted separately to the donor chimerism in the thymic FoxP3 negative SP4 T cells to estimate the parameters
    real value;  real chiEst = 0.8;   real qEst = 0.1;
    if (time - 10 < 0){              // t0 = 14 -- assumption: no donor cells seen in FoxP3neg SP4 compartment for 2 weeks
      value = 0;
      } else {
        value = chiEst * (1 - exp(-qEst * (time - 10)));
        }
    return value;
  }

  // influx of donor cells into the naive donor T cell compartment from the thymus (cells/day)
  real theta_donor(real time, real[] parms){
    real value;
    real tBMT = parms[6];
    //value = theta_spline(time, parms) * (time/tC)^m;
    value = theta_spline(time, parms) * Chi_spline(time - tBMT);
    return value;
  }

  // influx of host derived cells into the naive host T cell compartment from the thymus (cells/day)
  real theta_host(real time, real[] parms){
    return theta_spline(time, parms) - theta_donor(time, parms);
  }

// spline3 --
// proportions of ki67hi cells in the donor-derived FoxP3 negative SP4 T cells -- varies with time
real donor_eps_spline(real time){
  real t0 = 66.0;     // mean mouse age at BMT for the first ageBMT bin
  //parameters estimated from spline fit to the timecourse of ki67 fraction in the donor-derived FoxP3 negative SP4 T cells
  real b0 = 0.37; real b1= 0.27; real eps_f = 16.2;
  return b0 + (b1/(1 + (time/eps_f)^2));
}

// Ki67 distribution within the thymic influx -- varies with time
real ki_dist_theta(real ki, real time, int subpop){
  real k_bar = 1/exp(1.0);
  real value;
  real t0 = 1.0;
  
 if(subpop == 1){  // subpop = 1 is total thymic nai Treg pool
    if(ki >= 0 && ki < k_bar){
      value = (1 - 0.13)/k_bar;
      } else if (ki >= k_bar && ki <= 1.0){
          value = (0.13/(1 - k_bar));
          } else {
              value = 0.0;
              }           
  } else if (subpop == 2) {  // subpop = 1 is host thymic nai Treg pool
     if(ki >= 0 && ki < k_bar){
      value = (1 - 0.16)/k_bar;
      } else if (ki >= k_bar && ki <= 1.0){
          value = (0.16/(1 - k_bar));
          } else {
              value = 0.0;
              }               
   } else {  // subpop != 1 is donor thymic nai Treg pool
      if(ki >= 0 && ki < k_bar){
          value = (1 - donor_eps_spline(time))/k_bar;
          } else if (ki >= k_bar && ki <= 1.0){
              value = (donor_eps_spline(time)/(1 - k_bar));
              } else {
                  value = 0.0;
                  }     
   } 
  return value;
}

// Ki67 distribution of cells exisiting in the periphery at t0
real ki_dist_init(real ki){
  real value;
  real r_ki_init  = 5.69;         // parameter that shapes ki distribution within the init cohort
  
  if(ki >= 0.0 && ki <= 1.0){
    value = exp(-ki * r_ki_init)/((1 - exp(-r_ki_init))/r_ki_init);
    }  else {
      value = 0.0;
      }
  return value;
}

  // rate of cell division depending on cell age
  real rho_age(real age, real[] parms){
    real rho   = parms[4];
    real r_rho = parms[5];

    real value  = rho * exp(-r_rho * age);

    return value;
  }

  // function that calculates intgral of net loss rate --  solved analytically to speed up the numerical integration
  real lambda_integ(real lo_lim, real up_lim, real[] parms){
    real delta = parms[3];
    real rho   = parms[4];
    real r_rho = parms[5];

    real value = (delta * (up_lim - lo_lim)) + ((rho/r_rho) * (exp(-r_rho * up_lim) - exp(-r_rho * lo_lim)));
    return value;
  }

  // function that calculates intgral of net process rate --  solved analytically to speed up the numerical integration
  real alpha_integ(real lo_lim, real up_lim, real[] parms){
    real delta = parms[3];
    real rho   = parms[4];
    real r_rho = parms[5];

    real value = (delta * (up_lim - lo_lim)) + ((rho/r_rho) * (exp(-r_rho * lo_lim) - exp(-r_rho * up_lim)));
    return value;
  }

  
  // Cell age distribution of the initial cohort
  real Asm_init_age(real age, real time, real[] parms) {
    real t0 = 40.0;

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
    real t0 = 40.0;

    if(age < (time - t0)) {
      value =  theta_spline(time - age, parms) * exp(- lambda_integ(0, age, parms));
      } else {
        value = g_age(t0, parms) * exp(- lambda_integ(age - time + t0, age, parms));
    }

    return value;
  }

  // Cell age distribution of the initial cohort
  real Asm_Host_init_age(real age, real time, real[] parms) {
    real tBMT = parms[6];
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
     real tBMT = parms[6];
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
     real tBMT = parms[6];
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
     real tBMT = parms[6];
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
     real tBMT = parms[6];
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
     real tBMT = parms[6];
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
    real time = parms[6];   // time (host age) input  as a param

    value = Asm_total_age(age, time, parms);

    return {value};
  }

  real[] Asm_pooled_ode(real age,  real[] y, real[] parms, real[] x_r,  int[] x_i) {
    real value;
    real time = parms[7];   // time (host age) input  as a param

    value = Asm_pooled_age(age, time, parms);

    return {value};
  }

  real[] Asm_host_ode(real age,  real[] y, real[] parms, real[] x_r,  int[] x_i) {
    real value;
    real time = parms[7];   // time (host age) input  as a param

    value = Asm_host_age(age, time, parms);

    return {value};
  }

  real[] Asm_donor_ode(real age,  real[] y, real[] parms, real[] x_r,  int[] x_i) {
    real value;
    real time = parms[7];   // time (host age) input  as a param

    value = Asm_donor_age(age, time, parms);

    return {value};
  }

  real solve_total_counts(real[] parms) {
    int x_i[0];
    real value;
    real time = parms[6];   // time (host age) input  as a param

    // integrate_ode_rk45(function, y0, t0, t, theta, x_r, x_i);
    value = integrate_ode_rk45(Asm_total_ode, {0.0}, 0.0, rep_array(time, 1), parms, {0.0}, x_i)[1, 1];
    return value;
  }

  real solve_pooled_counts(real[] parms) {
    int x_i[0];
    real value;
    real time = parms[7];   // time (host age) input  as a param

    value = integrate_ode_rk45(Asm_pooled_ode, {0.0}, 0.0, rep_array(time, 1), parms, {0.0}, x_i)[1, 1];
    return value;
  }

  real solve_host_counts(real[] parms) {
    int x_i[0];
    real value;
    real time = parms[7];   // time (host age) input  as a param

    value = integrate_ode_rk45(Asm_host_ode, {0.0}, 0.0, rep_array(time, 1), parms, {0.0}, x_i)[1, 1];
    return value;
  }

  real solve_donor_counts(real[] parms) {
    int x_i[0];
    real value;
    real time = parms[7];   // time (host age) input  as a param

    value = integrate_ode_rk45(Asm_donor_ode, {0.0}, 0.0, rep_array(time, 1), parms, {0.0}, x_i)[1, 1];
    return value;
  }

  // Vectorised function for total pool size
  real[] N_total_time(real[] time, real[] parms){
   int ndim = size(time);
   real y_solve[ndim];
   real params[6];
   params[1:5] = parms[1:5];

   for (i in 1:ndim){
     params[6] = time[i];
     y_solve[i] = solve_total_counts(params);
   }
   return y_solve;
  }

  // Vectorised function for total pool size
  real[] N_pooled_time(real[] time, real[] tBMT, real[] parms){
    int ndim = size(time);
    real y_solve[ndim];
    real params[7];
    params[1:5] = parms[1:5];

    for (i in 1:ndim){
      params[6] = tBMT[i];
      params[7] = time[i];
      y_solve[i] = solve_pooled_counts(params);
    }
    return y_solve;
  }

  // Vectorised function for total pool size
  real[] N_host_time(real[] time, real[] tBMT, real[] parms){
    int ndim = size(time);
    real y_solve[ndim];
    real params[7];
    params[1:5] = parms[1:5];

    for (i in 1:ndim){
      params[6] = tBMT[i];
      params[7] = time[i];
      y_solve[i] = solve_host_counts(params);
    }
    return y_solve;
  }

  // Vectorised function for total pool size
  real[] N_donor_time(real[] time, real[] tBMT, real[] parms){
    int ndim = size(time);
    real y_solve[ndim];
    real params[7];
    params[1:5] = parms[1:5];

    for (i in 1:ndim){
      params[6] = tBMT[i];
      params[7] = time[i];
      y_solve[i] = solve_donor_counts(params);
    }
    return y_solve;
  }

  // init cohort -- The distribution of 'ki' for ages 'age' at times 'time'
  real U_init_ki_age(real ki, real age, real time, real[] parms){
    real t0 = 40.0;
    real beta  = 1/3.5;             // rate of loss of ki67 expression.
    real tau = -log(ki)/beta;     // time since cell divisions
    real value;

    if (ki <= exp(-beta * (time - t0)) ){
      value = g_age(t0, parms) * ki_dist_init(ki * exp(beta * (time - t0))) * exp(beta * (time - t0)) * exp(- alpha_integ(age - time + t0, age, parms));
    } else {
      value = 2.0 * rho_age(age, parms) * Asm_init_age(age - tau, time - tau, parms) * (1/(beta * ki)) * exp(- alpha_integ(age - tau, age, parms));
    }
    return value;
  }

  // theta cohort -- The distribution of 'ki' for ages 'age' at times 'time'
  real U_theta_ki_age(real ki, real age, real time, real[] parms){
    real t0 = 40.0;
    real beta  = 1/3.5;             // rate of loss of ki67 expression.
    real tau = -log(ki)/beta;     // time since cell divisions
    real value;

    if (ki <= exp(-beta * age)) {
      value = theta_spline(time - age, parms) * ki_dist_theta(ki * exp(beta * age), time-age, 1) * exp(beta * age) * exp(- alpha_integ(0.0, age, parms));
    } else {
      value = 2.0 * rho_age(age, parms) * Asm_theta_age(age - tau, time - tau, parms) * (1/(beta * ki)) * exp(- alpha_integ(age - tau, age, parms));
    }
    return value;
  }

  // theta cohort -- The distribution of 'ki' for ages 'age' at times 'time'
  real U_total_ki_age(real ki, real age, real time, real[] parms){
    real t0 = 40.0;
    real value;

    if (age < time - t0){
      value = U_theta_ki_age(ki, age, time, parms);
    } else {
      value = U_init_ki_age(ki, age, time, parms);
    }
    return value;
  }

  // init cohort -- The distribution of 'ki' for ages 'age' at times 'time'
  real host_init_ki_age(real ki, real age, real time, real[] parms){
    real tBMT = parms[6];
    real value;
    real beta  = 1/3.5;             // rate of loss of ki67 expression.
    real tau = -log(ki)/beta;     // time since cell divisions

    if (ki <= exp(-beta * (time - tBMT)) ){
      value = U_total_ki_age(ki * exp(beta * (time - tBMT)), age - time + tBMT, time, parms) * exp(beta * (time - tBMT)) * exp(- alpha_integ(age - time + tBMT, age, parms));
    } else {
      value = 2.0 * rho_age(age, parms) * Asm_Host_init_age(age - tau, time - tau, parms) * (1/(beta * ki)) * exp(- alpha_integ(age - tau, age, parms));
    }
    return value;
  }

  // theta cohort -- The distribution of 'ki' for ages 'age' at times 'time'
  real host_theta_ki_age(real ki, real age, real time, real[] parms){
    real t0 = 40.0;
    real beta  = 1/3.5;             // rate of loss of ki67 expression.
    real tau = -log(ki)/beta;     // time since cell divisions
    real value;

    if (ki <= exp(-beta * age)) {
      value = theta_host(time - age, parms) * ki_dist_theta(ki * exp(beta * age), time-age, 2) * exp(beta * age) * exp(- alpha_integ(0.0, age, parms));
    } else {
      value = 2.0 * rho_age(age, parms) * Asm_Host_theta_age(age - tau, time - tau, parms) * (1/(beta * ki)) * exp(- alpha_integ(age - tau, age, parms));
    }
    return value;
  }

  // init cohort -- The distribution of 'ki' for ages 'age' at times 'time'
  real donor_theta_ki_age(real ki, real age, real time, real[] parms){
    real t0 = 40.0;
    real beta  = 1/3.5;             // rate of loss of ki67 expression.
    real tau = -log(ki)/beta;       // time since cell divisions
    real value;

    if (ki <= exp(-beta * age)) {
      value = theta_donor(time - age, parms) * ki_dist_theta(ki * exp(beta * age), time-age, 3) * exp(beta * age) * exp(- alpha_integ(0.0, age, parms));
    } else {
      value = 2.0 * rho_age(age, parms) * Asm_Donor_theta_age(age - tau, time - tau, parms) * (1/(beta * ki)) * exp(- alpha_integ(age - tau, age, parms));
    }

    return value;
  }

  // The integrand function for age distribution of cells of ki intenisty 'ki'
  real[] U_total_kat(real ki, real[] y, real[] parms, real[] x_r, int[] x_i){
    real time = parms[6];
    real age = parms[7];
    real t0 = 40.0;
    real value;

    if (age < time - t0){
      value = U_theta_ki_age(ki, age, time, parms);
    } else {
      value = U_init_ki_age(ki, age, time, parms);
    }
    return {value};
  }

  // The integrand function for age distribution of cells of ki intenisty 'ki' for the pooled donor and host compartments
  real[] U_Pooled_kat(real ki, real[] y, real[] parms, real[] x_r, int[] x_i){
    real time = parms[7];
    real age = parms[8];
    real tBMT = parms[6];
    real value;

    if (age < time - tBMT){
      value = donor_theta_ki_age(ki, age, time, parms) + host_theta_ki_age(ki, age, time, parms);
    } else {
      value = host_init_ki_age(ki, age, time, parms);
    }
    return {value};
  }

  // The integrand function for age distribution of cells of ki intenisty 'ki' for the pooled donor and host compartments
  real[] U_host_kat(real ki, real[] y, real[] parms, real[] x_r, int[] x_i){
    real time = parms[7];
    real age = parms[8];
    real tBMT = parms[6];
    real value;

    if (age < time - tBMT){
      value = host_theta_ki_age(ki, age, time, parms);
    } else {
      value = host_init_ki_age(ki, age, time, parms);
    }
    return {value};
  }

  // The integrand function for age distribution of cells of ki intenisty 'ki' for the pooled donor and host compartments
  real[] U_donor_kat(real ki, real[] y, real[] parms, real[] x_r, int[] x_i){
    real time = parms[7];
    real age = parms[8];
    real tBMT = parms[6];
    real value;

    if (age < time - tBMT){
      value = donor_theta_ki_age(ki, age, time, parms);
    } else {
      value = 0.0;
    }
    return {value};
  }

  // integral across ki values -- kbar to 1.0
  real U_total_at(real[] parms){
    int x_i[0];
    real k_bar = 1/exp(1);

    real y_solve = integrate_ode_rk45(U_total_kat, {0.0}, 1/exp(1), rep_array(1.0, 1), parms, {0.0}, x_i)[1, 1];
    return y_solve;
  }

  // integral across ki values -- kbar to 1.0 for the pooled donor and host compartments
  real U_Pooled_at(real[] parms){
    int x_i[0];
    real k_bar = 1/exp(1);

    real y_solve = integrate_ode_rk45(U_Pooled_kat, {0.0}, 1/exp(1), rep_array(1.0, 1), parms, {0.0}, x_i)[1, 1];
    return y_solve;
  }

  // integral across ki values -- kbar to 1.0 for the pooled donor and host compartments
  real U_host_at(real[] parms){
    int x_i[0];
    real k_bar = 1/exp(1);

    real y_solve = integrate_ode_rk45(U_host_kat, {0.0}, 1/exp(1), rep_array(1.0, 1), parms, {0.0}, x_i)[1, 1];
    return y_solve;
  }

  // integral across ki values -- kbar to 1.0 for the pooled donor and host compartments
  real U_donor_at(real[] parms){
    int x_i[0];
    real k_bar = 1/exp(1);

    real y_solve = integrate_ode_rk45(U_donor_kat, {0.0}, 1/exp(1), rep_array(1.0, 1), parms, {0.0}, x_i)[1, 1];
    return y_solve;
  }

  // the age distribution
  real U_total_age(real age, real[] parms){
    real params[7];
    real value;
    params[1:6] = parms[1:6];
    params[7] = age;

    return U_total_at(params);
  }

  // the age distribution
  real U_Pooled_age(real age, real[] parms){
    real params[8];
    real value;
    params[1:7] = parms[1:7];
    params[8] = age;

    return U_Pooled_at(params);
  }

  // the age distribution
  real U_host_age(real age, real[] parms){
    real params[8];
    real value;
    params[1:7] = parms[1:7];
    params[8] = age;

    return U_host_at(params);
  }

  // the age distribution
  real U_donor_age(real age, real[] parms){
    real params[8];
    real value;
    params[1:7] = parms[1:7];
    params[8] = age;

    return U_donor_at(params);
  }

  // The integrand function for age distribution at time 'time'
  real[] U_total_ode(real age, real[] y, real[] parms, real[] x_r, int[] x_i){
    real time = parms[6];
    real value = U_total_age(age, parms);

    return {value};
  }

  // The integrand function for age distribution at time 'time'
  real[] U_Pooled_ode(real age, real[] y, real[] parms, real[] x_r, int[] x_i){
    real time = parms[7];
    real value = U_Pooled_age(age, parms);

    return {value};
  }

  // The integrand function for age distribution at time 'time'
  real[] U_host_ode(real age, real[] y, real[] parms, real[] x_r, int[] x_i){
    real time = parms[7];
    real value = U_host_age(age, parms);

    return {value};
  }

  // The integrand function for age distribution at time 'time'
  real[] U_donor_ode(real age, real[] y, real[] parms, real[] x_r, int[] x_i){
    real time = parms[7];
    real value = U_donor_age(age, parms);

    return {value};
  }

  // integral across age values -- 0.0 to time
  real U_total_t(real[] parms){
    int x_i[0];
    real time = parms[6];

    real y_solve = integrate_ode_rk45(U_total_ode, {0.0}, 0.0, rep_array(time, 1), parms, {0.0}, x_i)[1, 1];
    return y_solve;
  }

  // integral across age values -- 0.0 to time
  real U_Pooled_t(real[] parms){
    int x_i[0];
    real time = parms[7];

    real y_solve = integrate_ode_rk45(U_Pooled_ode, {0.0}, 0.0, rep_array(time, 1), parms, {0.0}, x_i)[1, 1];
    return y_solve;
  }

  // integral across age values -- 0.0 to time
  real U_host_t(real[] parms){
    int x_i[0];
    real time = parms[7];

    real y_solve = integrate_ode_rk45(U_host_ode, {0.0}, 0.0, rep_array(time, 1), parms, {0.0}, x_i)[1, 1];
    return y_solve;
  }

  // integral across age values -- 0.0 to time
  real U_donor_t(real[] parms){
    int x_i[0];
    real time = parms[7];

    real y_solve = integrate_ode_rk45(U_donor_ode, {0.0}, 0.0, rep_array(time, 1), parms, {0.0}, x_i)[1, 1];
    return y_solve;
  }

  // Vectorised function for theta div pool size
  real[] U_total_time(real[] time, real[] parms){
   int ndim = size(time);
   real y_solve[ndim];
   real params[6];
   params[1:5] = parms[1:5];

   for (i in 1:ndim){
     params[6] = time[i];
     y_solve[i] = U_total_t(params);
   }
   return y_solve;
  }

  // Vectorised function for theta div pool size
  real[] U_Pooled_time(real[] time, real[] tBMT, real[] parms){
   int ndim = size(time);
   real y_solve[ndim];
   real params[7];
   params[1:5] = parms[1:5];

   for (i in 1:ndim){
     params[6] = tBMT[i];
     params[7] = time[i];
     y_solve[i] = U_Pooled_t(params);
   }
   return y_solve;
  }

  // Vectorised function for theta div pool size
  real[] U_host_time(real[] time, real[] tBMT, real[] parms){
   int ndim = size(time);
   real y_solve[ndim];
   real params[7];
   params[1:5] = parms[1:5];

   for (i in 1:ndim){
     params[6] = tBMT[i];
     params[7] = time[i];
     y_solve[i] = U_host_t(params);
   }
   return y_solve;
  }

  // Vectorised function for theta div pool size
  real[] U_donor_time(real[] time, real[] tBMT, real[] parms){
   int ndim = size(time);
   real y_solve[ndim];
   real params[7];
   params[1:5] = parms[1:5];

   for (i in 1:ndim){
     params[6] = tBMT[i];
     params[7] = time[i];
     y_solve[i] = U_donor_t(params);
   }
   return y_solve;
  }

  vector math_reduce(vector global_params, vector local_params, real[] x_r, int[] x_i){
    // data for each shard
    int n = size(x_i); // n = 1
    int dat_t0 = x_i[1];                          // time zero -- for chimeras age at BMT

    //PDE solution
    real chi_counts_mean[n];
    real host_counts_mean[n];
    real donor_counts_mean[n];
    real host_ki_counts[n];
    real host_ki_mean;
    real donor_ki_counts[n];
    real donor_ki_mean;

    vector[4*n] y_mean_stacked;
    // each shard has a single datpoint so its unique ****
    // PDE solution for chimera dataset -- x_r = data time and x_i = time at BMT
    chi_counts_mean = N_pooled_time(x_r,  to_array_1d(to_vector(x_i)/1.0), to_array_1d(global_params));
    host_counts_mean = N_host_time(x_r, to_array_1d(to_vector(x_i)/1.0), to_array_1d(global_params));
    donor_counts_mean = N_donor_time(x_r, to_array_1d(to_vector(x_i)/1.0), to_array_1d(global_params));
    host_ki_counts = U_host_time(x_r,  to_array_1d(to_vector(x_i)/1.0), to_array_1d(global_params));
    host_ki_mean = host_ki_counts[1]/host_counts_mean[1];
    donor_ki_counts = U_donor_time(x_r,  to_array_1d(to_vector(x_i)/1.0), to_array_1d(global_params));
    donor_ki_mean = donor_ki_counts[1]/donor_counts_mean[1];
    y_mean_stacked[1] = chi_counts_mean[1];
    y_mean_stacked[2] = donor_counts_mean[1]/(chi_counts_mean[1] * Chi_spline(x_r[1] - (x_i[1]/1.0)));
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

parameters{
  real<lower=1E5, upper=3E6> N0;                  // total cells counts at t0
  real p_age;
  real<lower=0.0, upper=1.0> delta;
  real<lower=0.0, upper=1.0> rho;
  real r_del;

  real<lower=0> sigma_counts;
  real<lower=0> sigma_Nfd;
  real<lower=0> sigma_donor_ki;
  real<lower=0> sigma_host_ki;
}

transformed parameters{
  vector[5] global_params;
  vector[n_solve] y3_solve;               // PDE prediction for counts from chimera data
  vector[n_solve] y4_solve;               // PDE prediction for Nfd from chimera data
  vector[n_solve] y5_solve;               // PDE prediction for ki proportions in donor compartment from chimera data
  vector[n_solve] y6_solve;               // PDE prediction for ki proportions in host compartment from chimera data
  vector[(4*n_solve)] y_mean_stacked;        // compliled output across all nodes

  vector[numObs1] counts_naive_mean;               // ODE predictions for naive Treg counts in thymus
  vector[numObs2] Nfd_naive_mean;                  // ODE predictions for naive Treg Nfd in thymus
  vector[numObs3] ki_donor_naive_mean;             // ODE predictions for naive Treg donor ki67 proportions in thymus
  vector[numObs4] ki_host_naive_mean;              // ODE predictions for naive Treg host ki67 proportions in thymus

  global_params[1] = N0;
  global_params[2] = p_age;
  global_params[3] = delta;
  global_params[4] = rho;
  global_params[5] = r_del;

  // combining the output from all the shards
  y_mean_stacked = map_rect(math_reduce, global_params, local_params, x_r, x_i);

  for (i in 1:n_solve){
    y3_solve[i] = y_mean_stacked[4*i - 3];
    y4_solve[i] = y_mean_stacked[4*i - 2];
    y5_solve[i] = y_mean_stacked[4*i - 1];
    y6_solve[i] = y_mean_stacked[4*i];
  }

  for (i in 1:numObs1){
    counts_naive_mean[i] = y3_solve[time_index_counts[i]];
  }
  for (i in 1:numObs2){
    Nfd_naive_mean[i] = y4_solve[time_index_chi[i]];
  }
  for (i in 1:numObs3){
    ki_donor_naive_mean[i] = y5_solve[time_index_donorki[i]];
  }
  for (i in 1:numObs4){
    ki_host_naive_mean[i] = y6_solve[time_index_hostki[i]];
  }
}

model{
  N0 ~ normal(1e6, 3E5);
  p_age ~ normal(0.0, 0.3);
  delta ~ normal(0.05, 0.3);
  rho ~ normal(0.005, 0.3);
  r_del ~ normal(0.0, 0.3);

  sigma_counts ~ normal(0.4, 0.1);
  sigma_Nfd ~ normal(0.2, 0.05);
  sigma_donor_ki ~ normal(0.03, 0.05);
  sigma_host_ki ~ normal(0.03, 0.05);

  log(counts_naive) ~ normal(log(counts_naive_mean), sigma_counts);
  asinsqrt_array(Nfd_naive) ~ normal(asinsqrt_array(to_array_1d(Nfd_naive_mean)), sigma_Nfd);
  asinsqrt_array(ki_donor_naive) ~ normal(asinsqrt_array(to_array_1d(ki_donor_naive_mean)), sigma_donor_ki);
  asinsqrt_array(ki_host_naive) ~ normal(asinsqrt_array(to_array_1d(ki_host_naive_mean)), sigma_host_ki);
}


generated quantities{
  real y_chi_pred1[numPred, 4];
  real y_chi_pred2[numPred, 4];
  real y_chi_pred3[numPred, 4];
  real y_chi_pred4[numPred, 4];

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

  real host_counts_pred1[numPred]; real host_counts_pred2[numPred];  real host_counts_pred3[numPred]; real host_counts_pred4[numPred];
  real donor_counts_pred1[numPred]; real donor_counts_pred2[numPred];  real donor_counts_pred3[numPred]; real donor_counts_pred4[numPred];

  real host_ki_pred1[numPred]; real host_ki_pred2[numPred];  real host_ki_pred3[numPred];  real host_ki_pred4[numPred];
  real donor_ki_pred1[numPred]; real donor_ki_pred2[numPred];  real donor_ki_pred3[numPred]; real donor_ki_pred4[numPred];

  vector[numObs1] log_lik_chi_counts;
  vector[numObs2] log_lik_Nfd;
  vector[numObs3] log_lik_donor_ki;
  vector[numObs4] log_lik_host_ki;

  // PDE solution -- predictions for total counts, Nfd, donor_ki, host_ki
  counts_naive_mean_pred1 = N_pooled_time(ts_pred1,  tb_pred1, to_array_1d(global_params));
  counts_naive_mean_pred2 = N_pooled_time(ts_pred2,  tb_pred2, to_array_1d(global_params));
  counts_naive_mean_pred3 = N_pooled_time(ts_pred3,  tb_pred3, to_array_1d(global_params));
  counts_naive_mean_pred4 = N_pooled_time(ts_pred4,  tb_pred4, to_array_1d(global_params));

  host_counts_pred1 = N_host_time(ts_pred1, tb_pred1, to_array_1d(global_params));
  host_counts_pred2 = N_host_time(ts_pred2, tb_pred2, to_array_1d(global_params));
  host_counts_pred3 = N_host_time(ts_pred3, tb_pred3, to_array_1d(global_params));
  host_counts_pred4 = N_host_time(ts_pred4, tb_pred4, to_array_1d(global_params));

  donor_counts_pred1 = N_donor_time(ts_pred1, tb_pred1, to_array_1d(global_params));
  donor_counts_pred2 = N_donor_time(ts_pred2, tb_pred2, to_array_1d(global_params));
  donor_counts_pred3 = N_donor_time(ts_pred3, tb_pred3, to_array_1d(global_params));
  donor_counts_pred4 = N_donor_time(ts_pred4, tb_pred4, to_array_1d(global_params));

  host_ki_pred1 = U_host_time(ts_pred1, tb_pred1, to_array_1d(global_params));
  host_ki_pred2 = U_host_time(ts_pred2, tb_pred2, to_array_1d(global_params));
  host_ki_pred3 = U_host_time(ts_pred3, tb_pred3, to_array_1d(global_params));
  host_ki_pred4 = U_host_time(ts_pred4, tb_pred4, to_array_1d(global_params));

  donor_ki_pred1 = U_donor_time(ts_pred1, tb_pred1, to_array_1d(global_params));
  donor_ki_pred2 = U_donor_time(ts_pred2, tb_pred2, to_array_1d(global_params));
  donor_ki_pred3 = U_donor_time(ts_pred3, tb_pred3, to_array_1d(global_params));
  donor_ki_pred4 = U_donor_time(ts_pred4, tb_pred4, to_array_1d(global_params));


  for (i in 1:numPred){
    Nfd_naive_mean_pred1[i] = donor_counts_pred1[i]/(counts_naive_mean_pred1[i] * Chi_spline(ts_pred1[i] - tb_pred1[i]));
    Nfd_naive_mean_pred2[i] = donor_counts_pred2[i]/(counts_naive_mean_pred2[i] * Chi_spline(ts_pred2[i] - tb_pred2[i]));
    Nfd_naive_mean_pred3[i] = donor_counts_pred3[i]/(counts_naive_mean_pred3[i] * Chi_spline(ts_pred3[i] - tb_pred3[i]));
    Nfd_naive_mean_pred4[i] = donor_counts_pred4[i]/(counts_naive_mean_pred4[i] * Chi_spline(ts_pred4[i] - tb_pred4[i]));

    ki_donor_naive_mean_pred1[i] = donor_ki_pred1[i]/donor_counts_pred1[i];
    ki_donor_naive_mean_pred2[i] = donor_ki_pred2[i]/donor_counts_pred2[i];
    ki_donor_naive_mean_pred3[i] = donor_ki_pred3[i]/donor_counts_pred3[i];
    ki_donor_naive_mean_pred4[i] = donor_ki_pred4[i]/donor_counts_pred4[i];

    ki_host_naive_mean_pred1[i] = host_ki_pred1[i]/host_counts_pred1[i];
    ki_host_naive_mean_pred2[i] = host_ki_pred2[i]/host_counts_pred2[i];
    ki_host_naive_mean_pred3[i] = host_ki_pred3[i]/host_counts_pred3[i];
    ki_host_naive_mean_pred4[i] = host_ki_pred4[i]/host_counts_pred4[i];

    counts_naive_pred1[i] = exp(normal_rng(log(counts_naive_mean_pred1[i]), sigma_counts));
    counts_naive_pred2[i] = exp(normal_rng(log(counts_naive_mean_pred2[i]), sigma_counts));
    counts_naive_pred3[i] = exp(normal_rng(log(counts_naive_mean_pred3[i]), sigma_counts));
    counts_naive_pred4[i] = exp(normal_rng(log(counts_naive_mean_pred4[i]), sigma_counts));

    Nfd_naive_pred1[i] = asinsqrt_inv(normal_rng(asinsqrt_real(Nfd_naive_mean_pred1[i]), sigma_Nfd));
    Nfd_naive_pred2[i] = asinsqrt_inv(normal_rng(asinsqrt_real(Nfd_naive_mean_pred2[i]), sigma_Nfd));
    Nfd_naive_pred3[i] = asinsqrt_inv(normal_rng(asinsqrt_real(Nfd_naive_mean_pred3[i]), sigma_Nfd));
    Nfd_naive_pred4[i] = asinsqrt_inv(normal_rng(asinsqrt_real(Nfd_naive_mean_pred4[i]), sigma_Nfd));

    ki_donor_naive_pred1[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_donor_naive_mean_pred1[i]), sigma_donor_ki));
    ki_donor_naive_pred1[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_donor_naive_mean_pred2[i]), sigma_donor_ki));
    ki_donor_naive_pred1[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_donor_naive_mean_pred3[i]), sigma_donor_ki));
    ki_donor_naive_pred1[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_donor_naive_mean_pred4[i]), sigma_donor_ki));

    ki_host_naive_pred1[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_host_naive_mean_pred1[i]), sigma_host_ki));
    ki_host_naive_pred1[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_host_naive_mean_pred2[i]), sigma_host_ki));
    ki_host_naive_pred1[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_host_naive_mean_pred3[i]), sigma_host_ki));
    ki_host_naive_pred1[i] = asinsqrt_inv(normal_rng(asinsqrt_real(ki_host_naive_mean_pred4[i]), sigma_host_ki));
  }

  // calculating log likelihoods
  for (i in 1:numObs1) {
    log_lik_chi_counts[i] = normal_lpdf(log(counts_naive[i]) | log(counts_naive_mean[i]), sigma_counts);
  }
  for (i in 1:numObs2) {
    log_lik_Nfd[i]        = normal_lpdf(asinsqrt_real(Nfd_naive[i]) | asinsqrt_real(Nfd_naive_mean[i]), sigma_Nfd);
  }
  for (i in 1:numObs4) {
    log_lik_host_ki[i]    = normal_lpdf(asinsqrt_real(ki_host_naive[i]) | asinsqrt_real(ki_host_naive_mean[i]), sigma_host_ki);
  }
  for (i in 1:numObs3) {
    log_lik_donor_ki[i]   = normal_lpdf(asinsqrt_real(ki_donor_naive[i]) | asinsqrt_real(ki_donor_naive_mean[i]), sigma_donor_ki);
  }
}
