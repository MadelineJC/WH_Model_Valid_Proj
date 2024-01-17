// I have no idea what I'm doing but I'm TRYING, OKAY

functions{
  real[] dZ_dt(real t, // Time
               real[] Z, // System state {Parasite, Host}
               real[] alpha, // Parameters
               real[] x_r, // Real data; below is integer data 
               int[] x_i){
    real P = Z[1]; // System state coded as an array, such that Z = (P,H)
    real H = Z[2];
  
    real r = alpha[1]; // Parameters of the system, in order they appear
    real O = alpha[2];
    real h = alpha[3];
    real b = alpha[4];
    real c = alpha[5];
    real u = alpha[6];
  
    real dP_dt = P*r - H*(O*P/1 + O*P*h); // Mechanistic model
    real dH_dt = b + H*(c*(O*P/1 + O*P*h)-u);
    return{dP_dt,dH_dt}; // Return the system state
  }
}

data{
  int<lower=0>N; // Define N as non-negative integer
  real ts[N]; // Assigns time points to N 
  real y_init[2]; // Initial conditions for ODE
  real<lower=0>y[N,2]; // Define y as real and non-negative
}

parameters{
  real<lower=0>alpha[6]; // Make all items in alpha non-neg
  real<lower=0>Z_init[2]; // Initial population size non-neg
  real<lower=0>sigma[2]; // Error term non-negative
}

// ODE solver; uses the Runge Kutta Dopri algorithm
transformed parameters{
  real Z[N,2]
  = integrate_ode_rk45(dZ_dt,Z_init,0,ts,alpha,rep_array(0.0,0),rep_array(0,0),1e-6,1e-5,2000);
}

model{
  alpha[{1}]~normal(0,1); // r
  alpha[{2}]~normal(0,1); // O
  alpha[{3}]~normal(0,1); // h
  alpha[{4}]~normal(0,1); // b
  alpha[{5}]~normal(0,1); // c
  alpha[{6}]~normal(0,1); // u
  sigma~lognormal(-1,1);
  Z_init~lognormal(log(10),1);
  for (k in 1:2){
    y_init[k]~normal(log(Z_init[k]),sigma[k]);
    y[ ,k]~normal(log(Z[ ,k]),sigma[k]);
  }
}
