// THIS IS WRITTEN IN C++ AND IS SYNTACTICALLY CORRECT

// Need to do model fitting, model comparison, and parameter estimation, in the order.

// I have no idea how to do any of this, so here we go lmao

functions{
  real[] dZ_dt(real t, // Time
               real[] Z, // System state {Parasite, Host}
               real[] alpha, // Parameters
               real[] x_r, // Unused data (?)
               int[] x_i){
    real P = Z[1]; // System state coded as an array, such that Z = (P,H)
    real H = Z[2];
  
    real r = alpha[1]; // Parameters of the system, in order they appear
    real O = alpha[2];
    real h = alpha[3];
    real b = alpha[4];
    real c = alpha[5];
    real u = alpha[6];
  
    real dP_dt = P*r - H*(O*P/1 + O*P*h); // Deterministic mechanistic model
    real dH_dt = b + H*(c*(O*P/1 + O*P*h)-u);
    return{dP_dt,dH_dt}; // Return the system state
  }
}

data{
  int<lower=0>N; // Define N as non-negative integer
  real ts[N]; // Assigns time points to N (I think?)
  real y_init[2];
  real<lower=0>y[N,2]; // Define y as real and non-negative
}

parameters{
  real<lower=0>alpha[6]; // Make all items in alpha non-neg
  real<lower=0>Z_init[2]; // Initial population size non-neg
  real<lower=0>sigma[2]; // Error term non-negative
}

transformed parameters{
  real Z[N,2]
  = integrate_ode_rk45(dZ_dt,Z_init,0,ts,alpha,rep_array(0.0,0),rep_array(0,0),1e-6,1e-5,2000);
}

model{
  alpha[{1}]~uniform(0,10);
  alpha[{2}]~uniform(0,1);
  alpha[{3}]~uniform(0,60);
  alpha[{4}]~uniform(0,100);
  alpha[{5}]~uniform(0,1);
  alpha[{6}]~uniform(0,1);
  sigma~lognormal(-1,1);
  Z_init~lognormal(log(10),1);
  for (k in 1:2){
    y_init[k]~lognormal(log(Z_init[k]),sigma[k]);
    y[ ,k]~lognormal(log(Z[ ,k]),sigma[k]);
  }
}
