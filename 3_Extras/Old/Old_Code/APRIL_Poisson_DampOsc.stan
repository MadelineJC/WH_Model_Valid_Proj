functions{
  real[] dZ_dt(
    real t, // Time
    real[] Z, // System state {Parasite, Host}
    real[] theta, // Parms
    real[] x_r, // Real data; below is integer data
    int[] x_i){
    real P = Z[1]; // System state coded as an array, such that Z = (P,H)
    real H = Z[2];

    real r = theta[1];  // Parameters of the system, in order they appear
    real O = theta[2];
    real h = theta[3];
    real b = theta[4];
    real c = theta[5];
    real u = theta[6];

    real dP_dt = P*r - H*(O*P/(1 + O*P*h)); // Mechanistic model
    real dH_dt = b + H*(c*(O*P/(1 + O*P*h))-u);
    return({dP_dt,dH_dt}); // Return the system state
  }
}

data{
  int<lower=0>N; // Define N as non-negative integer
  real ts[N]; // Assigns time points to N
  real y0[2]; // Initial conditions for ODE
  real t0; // Initial time point
  int<lower=0>y[N,2]; // Define y as real and non-negative
}

parameters{
  real<lower=0>r;
  real<lower=0,upper=1>O;
  real<lower=0>h;
  real<lower=0>b;
  real<lower=0>c;
  real<lower=0,upper=1>u;
}

transformed parameters{
  // Stiff solver for ODE
  real Z[N,2]
  = integrate_ode_bdf(dZ_dt,y0,t0,ts,{r, O, h, b, c, u},rep_array(0.0,2),rep_array(0,2),1e-10,1e-10,2e4);
  for (i in 1:N) {
    Z[i,1] = Z[i,1] + 1e-6;
    Z[i,2] = Z[i,2] + 1e-6;
  }
}

model{ // 
  r~normal(2.5,1); // r
  O~beta(2,2); // O; bounded between 0 and 1
  h~normal(0.06,1); // h
  b~normal(35,1); // b
  c~normal(0.2,1); // c
  u~beta(2,2); // u; bounded between 0 and 1
  for (k in 1:2){
    y[ ,k]~poisson(Z[ ,k]);
  }
}

generated quantities {
  int y_rep[N, 2];
  for (k in 1:2) {
    for (n in 1:N)
      y_rep[n,k] = poisson_rng(Z[n,k]);
  }
}
