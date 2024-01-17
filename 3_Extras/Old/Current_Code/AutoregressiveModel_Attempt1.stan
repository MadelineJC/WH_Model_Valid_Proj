// This script is syntactically correct, but re-estimates the parms at every time step
// y are data, Z are estimates
// This is syntactically, very similar to the global version of this script

functions{
  real[] dZ_dt(
    real t,
    real[] Z,
    real[] theta,
    real[] x_r,
    int[] x_i){
    real P = Z[1];
    real H = Z[2];

    real r = theta[1];
    real O = theta[2];
    real h = theta[3];
    real b = theta[4];
    real c = theta[5];
    real u = theta[6];

    real dP_dt = P*r - H*(O*P/(1 + O*h*P));
    real dH_dt = b + H*(c*(O*P/(1 + O*h*P))-u);
    return({dP_dt,dH_dt});
  }
}
data {
  int<lower=0>N; 
  real ts[N]; 
  real<lower=0>y[N,2]; 
}
parameters {
  real<lower=0>r;
  real<lower=0,upper=1>O;
  real<lower=0>h;
  real<lower=0>b;
  real<lower=0>c;
  real<lower=0,upper=1>u;
  real<lower=0>sigma[2];
}
transformed parameters{
real Z[N, 2];
for(i in 2:N){
  Z[i:i,:] = integrate_ode_rk45(dZ_dt, // Function
  y[i-1], // Initial value (empirical data point at previous time step)
  ts[i-1], // Initial time step
  ts[i:i], // Next time step (time step to be solved/estimated)
  {r, O, h, b, c, u},
  rep_array(0.0,2),rep_array(0,2),1e-10,1e-10,2e4);
  }
}
model {
  r~normal(2.5,1);
  O~normal(0.01,2);
  h~normal(0.07,2);
  b~normal(35,1);
  c~normal(0.3,1);
  u~normal(0.4,1);
  sigma~lognormal(-1, 1);
  for (k in 1:2) {
    y[ , k] ~ lognormal(log(Z[ , k]), sigma[k]);
  }
}
generated quantities {
  real y_rep[N, 2];
  for (k in 1:2) {
    for (n in 1:N)
      y_rep[n, k] = lognormal_rng(log(Z[n, k]), sigma[k]);
  }
}
