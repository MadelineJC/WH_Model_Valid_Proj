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
    real dP_dt = P*r - H*(O*P/(1 + O*P*h));
    real dH_dt = b + H*(c*(O*P/(1 + O*P*h))-u);
    return({dP_dt,dH_dt});
  }
}
data {
  int<lower=0> N;           
  real ts[N];                 
  real y_init[2];             
  real<lower=0> y[N, 2];    
}
parameters {
  real<lower=0>r;
  real<lower=0,upper=1>O;
  real<lower=0>h;
  real<lower=0>b;
  real<lower=0>c;
  real<lower=0,upper=1>u;
  real<lower=0>Z_init[2];  
  real<lower=0>sigma[2];   
}
transformed parameters {
  real Z[N, 2]
    = integrate_ode_rk45(dZ_dt, Z_init, 0, ts, {r, O, h, b, c, u},
                         rep_array(0.0, 2), rep_array(0, 2),
                         1e-10, 1e-10, 5e2);
}
model {
  r~normal(2.5,1);
  O~normal(0.015,0.5); 
  h~normal(0.075,0.5); 
  b~normal(35,1); 
  c~normal(0.3,0.5); 
  u~normal(0.4,0.5);
  sigma~lognormal(-1, 1);
  Z_init~lognormal(log(3), 1);
  for (k in 1:2) {
    y_init[k] ~ lognormal(log(Z_init[k]), sigma[k]);
    y[ , k] ~ lognormal(log(Z[, k]), sigma[k]);
  }
}
generated quantities {
  real y_init_rep[2];
  real y_rep[N, 2];
  for (k in 1:2) {
    y_init_rep[k] = lognormal_rng(log(Z_init[k]), sigma[k]);
    for (n in 1:N)
      y_rep[n, k] = lognormal_rng(log(Z[n, k]), sigma[k]);
  }
}
