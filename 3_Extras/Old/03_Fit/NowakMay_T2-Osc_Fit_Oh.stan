
      functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];
    
    real r = theta[1];  
    // real O = theta[2];
    // real h = theta[3];
    real b = theta[2];
    real c = theta[3];
    real u = theta[4];
    
    real dP_dt = P*r - H*(0.008*P/(1 + 0.008*0.06*P));
    real dH_dt = b + H*(c*(0.008*P/(1 + 0.008*0.06*P))-u);
    return { dP_dt, dH_dt };
  }
}
data {
  int<lower = 0> N;           
  real ts[N];                 
  real y_init[2];             
  real<lower = 0> y[N, 2];    
}
parameters {
  real<lower = 0> r;
  // real<lower = 0> O;
  // real<lower = 0> h;
  real<lower = 0> b;
  real<lower = 0> c;
  real<lower = 0> u;
  real<lower = 0> z_init[2];  
  real<lower = 0> sigma[2];   
}
transformed parameters {
  real z[N, 2]
  = integrate_ode_bdf(dz_dt, z_init, 0, ts, {r, b, c, u},
                       rep_array(0.0, 0), rep_array(0, 0));
}
model {
  r ~ normal(1, 3); // r = 2.5
  // O ~ normal(0, 1); // O = 0.012
  // h ~ normal(0, 1); // h = 0.075
  b ~ uniform(0, 100); // b = 35
  c ~ normal(0, 1); // c = 0.3
  u ~ normal(0, 1); // u = 0.41
  sigma ~ lognormal(-1, 1);
  z_init[1] ~ normal(80, 1); // 80
  z_init[2] ~ normal(200, 1); // 200
  for (k in 1:2) {
    y_init[k] ~ lognormal(log(z_init[k]), sigma[k]);
    y[ , k] ~ lognormal(log(z[, k]), sigma[k]);
  }
}
generated quantities {
  real y_init_rep[2];
  real y_rep[N, 2];
  for (k in 1:2) {
    y_init_rep[k] = lognormal_rng(log(z_init[k]), sigma[k]);
    for (n in 1:N)
      y_rep[n, k] = lognormal_rng(log(z[n, k]), sigma[k]);
  }
}
