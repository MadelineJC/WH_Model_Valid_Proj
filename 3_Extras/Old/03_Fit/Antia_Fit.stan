functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real I = z[2];

    real r = theta[1];  
    real k = theta[2];
    real p = theta[3];
    real o = theta[4];

    real dP_dt = r*P - k*P*I;
    real dI_dt = p*I*(P/(P + o));
    return { dP_dt, dI_dt };
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
  real<lower = 0, upper = 1> k;
  real<lower = 0> p;
  real<lower = 0> o;
  real<lower = 0> z_init[2];  
  real<lower = 0> sigma[2];   
}
transformed parameters {
  real z[N, 2]
    = integrate_ode_bdf(dz_dt, z_init, 0, ts, {r, k, p, o},
                         rep_array(0.0, 0), rep_array(0, 0));
}
model {
  r ~ normal(1, 3); // r = 0.2
  k ~ lognormal(log(0.1), 1); // k = 0.01
  p ~ normal(1, 1); // p = 1
  o ~ normal(1000, 10); // o = 1000
  sigma ~ lognormal(-1, 1);
  z_init[1] ~ normal(1, 1); // Initial value of P 
  z_init[2] ~ normal(1, 1); // Initial value of I
  for (m in 1:2) {
    y_init[m] ~ lognormal(log(z_init[m]), sigma[m]);
    y[ , m] ~ lognormal(log(z[, m]), sigma[m]);
  }
}
generated quantities {
  real y_init_rep[2];
  real y_rep[N, 2];
  for (m in 1:2) {
    y_init_rep[m] = lognormal_rng(log(z_init[m]), sigma[m]);
    for (n in 1:N)
      y_rep[n, m] = lognormal_rng(log(z[n, m]), sigma[m]);
  }
}
