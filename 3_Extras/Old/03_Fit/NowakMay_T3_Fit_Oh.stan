functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];
    
    real r = theta[1];  
    //real O = theta[2];
    //real h = theta[3];
    real b = theta[2];
    real c = theta[3];
    real u = theta[4];
    
    real dP_dt = P*r - H*(0.0015*(P^2)/(1 + 0.0015*0.2*(P^2)));
    real dH_dt = b + H*(c*(0.0015*(P^2)/(1 + 0.0015*0.2*(P^2)))-u);
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
  real<lower = 0> theta[4];   
  real<lower = 0> z_init[2];  
  real<lower = 0> sigma[2];   
}
transformed parameters {
  real z[N, 2]
  = integrate_ode_rk45(dz_dt, z_init, 0, ts, theta,
                       rep_array(0.0, 0), rep_array(0, 0),
                       1e-5, 1e-3, 5e2);
}
model {
  theta[{1}] ~ normal(3, 1); // r = 3
  theta[{2}] ~ uniform(0,100); // b = 35
  theta[{3}] ~ normal(0.5,0.5); // c = 0.2
  theta[{4}] ~ normal(0.5,0.5); // u = 0.5
  sigma ~ lognormal(-1, 1);
  z_init ~ lognormal(log(140), 1);
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
