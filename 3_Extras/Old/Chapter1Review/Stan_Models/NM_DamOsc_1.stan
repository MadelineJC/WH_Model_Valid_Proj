
functions {
  real[] dz_dt(real t,       // time
               real[] z,     // system state {prey, predator}
               real[] theta, // parameters
               real[] x_r,   // unused data
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real O = theta[2];
    real h = theta[3];
    real b = theta[4];
    real c = theta[5];
    real u = theta[6];

    real dP_dt = P*r - H*(O*P/(1 + O*h*P));
    real dH_dt = b + H*(c*(O*P/(1 + O*h*P))-u);

    return { dP_dt, dH_dt };
  }
}
data {
  int<lower = 0> N;          // number of measurement times
  real ts[N];                // measurement times > 0
  real y_init[2];            // initial measured populations
  real<lower = 0> y[N, 2];   // measured populations
}
parameters {
  real<lower = 0> theta[6];   // {r, O, h, b, c, u}
  real<lower = 0> z_init[2];  // initial population
  real<lower = 0> sigma[2];   // measurement errors
}
transformed parameters {
  real z[N, 2]
    = integrate_ode_rk45(dz_dt, z_init, 0, ts, theta,
                         rep_array(0.0, 0), rep_array(0, 0),
                         1e-6, 1e-6, 1e6);
}
model {
  theta[{1}] ~ uniform(0, 10); // r = 2.5
  theta[{2, 3}] ~ uniform(0, 1); // O = 0.008, h = 0.06
  theta[{4}] ~ uniform(0, 1000); // b = 35
  theta[{5, 6}] ~ uniform(0, 1); // c = 0.2, u = 0.2
  sigma ~ lognormal(-1, 1);
  z_init ~ lognormal(log(10), 1);
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
