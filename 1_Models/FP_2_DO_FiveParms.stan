
functions {
  real[] dz_dt(real t,       // time
               real[] z,     // system state {prey, predator}
               real[] theta, // parameters
               real[] x_r,   // unused data
               int[] x_i) {
    real P = z[1];
    real I = z[2];

    real r = theta[1];  
    real h = theta[2];
    real b = theta[3];
    real e1 = theta[4];
    real delta = theta[5];

    real dP_dt = P*r - I*(0.008*P/(1 + 0.008*h*P));
    real dI_dt = b + I*(e1*(0.008*P/(1 + 0.008*h*P)) - delta);

    return { dP_dt, dI_dt };
  }
}
data {
  int<lower = 0> N;          // number of measurement times
  real ts[N];                // measurement times > 0
  real y_init[2];            // initial measured populations
  real<lower = 0> y[N, 2];   // measured populations
}
parameters {
  real<lower = 0> theta[5];   // {r, h, b, e1, delta}
  real<lower = 0> z_init[2];  // initial population
  real<lower = 0> sigma[2];   // measurement errors
}
transformed parameters {
  real z[N, 2]
    = integrate_ode_bdf(dz_dt, z_init, 0, ts, theta,
                         rep_array(0.0, 0), rep_array(0, 0));
}
model {
  theta[{1}] ~ uniform(0, 10); // r
  theta[{2}] ~ uniform(0, 1); // h
  theta[{3}] ~ uniform(0, 1000); // b
  theta[{4, 5}] ~ uniform(0, 1); // e1, delta
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

