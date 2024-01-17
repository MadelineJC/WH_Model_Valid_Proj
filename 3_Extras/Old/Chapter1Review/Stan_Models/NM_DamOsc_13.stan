
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
  int y_init[2];            // initial measured populations
  int<lower = 0> y[N, 2];   // measured populations
}
parameters {
  real<lower = 0> r;
  real<lower = 0> O;
  real<lower = 0> h;
  real<lower = 0> b;
  real<lower = 0> c;
  real<lower = 0> u;
  real<lower = 0> z_init[2];  // initial population
}
transformed parameters {
  real z[N, 2]
    = integrate_ode_bdf(dz_dt, z_init, 0, ts, {r, O, h, b, c, u},
                         rep_array(0.0, 0), rep_array(0, 0),
                         1e-10, 1e-10, 1e8);
}
model {
  r ~ uniform(0, 10); // r = 2.5
  O ~ uniform(0, 1); // O = 0.008
  h ~ uniform(0, 1); // h = 0.06
  b ~ uniform(0, 1000); // b = 35
  c ~ uniform(0, 1); // c = 0.2, u = 0.2
  u ~ uniform(0, 1); // c = 0.2, u = 0.2
  z_init ~ lognormal(log(10), 1);
  for (k in 1:2) {
    y_init[k] ~ poisson(z_init[k]);
    y[ ,k] ~ poisson(z[ ,k]);
  }
}
generated quantities {
  real y_init_rep[2];
  real y_rep[N, 2];
  for (k in 1:2) {
    y_init_rep[k] = poisson_rng(z_init[k]);
    for (n in 1:N)
      y_rep[n, k] = poisson_rng(z[n, k]);
  }
}      

