
data {
  int<lower=0> N; // Number of observations
  //real ts[N];
  //real y_init[1];             
  //real<lower = 0> y[N, 1]; // Outcome
  vector[N] y;
}
parameters {
  real<lower = 0> r; // Growth rate
  real<lower = 0> K; // Carrying capacity
  real<lower=0> sigma; // Error
}
model {
  r ~ normal(0, 10);
  K ~ normal(1000,10);
  //sigma ~ lognormal(-1, 1);
  for (t in 2:N) {
    y[t] ~ normal(y[t-1] + r*y[t-1]*(1-y[t-1]/K), sigma);
    //real mu = y[t-1] + r*y[t-1]*(1-y[t-1]/K);
      //y[t] ~ normal(mu, sigma);
  }
}
