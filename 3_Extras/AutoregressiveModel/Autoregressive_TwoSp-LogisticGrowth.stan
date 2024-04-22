
data {
  int<lower=0> N; // Number of observations
  vector [N] Y;
  vector [N] y;
}
parameters {
  real<lower = 0> R; // Growth rate of pop'n 1 (Y)
  real<lower = 0> K; // Carrying capacity of pop'n 1 (Y)
  real<lower = 0> r; // Growth rate of pop'n 2 (y)
  real<lower = 0> k; // Carrying capacity of pop'n 2 (y)
  real<lower=0> sigma; // Error
}
model {
  R ~ normal(0.1, 1);
  K ~ normal(1000, 1);
  r ~ normal(0.2, 1);
  k ~ normal(1200, 1);
  //sigma ~ lognormal(-1, 1);
  for (t in 2:N) {
    Y[t] ~ normal(Y[t-1] + R*Y[t-1]*(1-Y[t-1]/K), sigma); // Pop'n 1
    y[t] ~ normal(y[t-1] + r*y[t-1]*(1-y[t-1]/k), sigma); // Pop'n 2
  }
}
