
data {
  int<lower = 0> N; // Number of observations
  real<lower = 0> P[N];
  real<lower = 0> I[N];
}
parameters {
  real<lower = 0> r; // Replication rate
  real<lower = 0> k; // Rate of destruction of parasite
  real<lower = 0> p; // Max. growth rate of immunity
  real<lower = 0> o; // Parasite density at which immunity slows
  real<lower = 0> sigma; // Error
}
model {
  r ~ normal(0, 1); // 0.2
  k ~ normal(0, 1); // 0.01
  p ~ normal(1, 1); // 1
  o ~ normal(1000, 1); // 1000
  // sigma ~ lognormal(-1, 1);
  for (t in 2:N) {
    P[t] ~ normal(P[t-1] + r*P[t-1] - (1 - exp(-k*I[t-1]))*P[t - 1], sigma);
    I[t] ~ normal(I[t-1] + p*I[t-1]*(P[t-1]/(P[t-1] + o)), sigma);
  }
}
