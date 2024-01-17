
data {
  int<lower=0> N;
  real P[N];
  real H[N];
}
parameters {
  real<lower = 0> r; // Replication rate of parasite
  // real<lower = 0> O; // Recognition rate of parasite by host immune system
  // real<lower = 0> h; // Handling time of parasite by immune system
  real<lower = 0> b; // Immigration rate of immune cells in absence of infection
  real<lower = 0> c; // Activation/proliferation rate of host immune system
  real<lower = 0> u; // Natural mortality rate of host immune cells
  real<lower=0> sigma; // Error
}
model {
  r ~ uniform(0, 10); // 2.5
  // O ~ uniform(0, 1); // 0.008
  // h ~ uniform(0, 1); // 0.06
  b ~ uniform(0, 1000); // 35
  c ~ uniform(0, 1); // 0.2
  u ~ uniform(0, 1); // 0.2
  sigma ~ lognormal(-1, 1);
  for (t in 2:N) {
    P[t] ~ normal(P[t-1] + r*P[t-1] - (1-exp((-0.008/(1 + 0.008*0.06*P[t-1]))*H[t-1]))*P[t-1], sigma);
    H[t] ~ normal(H[t-1] + b + c*(0.008*P[t-1]/(1 + 0.008*P[t-1]*0.06))*H[t-1] - (1-exp(-u))*H[t-1], sigma);
  }
}

