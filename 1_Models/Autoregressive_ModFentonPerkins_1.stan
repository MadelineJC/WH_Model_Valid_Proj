
data {
  int<lower = 0> N;
  real P[N];
  real I[N];
}
parameters {
  real<lower = 0> r; // Replication rate of parasite
  real<lower = 0> B; // Recognition rate of parasite by host immune system
  real<lower = 0> b; // Immigration rate of immune cells in absence of infection
  real<lower = 0> e; // Activation/proliferation rate of host immune system
  real<lower = 0> delta; // Natural mortality rate of host immune cells
  real<lower = 0> sigma; // Error
}
model {
  r ~ uniform(0, 2); // 1.5
  B ~ uniform(0, 1); // 0.001
  b ~ uniform(0, 300); // 200
  e ~ uniform(0, 1); // 0.9
  delta ~ uniform(0, 1); // 0.41
  sigma ~ normal(0, 1);
  for (t in 2:N) {
    P[t] ~ normal((1 + r*0.1 - (1-exp(-B*I[t - 1]*0.1)))*P[t - 1], sigma);
    I[t] ~ normal((1 + b*0.1*(1/I[t - 1]) + e*(B*P[t - 1])*0.1 - (1-exp(-delta*0.1)))*I[t - 1], sigma);
  }
}
