
data {
  int<lower = 0> N;
  real P[N];
  real I[N];
}
parameters {
  real<lower = 0> r; // Replication rate of parasite
  real<lower = 0> B; // Recognition rate of parasite by host immune system
  real<lower = 0> h; // Handling time of parasite by immune system
  real<lower = 0> b; // Immigration rate of immune cells in absence of infection
  real<lower = 0> e; // Activation/proliferation rate of host immune system
  real<lower = 0> delta; // Natural mortality rate of host immune cells
  real<lower=0> sigma; // Error
}
model {
  r ~ normal(3, 1); // 3.0
  B ~ normal(0, 1); // 0.0015
  h ~ normal(0, 1); // 0.2
  b ~ normal(35, 1); // 35
  e ~ normal(0, 1); // 0.2
  delta ~ normal(0, 1); // 0.5
  sigma ~ lognormal(-1, 1);
  for (t in 2:N) {
    P[t] ~ normal((1 + r*0.1 - (1-exp((-B*P[t - 1]/(1 + B*h*P[t - 1]^2))*I[t - 1]*0.1)))*P[t - 1], sigma);
    I[t] ~ normal((1 + b*0.1*(1/I[t - 1]) + e*(B*P[t - 1]^2/(1 + B*P[t - 1]^2*h))*0.1 - (1-exp(-delta*0.1)))*I[t - 1], sigma);
  }
}
