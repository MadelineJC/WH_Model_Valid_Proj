library(rstan)
library(gdata)
library(bayesplot)
library(tidyverse)
library(brms)
library(beepr)
library(MCMCvis)

#### Data org. ####
NewData = read_csv('/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Data/NewData.csv')
N <- length(NewData$t) - 1
ts <- 1:N
y_init <- c(NewData$P[1], NewData$H[1])
y <- as.matrix(NewData[2:(N + 1), 2:3])
y <- cbind(y[ , 1], y[ , 2]); 
Data <- list(N = N, ts = ts, y_init = y_init, y = y)

#### Est. all ####

write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
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
  int<lower = 0> N;           
  real ts[N];                 
  real y_init[2];             
  real<lower = 0> y[N, 2];    
}
parameters {
  real<lower = 0> theta[6];   
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
  theta[{1}] ~ normal(2.5, 1);
  log(theta[{2}]) ~ normal(0.5,0.5);
  log(theta[{3}]) ~ normal(0.5,0.5);
  theta[{4}] ~ normal(35,1);
  theta[{5, 6}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model.stan")

Model <- stan_model("Global_Stan_Model.stan")
Fit_Stan <- sampling(Model, data = Data,
                     chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, 
                     seed = 123); beep(3)
Fit_Stan_Summ <- print(Fit_Stan, pars=c("theta", "sigma", "z_init"),
                         probs=c(0.1, 0.5, 0.9), digits = 3)

# chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4,

#### Give 1 parm ####

#### Give r ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real O = theta[1];  
    real h = theta[2];
    real b = theta[3];
    real c = theta[4];
    real u = theta[5];

    real dP_dt = P*2.5 - H*(O*P/(1 + O*h*P));
    real dH_dt = b + H*(c*(O*P/(1 + O*h*P))-u);
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
  real<lower = 0> theta[5];   
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
  //theta[{1}] ~ normal(2.5, 1);
  log(theta[{1}]) ~ normal(0.5,0.5);
  log(theta[{2}]) ~ normal(0.5,0.5);
  theta[{3}] ~ normal(35,1);
  theta[{4, 5}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-r.stan")

Model <- stan_model("Global_Stan_Model_Give-r.stan")
Fit_r_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 123); beep(3)
Fit_r_Stan_Summ <- print(Fit_r_Stan, pars=c("theta", "sigma", "z_init"),
                  probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give O ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real h = theta[2];
    real b = theta[3];
    real c = theta[4];
    real u = theta[5];

    real dP_dt = P*r - H*(0.012*P/(1 + 0.012*h*P));
    real dH_dt = b + H*(c*(0.012*P/(1 + 0.012*h*P))-u);
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
  real<lower = 0> theta[5];   
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
  theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1}]) ~ normal(0.5,0.5);
  log(theta[{2}]) ~ normal(0.5,0.5);
  theta[{3}] ~ normal(35,1);
  theta[{4, 5}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-O.stan")

Model <- stan_model("Global_Stan_Model_Give-O.stan")
Fit_O_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_O_Stan_Summ <- print(Fit_O_Stan, pars=c("theta", "sigma", "z_init"),
                         probs=c(0.1, 0.5, 0.9), digits = 3)


#### Give h ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real O = theta[2];
    real b = theta[3];
    real c = theta[4];
    real u = theta[5];

    real dP_dt = P*r - H*(O*P/(1 + O*0.075*P));
    real dH_dt = b + H*(c*(O*P/(1 + O*0.075*P))-u);
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
  real<lower = 0> theta[5];   
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
  theta[{1}] ~ normal(2.5, 1);
  log(theta[{2}]) ~ normal(0.5,0.5);
  theta[{3}] ~ normal(35,1);
  theta[{4, 5}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-h.stan")

Model <- stan_model("Global_Stan_Model_Give-h.stan")
Fit_h_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_h_Stan_Summ <- print(Fit_h_Stan, pars=c("theta", "sigma", "z_init"),
                         probs=c(0.1, 0.5, 0.9), digits = 3)

Y#### Give b ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real O = theta[2];
    real h = theta[3];
    real c = theta[4];
    real u = theta[5];

    real dP_dt = P*r - H*(O*P/(1 + O*h*P));
    real dH_dt = 35 + H*(c*(O*P/(1 + O*h*P))-u);
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
  real<lower = 0> theta[5];   
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
  theta[{1}] ~ normal(2.5, 1);
  log(theta[{2}]) ~ normal(0.5,0.5);
  log(theta[{3}]) ~ normal(0.5,0.5);
  //theta[{3}] ~ normal(35,1);
  theta[{4, 5}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-b.stan")

Model <- stan_model("Global_Stan_Model_Give-b.stan")
Fit_b_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 123); beep(3)
Fit_b_Stan_Summ <- print(Fit_b_Stan, pars=c("theta", "sigma", "z_init"),
                         probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give c ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real O = theta[2];
    real h = theta[3];
    real b = theta[4];
    real u = theta[5];

    real dP_dt = P*r - H*(O*P/(1 + O*h*P));
    real dH_dt = b + H*(0.3*(O*P/(1 + O*h*P))-u);
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
  real<lower = 0> theta[5];   
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
  theta[{1}] ~ normal(2.5, 1);
  log(theta[{2}]) ~ normal(0.5,0.5);
  log(theta[{3}]) ~ normal(0.5,0.5);
  theta[{4}] ~ normal(35,1);
  //theta[{4, 5}] ~ normal(0.5,0.5);
  theta[{5}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-c.stan")

Model <- stan_model("Global_Stan_Model_Give-c.stan")
Fit_c_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_c_Stan_Summ <- print(Fit_c_Stan, pars=c("theta", "sigma", "z_init"),
                         probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give u ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real O = theta[2];
    real h = theta[3];
    real b = theta[4];
    real c = theta[5];

    real dP_dt = P*r - H*(O*P/(1 + O*h*P));
    real dH_dt = b + H*(c*(O*P/(1 + O*h*P))-0.41);
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
  real<lower = 0> theta[5];   
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
  theta[{1}] ~ normal(2.5, 1);
  log(theta[{2}]) ~ normal(0.5,0.5);
  log(theta[{3}]) ~ normal(0.5,0.5);
  theta[{4}] ~ normal(35,1);
  //theta[{4, 5}] ~ normal(0.5,0.5);
  theta[{5}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-u.stan")

Model <- stan_model("Global_Stan_Model_Give-u.stan")
Fit_u_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 123); beep(3)
Fit_u_Stan_Summ <- print(Fit_u_Stan, pars=c("theta", "sigma", "z_init"),
                         probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give 2 parms ####

#### Give rO ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real h = theta[1];  
    real b = theta[2];
    real c = theta[3];
    real u = theta[4];

    real dP_dt = P*2.5 - H*(0.012*P/(1 + 0.012*h*P));
    real dH_dt = b + H*(c*(0.012*P/(1 + 0.012*h*P))-u);
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
  //theta[{1}] ~ normal(2.5, 1);
  log(theta[{1}]) ~ normal(0.5,0.5);
  //log(theta[{3}]) ~ normal(0.5,0.5);
  theta[{2}] ~ normal(35,1);
  theta[{3, 4}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rO.stan")

Model <- stan_model("Global_Stan_Model_Give-rO.stan")
Fit_rO_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rO_Stan_Summ <- print(Fit_rO_Stan, pars=c("theta", "sigma", "z_init"),
                         probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give rh ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real O = theta[1];  
    real b = theta[2];
    real c = theta[3];
    real u = theta[4];

    real dP_dt = P*2.5 - H*(O*P/(1 + O*0.075*P));
    real dH_dt = b + H*(c*(O*P/(1 + O*0.075*P))-u);
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
  //theta[{1}] ~ normal(2.5, 1);
  log(theta[{1}]) ~ normal(0.5,0.5);
  //log(theta[{3}]) ~ normal(0.5,0.5);
  theta[{2}] ~ normal(35,1);
  theta[{3, 4}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rh.stan")

Model <- stan_model("Global_Stan_Model_Give-rh.stan")
Fit_rh_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rh_Stan_Summ <- print(Fit_rh_Stan, pars=c("theta", "sigma", "z_init"),
                          probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give rb ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real O = theta[1];  
    real h = theta[2];
    real c = theta[3];
    real u = theta[4];

    real dP_dt = P*2.5 - H*(O*P/(1 + O*h*P));
    real dH_dt = 35 + H*(c*(O*P/(1 + O*h*P))-u);
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
  //theta[{1}] ~ normal(2.5, 1);
  log(theta[{1}]) ~ normal(0.5,0.5);
  log(theta[{2}]) ~ normal(0.5,0.5);
  //theta[{2}] ~ normal(35,1);
  theta[{3, 4}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rb.stan")

Model <- stan_model("Global_Stan_Model_Give-rb.stan")
Fit_rb_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rb_Stan_Summ <- print(Fit_rb_Stan, pars=c("theta", "sigma", "z_init"),
                          probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give rc ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real O = theta[1];  
    real h = theta[2];
    real b = theta[3];
    real u = theta[4];

    real dP_dt = P*2.5 - H*(O*P/(1 + O*h*P));
    real dH_dt = b + H*(0.3*(O*P/(1 + O*h*P))-u);
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
  //theta[{1}] ~ normal(2.5, 1);
  log(theta[{1}]) ~ normal(0.5,0.5);
  log(theta[{2}]) ~ normal(0.5,0.5);
  theta[{3}] ~ normal(35,1);
  theta[{4}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rc.stan")

Model <- stan_model("Global_Stan_Model_Give-rc.stan")
Fit_rc_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rc_Stan_Summ <- print(Fit_rc_Stan, pars=c("theta", "sigma", "z_init"),
                          probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give ru ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real O = theta[1];  
    real h = theta[2];
    real b = theta[3];
    real c = theta[4];

    real dP_dt = P*2.5 - H*(O*P/(1 + O*h*P));
    real dH_dt = b + H*(c*(O*P/(1 + O*h*P))-0.41);
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
  //theta[{1}] ~ normal(2.5, 1);
  log(theta[{1}]) ~ normal(0.5,0.5);
  log(theta[{2}]) ~ normal(0.5,0.5);
  theta[{3}] ~ normal(35,1);
  theta[{4}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-ru.stan")

Model <- stan_model("Global_Stan_Model_Give-ru.stan")
Fit_ru_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_ru_Stan_Summ <- print(Fit_ru_Stan, pars=c("theta", "sigma", "z_init"),
                          probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give Oh ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real b = theta[2];
    real c = theta[3];
    real u = theta[4];

    real dP_dt = P*r - H*(0.012*P/(1 + 0.012*0.075*P));
    real dH_dt = b + H*(c*(0.012*P/(1 + 0.012*0.075*P))-u);
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
  theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  theta[{2}] ~ normal(35,1);
  theta[{3,4}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-Oh.stan")

Model <- stan_model("Global_Stan_Model_Give-Oh.stan")
Fit_Oh_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_Oh_Stan_Summ <- print(Fit_Oh_Stan, pars=c("theta", "sigma", "z_init"),
                          probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give Ob ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real h = theta[2];
    real c = theta[3];
    real u = theta[4];

    real dP_dt = P*r - H*(0.012*P/(1 + 0.012*h*P));
    real dH_dt = 35 + H*(c*(0.012*P/(1 + 0.012*h*P))-u);
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
  theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1}]) ~ normal(0.5,0.5);
  log(theta[{2}]) ~ normal(0.5,0.5);
  //theta[{2}] ~ normal(35,1);
  theta[{3,4}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-Ob.stan")

Model <- stan_model("Global_Stan_Model_Give-Ob.stan")
Fit_Ob_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_Ob_Stan_Summ <- print(Fit_Ob_Stan, pars=c("theta", "sigma", "z_init"),
                          probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give Oc ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real h = theta[2];
    real b = theta[3];
    real u = theta[4];

    real dP_dt = P*r - H*(0.012*P/(1 + 0.012*h*P));
    real dH_dt = b + H*(0.3*(0.012*P/(1 + 0.012*h*P))-u);
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
  theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1}]) ~ normal(0.5,0.5);
  log(theta[{2}]) ~ normal(0.5,0.5);
  theta[{3}] ~ normal(35,1);
  theta[{4}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-Oc.stan")

Model <- stan_model("Global_Stan_Model_Give-Oc.stan")
Fit_Oc_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_Oc_Stan_Summ <- print(Fit_Oc_Stan, pars=c("theta", "sigma", "z_init"),
                          probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give Ou ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real h = theta[2];
    real b = theta[3];
    real c = theta[4];

    real dP_dt = P*r - H*(0.012*P/(1 + 0.012*h*P));
    real dH_dt = b + H*(c*(0.012*P/(1 + 0.012*h*P))-0.41);
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
  theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1}]) ~ normal(0.5,0.5);
  log(theta[{2}]) ~ normal(0.5,0.5);
  theta[{3}] ~ normal(35,1);
  theta[{4}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-Ou.stan")

Model <- stan_model("Global_Stan_Model_Give-Ou.stan")
Fit_Ou_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_Ou_Stan_Summ <- print(Fit_Ou_Stan, pars=c("theta", "sigma", "z_init"),
                          probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give hb ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real O = theta[2];
    real c = theta[3];
    real u = theta[4];

    real dP_dt = P*r - H*(O*P/(1 + O*0.075*P));
    real dH_dt = 35 + H*(c*(O*P/(1 + O*0.075*P))-u);
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
  theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1}]) ~ normal(0.5,0.5);
  log(theta[{2}]) ~ normal(0.5,0.5);
  //theta[{3}] ~ normal(35,1);
  theta[{3,4}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-hb.stan")

Model <- stan_model("Global_Stan_Model_Give-hb.stan")
Fit_hb_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_hb_Stan_Summ <- print(Fit_hb_Stan, pars=c("theta", "sigma", "z_init"),
                          probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give hc ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real O = theta[2];
    real b = theta[3];
    real u = theta[4];

    real dP_dt = P*r - H*(O*P/(1 + O*0.075*P));
    real dH_dt = b + H*(0.3*(O*P/(1 + O*0.075*P))-u);
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
  theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1}]) ~ normal(0.5,0.5);
  log(theta[{2}]) ~ normal(0.5,0.5);
  theta[{3}] ~ normal(35,1);
  theta[{4}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-hc.stan")

Model <- stan_model("Global_Stan_Model_Give-hc.stan")
Fit_hc_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_hc_Stan_Summ <- print(Fit_hc_Stan, pars=c("theta", "sigma", "z_init"),
                          probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give hu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real O = theta[2];
    real b = theta[3];
    real c = theta[4];

    real dP_dt = P*r - H*(O*P/(1 + O*0.075*P));
    real dH_dt = b + H*(c*(O*P/(1 + O*0.075*P))-0.41);
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
  theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1}]) ~ normal(0.5,0.5);
  log(theta[{2}]) ~ normal(0.5,0.5);
  theta[{3}] ~ normal(35,1);
  theta[{4}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-hu.stan")

Model <- stan_model("Global_Stan_Model_Give-hu.stan")
Fit_hu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_hu_Stan_Summ <- print(Fit_hu_Stan, pars=c("theta", "sigma", "z_init"),
                          probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give bc ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real O = theta[2];
    real h = theta[3];
    real u = theta[4];

    real dP_dt = P*r - H*(O*P/(1 + O*h*P));
    real dH_dt = 35 + H*(0.3*(O*P/(1 + O*h*P))-u);
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
  theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1}]) ~ normal(0.5,0.5);
  log(theta[{2,3}]) ~ normal(0.5,0.5);
  //theta[{3}] ~ normal(35,1);
  theta[{4}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-bc.stan")

Model <- stan_model("Global_Stan_Model_Give-bc.stan")
Fit_bc_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_bc_Stan_Summ <- print(Fit_bc_Stan, pars=c("theta", "sigma", "z_init"),
                          probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give bu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real O = theta[2];
    real h = theta[3];
    real c = theta[4];

    real dP_dt = P*r - H*(O*P/(1 + O*h*P));
    real dH_dt = 35 + H*(c*(O*P/(1 + O*h*P))-0.41);
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
  theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1}]) ~ normal(0.5,0.5);
  log(theta[{2,3}]) ~ normal(0.5,0.5);
  //theta[{3}] ~ normal(35,1);
  theta[{4}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-bu.stan")

Model <- stan_model("Global_Stan_Model_Give-bu.stan")
Fit_bu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_bu_Stan_Summ <- print(Fit_bu_Stan, pars=c("theta", "sigma", "z_init"),
                          probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give cu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real O = theta[2];
    real h = theta[3];
    real b = theta[4];

    real dP_dt = P*r - H*(O*P/(1 + O*h*P));
    real dH_dt = b + H*(0.3*(O*P/(1 + O*h*P))-0.41);
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
  theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1}]) ~ normal(0.5,0.5);
  log(theta[{2,3}]) ~ normal(0.5,0.5);
  theta[{4}] ~ normal(35,1);
  //theta[{4}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-cu.stan")

Model <- stan_model("Global_Stan_Model_Give-cu.stan")
Fit_cu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_cu_Stan_Summ <- print(Fit_cu_Stan, pars=c("theta", "sigma", "z_init"),
                          probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give 3 parms ####

#### Give bcu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real O = theta[2];
    real h = theta[3];

    real dP_dt = P*r - H*(O*P/(1 + O*h*P));
    real dH_dt = 35 + H*(0.3*(O*P/(1 + O*h*P))-0.41);
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
  real<lower = 0> theta[3];   
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
  theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1}]) ~ normal(0.5,0.5);
  log(theta[{2,3}]) ~ normal(0.5,0.5);
  //theta[{4}] ~ normal(35,1);
  //theta[{4}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-bcu.stan")

Model <- stan_model("Global_Stan_Model_Give-bcu.stan")
Fit_bcu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_bcu_Stan_Summ <- print(Fit_bcu_Stan, pars=c("theta", "sigma", "z_init"),
                          probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give hcu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real O = theta[2];
    real b = theta[3];

    real dP_dt = P*r - H*(O*P/(1 + O*0.075*P));
    real dH_dt = b + H*(0.3*(O*P/(1 + O*0.075*P))-0.41);
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
  real<lower = 0> theta[3];   
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
  theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1}]) ~ normal(0.5,0.5);
  log(theta[{2}]) ~ normal(0.5,0.5);
  theta[{3}] ~ normal(35,1);
  //theta[{4}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-hcu.stan")

Model <- stan_model("Global_Stan_Model_Give-hcu.stan")
Fit_hcu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_hcu_Stan_Summ <- print(Fit_hcu_Stan, pars=c("theta", "sigma", "z_init"),
                           probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give hbu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real O = theta[2];
    real c = theta[3];

    real dP_dt = P*r - H*(O*P/(1 + O*0.075*P));
    real dH_dt = 35 + H*(c*(O*P/(1 + O*0.075*P))-0.41);
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
  real<lower = 0> theta[3];   
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
  theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1}]) ~ normal(0.5,0.5);
  log(theta[{2}]) ~ normal(0.5,0.5);
  //theta[{3}] ~ normal(35,1);
  theta[{3}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rOc.stan")

Model <- stan_model("Global_Stan_Model_Give-hbu.stan")
Fit_hbu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_hbu_Stan_Summ <- print(Fit_hbu_Stan, pars=c("theta", "sigma", "z_init"),
                           probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give hbc ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real O = theta[2];
    real u = theta[3];

    real dP_dt = P*r - H*(O*P/(1 + O*0.075*P));
    real dH_dt = 35 + H*(0.3*(O*P/(1 + O*0.075*P))-u);
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
  real<lower = 0> theta[3];   
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
  theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1}]) ~ normal(0.5,0.5);
  log(theta[{2}]) ~ normal(0.5,0.5);
  //theta[{3}] ~ normal(35,1);
  theta[{3}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-hbc.stan")

Model <- stan_model("Global_Stan_Model_Give-hbc.stan")
Fit_hbc_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_hbc_Stan_Summ <- print(Fit_hbc_Stan, pars=c("theta", "sigma", "z_init"),
                           probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give Ocu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real h = theta[2];
    real b = theta[3];

    real dP_dt = P*r - H*(0.012*P/(1 + 0.012*h*P));
    real dH_dt = b + H*(0.3*(0.012*P/(1 + 0.012*h*P))-0.41);
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
  real<lower = 0> theta[3];   
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
  theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1}]) ~ normal(0.5,0.5);
  log(theta[{2}]) ~ normal(0.5,0.5);
  theta[{3}] ~ normal(35,1);
  //theta[{3}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-Ocu.stan")

Model <- stan_model("Global_Stan_Model_Give-Ocu.stan")
Fit_Ocu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_Ocu_Stan_Summ <- print(Fit_Ocu_Stan, pars=c("theta", "sigma", "z_init"),
                           probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give Obu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real h = theta[2];
    real c = theta[3];

    real dP_dt = P*r - H*(0.012*P/(1 + 0.012*h*P));
    real dH_dt = 35 + H*(c*(0.012*P/(1 + 0.012*h*P))-0.41);
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
  real<lower = 0> theta[3];   
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
  theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1}]) ~ normal(0.5,0.5);
  log(theta[{2}]) ~ normal(0.5,0.5);
  //theta[{3}] ~ normal(35,1);
  theta[{3}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-Obu.stan")

Model <- stan_model("Global_Stan_Model_Give-Obu.stan")
Fit_Obu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_Obu_Stan_Summ <- print(Fit_Obu_Stan, pars=c("theta", "sigma", "z_init"),
                           probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give Obc ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real h = theta[2];
    real u = theta[3];

    real dP_dt = P*r - H*(0.012*P/(1 + 0.012*h*P));
    real dH_dt = 35 + H*(0.3*(0.012*P/(1 + 0.012*h*P))-u);
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
  real<lower = 0> theta[3];   
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
  theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1}]) ~ normal(0.5,0.5);
  log(theta[{2}]) ~ normal(0.5,0.5);
  //theta[{3}] ~ normal(35,1);
  theta[{3}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-Obc.stan")

Model <- stan_model("Global_Stan_Model_Give-Obc.stan")
Fit_Obc_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_Obc_Stan_Summ <- print(Fit_Obc_Stan, pars=c("theta", "sigma", "z_init"),
                           probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give Ohu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real b = theta[2];
    real c = theta[3];

    real dP_dt = P*r - H*(0.012*P/(1 + 0.012*0.075*P));
    real dH_dt = b + H*(c*(0.012*P/(1 + 0.012*0.075*P))-0.41);
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
  real<lower = 0> theta[3];   
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
  theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  theta[{2}] ~ normal(35,1);
  theta[{3}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-Ohu.stan")

Model <- stan_model("Global_Stan_Model_Give-Ohu.stan")
Fit_Ohu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_Ohu_Stan_Summ <- print(Fit_Ohu_Stan, pars=c("theta", "sigma", "z_init"),
                           probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give Ohc ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real b = theta[2];
    real u = theta[3];

    real dP_dt = P*r - H*(0.012*P/(1 + 0.012*0.075*P));
    real dH_dt = b + H*(0.3*(0.012*P/(1 + 0.012*0.075*P))-u);
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
  real<lower = 0> theta[3];   
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
  theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  theta[{2}] ~ normal(35,1);
  theta[{3}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-Ohc.stan")

Model <- stan_model("Global_Stan_Model_Give-Ohc.stan")
Fit_Ohc_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_Ohc_Stan_Summ <- print(Fit_Ohc_Stan, pars=c("theta", "sigma", "z_init"),
                           probs=c(0.1, 0.5, 0.9), digits = 3)

Y#### Give Ohb ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real c = theta[2];
    real u = theta[3];

    real dP_dt = P*r - H*(0.012*P/(1 + 0.012*0.075*P));
    real dH_dt = 35 + H*(c*(0.012*P/(1 + 0.012*0.075*P))-u);
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
  real<lower = 0> theta[3];   
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
  theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  //theta[{2}] ~ normal(35,1);
  theta[{2,3}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-Ohb.stan")

Model <- stan_model("Global_Stan_Model_Give-Ohb.stan")
Fit_Ohb_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_Ohb_Stan_Summ <- print(Fit_Ohb_Stan, pars=c("theta", "sigma", "z_init"),
                           probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give rOh ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real b = theta[1];  
    real c = theta[2];
    real u = theta[3];

    real dP_dt = P*2.5 - H*(0.012*P/(1 + 0.012*0.075*P));
    real dH_dt = b + H*(c*(0.012*P/(1 + 0.012*0.075*P))-u);
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
  real<lower = 0> theta[3];   
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
  //theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  theta[{1}] ~ normal(35,1);
  theta[{2,3}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rOh.stan")

Model <- stan_model("Global_Stan_Model_Give-rOh.stan")
Fit_rOh_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rOh_Stan_Summ <- print(Fit_rOh_Stan, pars=c("theta", "sigma", "z_init"),
                           probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give rOb ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real h = theta[1];  
    real c = theta[2];
    real u = theta[3];

    real dP_dt = P*2.5 - H*(0.012*P/(1 + 0.012*h*P));
    real dH_dt = 35 + H*(c*(0.012*P/(1 + 0.012*h*P))-u);
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
  real<lower = 0> theta[3];   
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
  //theta[{1}] ~ normal(2.5, 1);
  log(theta[{1}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  //theta[{1}] ~ normal(35,1);
  theta[{2,3}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rOb.stan")

Model <- stan_model("Global_Stan_Model_Give-rOb.stan")
Fit_rOb_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rOb_Stan_Summ <- print(Fit_rOb_Stan, pars=c("theta", "sigma", "z_init"),
                           probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give rOc ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real h = theta[1];  
    real b = theta[2];
    real u = theta[3];

    real dP_dt = P*2.5 - H*(0.012*P/(1 + 0.012*h*P));
    real dH_dt = b + H*(0.3*(0.012*P/(1 + 0.012*h*P))-u);
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
  real<lower = 0> theta[3];   
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
  //theta[{1}] ~ normal(2.5, 1);
  log(theta[{1}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  theta[{2}] ~ normal(35,1);
  theta[{3}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rOc.stan")

Model <- stan_model("Global_Stan_Model_Give-rOc.stan")
Fit_rOc_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rOc_Stan_Summ <- print(Fit_rOc_Stan, pars=c("theta", "sigma", "z_init"),
                           probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give rOu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real h = theta[1];  
    real b = theta[2];
    real c = theta[3];

    real dP_dt = P*2.5 - H*(0.012*P/(1 + 0.012*h*P));
    real dH_dt = b + H*(c*(0.012*P/(1 + 0.012*h*P))-0.41);
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
  real<lower = 0> theta[3];   
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
  //theta[{1}] ~ normal(2.5, 1);
  log(theta[{1}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  theta[{2}] ~ normal(35,1);
  theta[{3}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rOu.stan")

Model <- stan_model("Global_Stan_Model_Give-rOu.stan")
Fit_rOu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rOu_Stan_Summ <- print(Fit_rOu_Stan, pars=c("theta", "sigma", "z_init"),
                           probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give rhb ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real O = theta[1];  
    real c = theta[2];
    real u = theta[3];

    real dP_dt = P*2.5 - H*(O*P/(1 + O*0.075*P));
    real dH_dt = 35 + H*(c*(O*P/(1 + O*0.075*P))-u);
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
  real<lower = 0> theta[3];   
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
  //theta[{1}] ~ normal(2.5, 1);
  log(theta[{1}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  //theta[{2}] ~ normal(35,1);
  theta[{2,3}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rhb.stan")

Model <- stan_model("Global_Stan_Model_Give-rhb.stan")
Fit_rhb_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rhb_Stan_Summ <- print(Fit_rhb_Stan, pars=c("theta", "sigma", "z_init"),
                           probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give rhc ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real O = theta[1];  
    real b = theta[2];
    real u = theta[3];

    real dP_dt = P*2.5 - H*(O*P/(1 + O*0.075*P));
    real dH_dt = b + H*(0.3*(O*P/(1 + O*0.075*P))-u);
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
  real<lower = 0> theta[3];   
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
  //theta[{1}] ~ normal(2.5, 1);
  log(theta[{1}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  theta[{2}] ~ normal(35,1);
  theta[{3}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rhc.stan")

Model <- stan_model("Global_Stan_Model_Give-rhc.stan")
Fit_rhc_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rhc_Stan_Summ <- print(Fit_rhc_Stan, pars=c("theta", "sigma", "z_init"),
                           probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give rhu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real O = theta[1];  
    real b = theta[2];
    real c = theta[3];

    real dP_dt = P*2.5 - H*(O*P/(1 + O*0.075*P));
    real dH_dt = b + H*(c*(O*P/(1 + O*0.075*P))-0.41);
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
  real<lower = 0> theta[3];   
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
  //theta[{1}] ~ normal(2.5, 1);
  log(theta[{1}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  theta[{2}] ~ normal(35,1);
  theta[{3}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rhu.stan")

Model <- stan_model("Global_Stan_Model_Give-rhu.stan")
Fit_rhu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rhu_Stan_Summ <- print(Fit_rhu_Stan, pars=c("theta", "sigma", "z_init"),
                           probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give rbc ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real O = theta[1];  
    real h = theta[2];
    real u = theta[3];

    real dP_dt = P*2.5 - H*(O*P/(1 + O*h*P));
    real dH_dt = 35 + H*(0.3*(O*P/(1 + O*h*P))-u);
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
  real<lower = 0> theta[3];   
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
  //theta[{1}] ~ normal(2.5, 1);
  log(theta[{1,2}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  //theta[{2}] ~ normal(35,1);
  theta[{3}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rbc.stan")

Model <- stan_model("Global_Stan_Model_Give-rbc.stan")
Fit_rbc_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rbc_Stan_Summ <- print(Fit_rbc_Stan, pars=c("theta", "sigma", "z_init"),
                           probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give rbu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real O = theta[1];  
    real h = theta[2];
    real c = theta[3];

    real dP_dt = P*2.5 - H*(O*P/(1 + O*h*P));
    real dH_dt = 35 + H*(c*(O*P/(1 + O*h*P))-0.41);
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
  real<lower = 0> theta[3];   
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
  //theta[{1}] ~ normal(2.5, 1);
  log(theta[{1,2}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  //theta[{2}] ~ normal(35,1);
  theta[{3}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rbu.stan")

Model <- stan_model("Global_Stan_Model_Give-rbu.stan")
Fit_rbu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rbu_Stan_Summ <- print(Fit_rbu_Stan, pars=c("theta", "sigma", "z_init"),
                           probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give rcu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real O = theta[1];  
    real h = theta[2];
    real b = theta[3];

    real dP_dt = P*2.5 - H*(O*P/(1 + O*h*P));
    real dH_dt = b + H*(0.3*(O*P/(1 + O*h*P))-0.41);
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
  real<lower = 0> theta[3];   
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
  //theta[{1}] ~ normal(2.5, 1);
  log(theta[{1,2}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  theta[{3}] ~ normal(35,1);
  //theta[{3}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rcu.stan")

Model <- stan_model("Global_Stan_Model_Give-rcu.stan")
Fit_rcu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rcu_Stan_Summ <- print(Fit_rcu_Stan, pars=c("theta", "sigma", "z_init"),
                           probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give 4 parms ####

#### Give rOhb ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real c = theta[1];  
    real u = theta[2];

    real dP_dt = P*2.5 - H*(0.012*P/(1 + 0.012*0.075*P));
    real dH_dt = 35 + H*(c*(0.012*P/(1 + 0.012*0.075*P))-u);
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
  real<lower = 0> theta[2];   
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
  //theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  //theta[{3}] ~ normal(35,1);
  theta[{1,2}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rOhb.stan")

Model <- stan_model("Global_Stan_Model_Give-rOhb.stan")
Fit_rOhb_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rOhb_Stan_Summ <- print(Fit_rOhb_Stan, pars=c("theta", "sigma", "z_init"),
                           probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give rOhc ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real b = theta[1];  
    real u = theta[2];

    real dP_dt = P*2.5 - H*(0.012*P/(1 + 0.012*0.075*P));
    real dH_dt = b + H*(0.3*(0.012*P/(1 + 0.012*0.075*P))-u);
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
  real<lower = 0> theta[2];   
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
  //theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  theta[{1}] ~ normal(35,1);
  theta[{2}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rOhc.stan")

Model <- stan_model("Global_Stan_Model_Give-rOhc.stan")
Fit_rOhc_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rOhc_Stan_Summ <- print(Fit_rOhc_Stan, pars=c("theta", "sigma", "z_init"),
                           probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give rOhu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real b = theta[1];  
    real c = theta[2];

    real dP_dt = P*2.5 - H*(0.012*P/(1 + 0.012*0.075*P));
    real dH_dt = b + H*(c*(0.012*P/(1 + 0.012*0.075*P))-0.41);
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
  real<lower = 0> theta[2];   
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
  //theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  theta[{1}] ~ normal(35,1);
  theta[{2}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rOhu.stan")

Model <- stan_model("Global_Stan_Model_Give-rOhu.stan")
Fit_rOhu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rOhu_Stan_Summ <- print(Fit_rOhu_Stan, pars=c("theta", "sigma", "z_init"),
                            probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give rObc ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real h = theta[1];  
    real u = theta[2];

    real dP_dt = P*2.5 - H*(0.012*P/(1 + 0.012*h*P));
    real dH_dt = 35 + H*(0.3*(0.012*P/(1 + 0.012*h*P))-u);
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
  real<lower = 0> theta[2];   
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
  //theta[{1}] ~ normal(2.5, 1);
  log(theta[{1}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  //theta[{1}] ~ normal(35,1);
  theta[{2}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rObc.stan")

Model <- stan_model("Global_Stan_Model_Give-rObc.stan")
Fit_rObc_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rObc_Stan_Summ <- print(Fit_rObc_Stan, pars=c("theta", "sigma", "z_init"),
                            probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give rObu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real h = theta[1];  
    real c = theta[2];

    real dP_dt = P*2.5 - H*(0.012*P/(1 + 0.012*h*P));
    real dH_dt = 35 + H*(c*(0.012*P/(1 + 0.012*h*P))-0.41);
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
  real<lower = 0> theta[2];   
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
  //theta[{1}] ~ normal(2.5, 1);
  log(theta[{1}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  //theta[{1}] ~ normal(35,1);
  theta[{2}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rObu.stan")

Model <- stan_model("Global_Stan_Model_Give-rObu.stan")
Fit_rObu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rObu_Stan_Summ <- print(Fit_rObu_Stan, pars=c("theta", "sigma", "z_init"),
                            probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give rOcu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real h = theta[1];  
    real b = theta[2];

    real dP_dt = P*2.5 - H*(0.012*P/(1 + 0.012*h*P));
    real dH_dt = b + H*(0.3*(0.012*P/(1 + 0.012*h*P))-0.41);
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
  real<lower = 0> theta[2];   
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
  //theta[{1}] ~ normal(2.5, 1);
  log(theta[{1}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  theta[{2}] ~ normal(35,1);
  //theta[{2}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rOcu.stan")

Model <- stan_model("Global_Stan_Model_Give-rOcu.stan")
Fit_rOcu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rOcu_Stan_Summ <- print(Fit_rOcu_Stan, pars=c("theta", "sigma", "z_init"),
                            probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give rhbc ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real O = theta[1];  
    real u = theta[2];

    real dP_dt = P*2.5 - H*(O*P/(1 + O*0.075*P));
    real dH_dt = 35 + H*(0.3*(O*P/(1 + O*0.075*P))-u);
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
  real<lower = 0> theta[2];   
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
  //theta[{1}] ~ normal(2.5, 1);
  log(theta[{1}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  //theta[{2}] ~ normal(35,1);
  theta[{2}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rhbc.stan")

Model <- stan_model("Global_Stan_Model_Give-rhbc.stan")
Fit_rhbc_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rhbc_Stan_Summ <- print(Fit_rhbc_Stan, pars=c("theta", "sigma", "z_init"),
                            probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give rhbu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real O = theta[1];  
    real c = theta[2];

    real dP_dt = P*2.5 - H*(O*P/(1 + O*0.075*P));
    real dH_dt = 35 + H*(c*(O*P/(1 + O*0.075*P))-0.41);
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
  real<lower = 0> theta[2];   
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
  //theta[{1}] ~ normal(2.5, 1);
  log(theta[{1}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  //theta[{2}] ~ normal(35,1);
  theta[{2}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rhbu.stan")

Model <- stan_model("Global_Stan_Model_Give-rhbu.stan")
Fit_rhbu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rhbu_Stan_Summ <- print(Fit_rhbu_Stan, pars=c("theta", "sigma", "z_init"),
                            probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give rhcu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real O = theta[1];  
    real b = theta[2];

    real dP_dt = P*2.5 - H*(O*P/(1 + O*0.075*P));
    real dH_dt = b + H*(0.3*(O*P/(1 + O*0.075*P))-0.41);
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
  real<lower = 0> theta[2];   
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
  //theta[{1}] ~ normal(2.5, 1);
  log(theta[{1}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  theta[{2}] ~ normal(35,1);
  //theta[{3}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rhcu.stan")

Model <- stan_model("Global_Stan_Model_Give-rhcu.stan")
Fit_rhcu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rhcu_Stan_Summ <- print(Fit_rhcu_Stan, pars=c("theta", "sigma", "z_init"),
                            probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give rbcu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real O = theta[1];  
    real h = theta[2];

    real dP_dt = P*2.5 - H*(O*P/(1 + O*h*P));
    real dH_dt = 35 + H*(0.3*(O*P/(1 + O*h*P))-0.41);
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
  real<lower = 0> theta[2];   
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
  //theta[{1}] ~ normal(2.5, 1);
  log(theta[{1,2}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  //theta[{2}] ~ normal(35,1);
  //theta[{2}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rbcu.stan")

Model <- stan_model("Global_Stan_Model_Give-rbcu.stan")
Fit_rbcu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rbcu_Stan_Summ <- print(Fit_rbcu_Stan, pars=c("theta", "sigma", "z_init"),
                            probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give Ohbc ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real u = theta[2];

    real dP_dt = P*r - H*(0.012*P/(1 + 0.012*0.075*P));
    real dH_dt = 35 + H*(0.3*(0.012*P/(1 + 0.012*0.075*P))-u);
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
  real<lower = 0> theta[2];   
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
  theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  //theta[{2}] ~ normal(35,1);
  theta[{2}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-Ohbc.stan")

Model <- stan_model("Global_Stan_Model_Give-Ohbc.stan")
Fit_Ohbc_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_Ohbc_Stan_Summ <- print(Fit_Ohbc_Stan, pars=c("theta", "sigma", "z_init"),
                            probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give Ohbu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real c = theta[2];

    real dP_dt = P*r - H*(0.012*P/(1 + 0.012*0.075*P));
    real dH_dt = 35 + H*(c*(0.012*P/(1 + 0.012*0.075*P))-0.41);
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
  real<lower = 0> theta[2];   
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
  theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1,2}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  //theta[{2}] ~ normal(35,1);
  theta[{2}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-Ohbu.stan")

Model <- stan_model("Global_Stan_Model_Give-Ohbu.stan")
Fit_Ohbu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_Ohbu_Stan_Summ <- print(Fit_Ohbu_Stan, pars=c("theta", "sigma", "z_init"),
                            probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give Ohcu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real b = theta[2];

    real dP_dt = P*r - H*(0.012*P/(1 + 0.012*0.075*P));
    real dH_dt = b + H*(0.3*(0.012*P/(1 + 0.012*0.075*P))-0.41);
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
  real<lower = 0> theta[2];   
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
  theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1,2}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  theta[{2}] ~ normal(35,1);
  //theta[{2}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-Ohcu.stan")

Model <- stan_model("Global_Stan_Model_Give-Ohcu.stan")
Fit_Ohcu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_Ohcu_Stan_Summ <- print(Fit_Ohcu_Stan, pars=c("theta", "sigma", "z_init"),
                            probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give Obcu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real h = theta[2];

    real dP_dt = P*r - H*(0.012*P/(1 + 0.012*h*P));
    real dH_dt = 35 + H*(0.3*(0.012*P/(1 + 0.012*h*P))-0.41);
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
  real<lower = 0> theta[2];   
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
  theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1,2}]) ~ normal(0.5,0.5);
  log(theta[{2}]) ~ normal(0.5,0.5);
  //theta[{2}] ~ normal(35,1);
  //theta[{2}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-Obcu.stan")

Model <- stan_model("Global_Stan_Model_Give-Obcu.stan")
Fit_Obcu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_Obcu_Stan_Summ <- print(Fit_Obcu_Stan, pars=c("theta", "sigma", "z_init"),
                            probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give hbcu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  
    real O = theta[2];

    real dP_dt = P*r - H*(O*P/(1 + O*0.075*P));
    real dH_dt = 35 + H*(0.3*(O*P/(1 + O*0.075*P))-0.41);
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
  real<lower = 0> theta[2];   
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
  theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1,2}]) ~ normal(0.5,0.5);
  log(theta[{2}]) ~ normal(0.5,0.5);
  //theta[{2}] ~ normal(35,1);
  //theta[{2}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-hbcu.stan")

Model <- stan_model("Global_Stan_Model_Give-hbcu.stan")
Fit_hbcu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_hbcu_Stan_Summ <- print(Fit_hbcu_Stan, pars=c("theta", "sigma", "z_init"),
                            probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give 5 parms ####

#### Give rOhbc ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real u = theta[1];  

    real dP_dt = P*2.5 - H*(0.012*P/(1 + 0.012*0.075*P));
    real dH_dt = 35 + H*(0.3*(0.012*P/(1 + 0.012*0.075*P))-u);
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
  real<lower = 0> theta[1];   
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
  //theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1,2}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  //theta[{2}] ~ normal(35,1);
  theta[{1}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rOhbc.stan")

Model <- stan_model("Global_Stan_Model_Give-rOhbc.stan")
Fit_rOhbc_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rOhbc_Stan_Summ <- print(Fit_rOhbc_Stan, pars=c("theta", "sigma", "z_init"),
                            probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give rOhbu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real c = theta[1];  

    real dP_dt = P*2.5 - H*(0.012*P/(1 + 0.012*0.075*P));
    real dH_dt = 35 + H*(c*(0.012*P/(1 + 0.012*0.075*P))-0.41);
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
  real<lower = 0> theta[1];   
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
  //theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1,2}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  //theta[{2}] ~ normal(35,1);
  theta[{1}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rOhbu.stan")

Model <- stan_model("Global_Stan_Model_Give-rOhbu.stan")
Fit_rOhbu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rOhbu_Stan_Summ <- print(Fit_rOhbu_Stan, pars=c("theta", "sigma", "z_init"),
                             probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give rOhcu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real b = theta[1];  

    real dP_dt = P*2.5 - H*(0.012*P/(1 + 0.012*0.075*P));
    real dH_dt = b + H*(0.3*(0.012*P/(1 + 0.012*0.075*P))-0.41);
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
  real<lower = 0> theta[1];   
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
  //theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1,2}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  theta[{1}] ~ normal(35,1);
  //theta[{1}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rOhcu.stan")

Model <- stan_model("Global_Stan_Model_Give-rOhcu.stan")
Fit_rOhcu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rOhcu_Stan_Summ <- print(Fit_rOhcu_Stan, pars=c("theta", "sigma", "z_init"),
                             probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give rObcu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real h = theta[1];  

    real dP_dt = P*2.5 - H*(0.012*P/(1 + 0.012*h*P));
    real dH_dt = 35 + H*(0.3*(0.012*P/(1 + 0.012*h*P))-0.41);
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
  real<lower = 0> theta[1];   
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
  //theta[{1}] ~ normal(2.5, 1);
  log(theta[{1}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  //theta[{1}] ~ normal(35,1);
  //theta[{1}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rObcu.stan")

Model <- stan_model("Global_Stan_Model_Give-rObcu.stan")
Fit_rObcu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rObcu_Stan_Summ <- print(Fit_rObcu_Stan, pars=c("theta", "sigma", "z_init"),
                             probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give rhbcu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real O = theta[1];  

    real dP_dt = P*2.5 - H*(O*P/(1 + O*0.075*P));
    real dH_dt = 35 + H*(0.3*(O*P/(1 + O*0.075*P))-0.41);
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
  real<lower = 0> theta[1];   
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
  //theta[{1}] ~ normal(2.5, 1);
  log(theta[{1}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  //theta[{1}] ~ normal(35,1);
  //theta[{1}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-rhbcu.stan")

Model <- stan_model("Global_Stan_Model_Give-rhbcu.stan")
Fit_rhbcu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_rhbcu_Stan_Summ <- print(Fit_rhbcu_Stan, pars=c("theta", "sigma", "z_init"),
                             probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give Ohbcu ####
write("
functions {
  real[] dz_dt(real t,       
               real[] z,     
               real[] theta, 
               real[] x_r,  
               int[] x_i) {
    real P = z[1];
    real H = z[2];

    real r = theta[1];  

    real dP_dt = P*r - H*(0.012*P/(1 + 0.012*0.075*P));
    real dH_dt = 35 + H*(0.3*(0.012*P/(1 + 0.012*0.075*P))-0.41);
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
  real<lower = 0> theta[1];   
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
  theta[{1}] ~ normal(2.5, 1);
  //log(theta[{1}]) ~ normal(0.5,0.5);
  //log(theta[{2}]) ~ normal(0.5,0.5);
  //theta[{1}] ~ normal(35,1);
  //theta[{1}] ~ normal(0.5,0.5);
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
}",
"Global_Stan_Model_Give-Ohbcu.stan")

Model <- stan_model("Global_Stan_Model_Give-Ohbcu.stan")
Fit_Ohbcu_Stan <- sampling(Model, data = Data, chains = 4, iter = 10000, warmup = 5000, thin = 1, cores = 4, seed = 12); beep(3)
Fit_Ohbcu_Stan_Summ <- print(Fit_Ohbcu_Stan, pars=c("theta", "sigma", "z_init"),
                             probs=c(0.1, 0.5, 0.9), digits = 3)

#### Make outputs in .csv ####
Fit_Stan_Summ <- data.frame(MCMCsummary(Fit_Stan))
write.csv(Fit_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Fit_Stan.csv",row.names=TRUE)

####

Fit_r_Stan_Summ <- data.frame(summary(Fit_r_Stan))
write.csv(Fit_r_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_r_Stan.csv",row.names=TRUE)

Fit_O_Stan_Summ <- data.frame(summary(Fit_O_Stan))
write.csv(Fit_O_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_O_Stan.csv",row.names=TRUE)

Fit_h_Stan_Summ <- data.frame(summary(Fit_h_Stan))
write.csv(Fit_h_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_h_Stan.csv",row.names=TRUE)

Fit_b_Stan_Summ <- data.frame(summary(Fit_b_Stan))
write.csv(Fit_b_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_b_Stan.csv",row.names=TRUE)

Fit_c_Stan_Summ <- data.frame(summary(Fit_c_Stan))
write.csv(Fit_c_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_c_Stan.csv",row.names=TRUE)

Fit_u_Stan_Summ <- data.frame(summary(Fit_u_Stan))
write.csv(Fit_u_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_u_Stan.csv",row.names=TRUE)

Fit_rO_Stan_Summ <- data.frame(summary(Fit_rO_Stan))
write.csv(Fit_rO_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rO_Stan.csv",row.names=TRUE)

Fit_rh_Stan_Summ <- data.frame(summary(Fit_rh_Stan))
write.csv(Fit_rh_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rh_Stan.csv",row.names=TRUE)

Fit_rb_Stan_Summ <- data.frame(summary(Fit_rb_Stan))
write.csv(Fit_rb_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rb_Stan.csv",row.names=TRUE)

Fit_rc_Stan_Summ <- data.frame(summary(Fit_rc_Stan))
write.csv(Fit_rc_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rc_Stan.csv",row.names=TRUE)

Fit_ru_Stan_Summ <- data.frame(summary(Fit_ru_Stan))
write.csv(Fit_ru_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_ru_Stan.csv",row.names=TRUE)

Fit_Oh_Stan_Summ <- data.frame(summary(Fit_Oh_Stan))
write.csv(Fit_Oh_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_Oh_Stan.csv",row.names=TRUE)

Fit_Ob_Stan_Summ <- data.frame(summary(Fit_Ob_Stan))
write.csv(Fit_Ob_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_Ob_Stan.csv",row.names=TRUE)

Fit_Oc_Stan_Summ <- data.frame(summary(Fit_Oc_Stan))
write.csv(Fit_Oc_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_Oc_Stan.csv",row.names=TRUE)

Fit_Ou_Stan_Summ <- data.frame(summary(Fit_Ou_Stan))
write.csv(Fit_Ou_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_Ou_Stan.csv",row.names=TRUE)

Fit_hb_Stan_Summ <- data.frame(summary(Fit_hb_Stan))
write.csv(Fit_hb_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_hb_Stan.csv",row.names=TRUE)

Fit_hc_Stan_Summ <- data.frame(summary(Fit_hc_Stan))
write.csv(Fit_hc_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_hc_Stan.csv",row.names=TRUE)

Fit_hu_Stan_Summ <- data.frame(summary(Fit_hu_Stan))
write.csv(Fit_hu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_hu_Stan.csv",row.names=TRUE)

Fit_bc_Stan_Summ <- data.frame(summary(Fit_bc_Stan))
write.csv(Fit_bc_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_bc_Stan.csv",row.names=TRUE)

Fit_bu_Stan_Summ <- data.frame(summary(Fit_bu_Stan))
write.csv(Fit_bu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_bu_Stan.csv",row.names=TRUE)

Fit_cu_Stan_Summ <- data.frame(summary(Fit_cu_Stan))
write.csv(Fit_cu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_cu_Stan.csv",row.names=TRUE)

Fit_rOh_Stan_Summ <- data.frame(summary(Fit_rOh_Stan))
write.csv(Fit_rOh_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rOh_Stan.csv",row.names=TRUE)

Fit_rOb_Stan_Summ <- data.frame(summary(Fit_rOb_Stan))
write.csv(Fit_rOb_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rOb_Stan.csv",row.names=TRUE)

Fit_rOc_Stan_Summ <- data.frame(summary(Fit_rOc_Stan))
write.csv(Fit_rOc_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rOc_Stan.csv",row.names=TRUE)

Fit_rOu_Stan_Summ <- data.frame(summary(Fit_rOu_Stan))
write.csv(Fit_rOu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rOu_Stan.csv",row.names=TRUE)

Fit_rhb_Stan_Summ <- data.frame(summary(Fit_rhb_Stan))
write.csv(Fit_rhb_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rhb_Stan.csv",row.names=TRUE)

Fit_rhc_Stan_Summ <- data.frame(summary(Fit_rhc_Stan))
write.csv(Fit_rhc_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rhc_Stan.csv",row.names=TRUE)

Fit_rhu_Stan_Summ <- data.frame(summary(Fit_rhu_Stan))
write.csv(Fit_rhu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rhu_Stan.csv",row.names=TRUE)

Fit_rbc_Stan_Summ <- data.frame(summary(Fit_rbc_Stan))
write.csv(Fit_rbc_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rbc_Stan.csv",row.names=TRUE)

Fit_rbu_Stan_Summ <- data.frame(summary(Fit_rbu_Stan))
write.csv(Fit_rbu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rbu_Stan.csv",row.names=TRUE)

Fit_rcu_Stan_Summ <- data.frame(summary(Fit_rcu_Stan))
write.csv(Fit_rcu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rcu_Stan.csv",row.names=TRUE)

Fit_Ohb_Stan_Summ <- data.frame(summary(Fit_Ohb_Stan))
write.csv(Fit_Ohb_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_Ohb_Stan.csv",row.names=TRUE)

Fit_Ohc_Stan_Summ <- data.frame(summary(Fit_Ohc_Stan))
write.csv(Fit_Ohc_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_Ohc_Stan.csv",row.names=TRUE)

Fit_Ohu_Stan_Summ <- data.frame(summary(Fit_Ohu_Stan))
write.csv(Fit_Ohu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_Ohu_Stan.csv",row.names=TRUE)

Fit_Obc_Stan_Summ <- data.frame(summary(Fit_Obc_Stan))
write.csv(Fit_Obc_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_Obc_Stan.csv",row.names=TRUE)

Fit_Obu_Stan_Summ <- data.frame(summary(Fit_Obu_Stan))
write.csv(Fit_Obu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_Obu_Stan.csv",row.names=TRUE)

Fit_Ocu_Stan_Summ <- data.frame(summary(Fit_Ocu_Stan))
write.csv(Fit_Ocu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_Ocu_Stan.csv",row.names=TRUE)

Fit_hbc_Stan_Summ <- data.frame(summary(Fit_hbc_Stan))
write.csv(Fit_hbc_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_hbc_Stan.csv",row.names=TRUE)

Fit_hbu_Stan_Summ <- data.frame(summary(Fit_hbu_Stan))
write.csv(Fit_hbu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_hbu_Stan.csv",row.names=TRUE)

Fit_hcu_Stan_Summ <- data.frame(summary(Fit_hcu_Stan))
write.csv(Fit_hcu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_hcu_Stan.csv",row.names=TRUE)

Fit_bcu_Stan_Summ <- data.frame(summary(Fit_bcu_Stan))
write.csv(Fit_bcu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_bcu_Stan.csv",row.names=TRUE)

Fit_rOhb_Stan_Summ <- data.frame(summary(Fit_rOhb_Stan))
write.csv(Fit_rOhb_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rOhb_Stan.csv",row.names=TRUE)

Fit_rOhc_Stan_Summ <- data.frame(summary(Fit_rOhc_Stan))
write.csv(Fit_rOhc_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rOhc_Stan.csv",row.names=TRUE)

Fit_rOhu_Stan_Summ <- data.frame(summary(Fit_rOhu_Stan))
write.csv(Fit_rOhu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rOhu_Stan.csv",row.names=TRUE)

Fit_rObc_Stan_Summ <- data.frame(summary(Fit_rObc_Stan))
write.csv(Fit_rObc_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rObc_Stan.csv",row.names=TRUE)

Fit_rObu_Stan_Summ <- data.frame(summary(Fit_rObu_Stan))
write.csv(Fit_rObu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rObu_Stan.csv",row.names=TRUE)

Fit_rOcu_Stan_Summ <- data.frame(summary(Fit_rOcu_Stan))
write.csv(Fit_rOcu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rOcu_Stan.csv",row.names=TRUE)

Fit_rhbc_Stan_Summ <- data.frame(summary(Fit_rhbc_Stan))
write.csv(Fit_rhbc_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rhbc_Stan.csv",row.names=TRUE)

Fit_rhbu_Stan_Summ <- data.frame(summary(Fit_rhbu_Stan))
write.csv(Fit_rhbu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rhbu_Stan.csv",row.names=TRUE)

Fit_rhcu_Stan_Summ <- data.frame(summary(Fit_rhcu_Stan))
write.csv(Fit_rhcu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rhcu_Stan.csv",row.names=TRUE)

Fit_rbcu_Stan_Summ <- data.frame(summary(Fit_rbcu_Stan))
write.csv(Fit_rbcu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rbcu_Stan.csv",row.names=TRUE)

Fit_Ohbc_Stan_Summ <- data.frame(summary(Fit_Ohbc_Stan))
write.csv(Fit_Ohbc_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_Ohbc_Stan.csv",row.names=TRUE)

Fit_Ohbu_Stan_Summ <- data.frame(summary(Fit_Ohbu_Stan))
write.csv(Fit_Ohbu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_Ohbu_Stan.csv",row.names=TRUE)

Fit_Ohcu_Stan_Summ <- data.frame(summary(Fit_Ohcu_Stan))
write.csv(Fit_Ohcu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_Ohcu_Stan.csv",row.names=TRUE)

Fit_Obcu_Stan_Summ <- data.frame(summary(Fit_Obcu_Stan))
write.csv(Fit_Obcu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_Obcu_Stan.csv",row.names=TRUE)

Fit_hbcu_Stan_Summ <- data.frame(summary(Fit_hbcu_Stan))
write.csv(Fit_hbcu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_hbcu_Stan.csv",row.names=TRUE)

Fit_rOhbc_Stan_Summ <- data.frame(summary(Fit_rOhbc_Stan))
write.csv(Fit_rOhbc_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rOhbc_Stan.csv",row.names=TRUE)

Fit_rOhbu_Stan_Summ <- data.frame(summary(Fit_rOhbu_Stan))
write.csv(Fit_rOhbu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rOhbu_Stan.csv",row.names=TRUE)

Fit_rOhcu_Stan_Summ <- data.frame(summary(Fit_rOhcu_Stan))
write.csv(Fit_rOhcu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rOhcu_Stan.csv",row.names=TRUE)

Fit_rObcu_Stan_Summ <- data.frame(summary(Fit_rObcu_Stan))
write.csv(Fit_rObcu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rObcu_Stan.csv",row.names=TRUE)

Fit_rhbcu_Stan_Summ <- data.frame(summary(Fit_rhbcu_Stan))
write.csv(Fit_rhbcu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_rhbcu_Stan.csv",row.names=TRUE)

Fit_Ohbcu_Stan_Summ <- data.frame(summary(Fit_Ohbcu_Stan))
write.csv(Fit_Ohbcu_Stan_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/OSA_JAGS_AddInfoOutputs/Fit_Ohbcu_Stan.csv",row.names=TRUE)

#### PDF Outputs ####
# To specify which parms; parms <- c("theta[1]","theta[2]","theta[3]","theta[4]","theta[5]")

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_Stan.pdf",onefile=TRUE)
stan_trace(Fit_Stan)
dev.off()

#### 1 ####

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_r_Stan.pdf",onefile=TRUE)
stan_trace(Fit_r_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_O_Stan.pdf",onefile=TRUE)
stan_trace(Fit_O_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_h_Stan.pdf",onefile=TRUE)
stan_trace(Fit_h_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_b_Stan.pdf",onefile=TRUE)
stan_trace(Fit_b_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_c_Stan.pdf",onefile=TRUE)
stan_trace(Fit_c_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_u_Stan.pdf",onefile=TRUE)
stan_trace(Fit_u_Stan)
dev.off()

#### 2 ####

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rO_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rO_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rh_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rh_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rb_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rb_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rc_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rc_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_ru_Stan.pdf",onefile=TRUE)
stan_trace(Fit_ru_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_Oh_Stan.pdf",onefile=TRUE)
stan_trace(Fit_Oh_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_Ob_Stan.pdf",onefile=TRUE)
stan_trace(Fit_Ob_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_Oc_Stan.pdf",onefile=TRUE)
stan_trace(Fit_Oc_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_Ou_Stan.pdf",onefile=TRUE)
stan_trace(Fit_Ou_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_hb_Stan.pdf",onefile=TRUE)
stan_trace(Fit_hb_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_hc_Stan.pdf",onefile=TRUE)
stan_trace(Fit_hc_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_hu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_hu_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_bc_Stan.pdf",onefile=TRUE)
stan_trace(Fit_bc_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_bu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_bu_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_cu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_cu_Stan)
dev.off()

#### 3 #### 

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rOh_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rOh_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rOb_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rOb_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rOc_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rOc_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rOu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rOu_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rhb_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rhb_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rhc_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rhc_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rhu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rhu_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rbc_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rbc_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rbu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rbu_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rcu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rcu_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_Ohb_Stan.pdf",onefile=TRUE)
stan_trace(Fit_Ohb_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_Ohc_Stan.pdf",onefile=TRUE)
stan_trace(Fit_Ohc_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_Ohu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_Ohu_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_Obc_Stan.pdf",onefile=TRUE)
stan_trace(Fit_Obc_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_Obu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_Obu_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_Ocu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_Ocu_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_hbc_Stan.pdf",onefile=TRUE)
stan_trace(Fit_hbc_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_hbu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_hbu_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_hcu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_hcu_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_bcu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_bcu_Stan)
dev.off()

#### 4 ####

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rOhb_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rOhb_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rOhc_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rOhc_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rOhu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rOhu_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rObc_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rObc_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rObu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rObu_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rOcu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rOcu_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rhbc_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rhbc_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rhbu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rhbu_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rhcu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rhcu_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rbcu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rbcu_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_Ohbc_Stan.pdf",onefile=TRUE)
stan_trace(Fit_Ohbc_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_Ohbu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_Ohbu_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_Ohcu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_Ohcu_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_Obcu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_Obcu_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_hbcu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_hbcu_Stan)
dev.off()

#### 5 ####

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rOhbcu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rOhbc_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rOhbu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rOhbu_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rOhcu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rOhcu_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rObcu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rObcu_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_rhbcu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_rhbcu_Stan)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Global_Stan_AddInfoOutputs/Plots/Fit_Ohbcu_Stan.pdf",onefile=TRUE)
stan_trace(Fit_Ohbcu_Stan)
dev.off()