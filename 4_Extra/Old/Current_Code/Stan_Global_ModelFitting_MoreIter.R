library(rstan)
library(gdata)
library(bayesplot)
library(tidyverse)
library(brms)
library(beepr)

## Data org.
NewData = read_csv('/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Data/NewData.csv')
N <- length(NewData$t) - 1
ts <- 1:N
y_init <- c(NewData$P[1], NewData$H[1])
y <- as.matrix(NewData[2:(N + 1), 2:3])
y <- cbind(y[ , 1], y[ , 2]); 
Data <- list(N = N, ts = ts, y_init = y_init, y = y)

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
Fit_r_Stan <- sampling(Model, data = Data, 
                       chains = 4, iter = 10000, warmup = 5000, thin = 1,
                       cores = 2, seed = 123); beep(2)
Fit_r_Stan_Summ <- print(Fit_r_Stan, pars=c("theta", "sigma", "z_init"),
                         probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give b ####
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
Fit_b_Stan <- sampling(Model, data = Data, 
                       chains = 4, iter = 10000, warmup = 5000, thin = 1,
                       cores = 2, seed = 123); beep(2)
Fit_b_Stan_Summ <- print(Fit_b_Stan, pars=c("theta", "sigma", "z_init"),
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
Fit_u_Stan <- sampling(Model, data = Data, 
                       chains = 4, iter = 10000, warmup = 5000, thin = 1,
                       cores = 2, seed = 123); beep(2)
Fit_u_Stan_Summ <- print(Fit_u_Stan, pars=c("theta", "sigma", "z_init"),
                         probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give 2 parms ####

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
Fit_rb_Stan <- sampling(Model, data = Data, 
                        chains = 4, iter = 10000, warmup = 5000, thin = 1,
                        cores = 2, seed = 123); beep(2)
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
Fit_rc_Stan <- sampling(Model, data = Data, 
                        chains = 4, iter = 10000, warmup = 5000, thin = 1,
                        cores = 2, seed = 123); beep(2)
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
Fit_ru_Stan <- sampling(Model, data = Data, 
                        chains = 4, iter = 10000, warmup = 5000, thin = 1,
                        cores = 2, seed = 123); beep(2)
Fit_ru_Stan_Summ <- print(Fit_ru_Stan, pars=c("theta", "sigma", "z_init"),
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
Fit_Ou_Stan <- sampling(Model, data = Data, 
                        chains = 4, iter = 10000, warmup = 5000, thin = 1,
                        cores = 2, seed = 123); beep(2)
Fit_Ou_Stan_Summ <- print(Fit_Ou_Stan, pars=c("theta", "sigma", "z_init"),
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
Fit_hu_Stan <- sampling(Model, data = Data, 
                        chains = 4, iter = 10000, warmup = 5000, thin = 1,
                        cores = 2, seed = 123); beep(2)
Fit_hu_Stan_Summ <- print(Fit_hu_Stan, pars=c("theta", "sigma", "z_init"),
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
Fit_bu_Stan <- sampling(Model, data = Data, 
                        chains = 4, iter = 10000, warmup = 5000, thin = 1,
                        cores = 2, seed = 123); beep(2)
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
Fit_cu_Stan <- sampling(Model, data = Data, 
                        chains = 4, iter = 10000, warmup = 5000, thin = 1,
                        cores = 2, seed = 123); beep(2)
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
Fit_bcu_Stan <- sampling(Model, data = Data, 
                         chains = 4, iter = 10000, warmup = 5000, thin = 1,
                         cores = 2, seed = 123); beep(2)
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
Fit_hcu_Stan <- sampling(Model, data = Data, 
                         chains = 4, iter = 10000, warmup = 5000, thin = 1,
                         cores = 2, seed = 123); beep(2)
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
Fit_hbu_Stan <- sampling(Model, data = Data, 
                         chains = 4, iter = 10000, warmup = 5000, thin = 1,
                         cores = 2, seed = 123); beep(2)
Fit_hbu_Stan_Summ <- print(Fit_hbu_Stan, pars=c("theta", "sigma", "z_init"),
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
Fit_Ohu_Stan <- sampling(Model, data = Data, 
                         chains = 4, iter = 10000, warmup = 5000, thin = 1,
                         cores = 2, seed = 123); beep(2)
Fit_Ohu_Stan_Summ <- print(Fit_Ohu_Stan, pars=c("theta", "sigma", "z_init"),
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
Fit_rOh_Stan <- sampling(Model, data = Data, 
                         chains = 4, iter = 10000, warmup = 5000, thin = 1,
                         cores = 2, seed = 123); beep(2)
Fit_rOh_Stan_Summ <- print(Fit_rOh_Stan, pars=c("theta", "sigma", "z_init"),
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
Fit_rhc_Stan <- sampling(Model, data = Data, 
                         chains = 4, iter = 10000, warmup = 5000, thin = 1,
                         cores = 2, seed = 123); beep(2)
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
Fit_rhu_Stan <- sampling(Model, data = Data, 
                         chains = 4, iter = 10000, warmup = 5000, thin = 1,
                         cores = 2, seed = 123); beep(2)
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
Fit_rbc_Stan <- sampling(Model, data = Data, 
                         chains = 4, iter = 10000, warmup = 5000, thin = 1,
                         cores = 2, seed = 123); beep(2)
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
Fit_rbu_Stan <- sampling(Model, data = Data, 
                         chains = 4, iter = 10000, warmup = 5000, thin = 1,
                         cores = 2, seed = 123); beep(2)
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
Fit_rcu_Stan <- sampling(Model, data = Data, 
                         chains = 4, iter = 10000, warmup = 5000, thin = 1,
                         cores = 2, seed = 123); beep(2)
Fit_rcu_Stan_Summ <- print(Fit_rcu_Stan, pars=c("theta", "sigma", "z_init"),
                           probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give 4 parms ####

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
Fit_rOhc_Stan <- sampling(Model, data = Data, 
                          chains = 4, iter = 10000, warmup = 5000, thin = 1,
                          cores = 2, seed = 123); beep(2)
Fit_rOhc_Stan_Summ <- print(Fit_rOhc_Stan, pars=c("theta", "sigma", "z_init"),
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
Fit_rhbu_Stan <- sampling(Model, data = Data, 
                          chains = 4, iter = 10000, warmup = 5000, thin = 1,
                          cores = 2, seed = 123); beep(2)
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
Fit_rhcu_Stan <- sampling(Model, data = Data, 
                          chains = 4, iter = 10000, warmup = 5000, thin = 1,
                          cores = 2, seed = 123); beep(2)
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
Fit_rbcu_Stan <- sampling(Model, data = Data, 
                          chains = 4, iter = 10000, warmup = 5000, thin = 1,
                          cores = 2, seed = 123); beep(2)
Fit_rbcu_Stan_Summ <- print(Fit_rbcu_Stan, pars=c("theta", "sigma", "z_init"),
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
Fit_Ohcu_Stan <- sampling(Model, data = Data, 
                          chains = 4, iter = 10000, warmup = 5000, thin = 1,
                          cores = 2, seed = 123); beep(2)
Fit_Ohcu_Stan_Summ <- print(Fit_Ohcu_Stan, pars=c("theta", "sigma", "z_init"),
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
Fit_hbcu_Stan <- sampling(Model, data = Data, 
                          chains = 4, iter = 10000, warmup = 5000, thin = 1,
                          cores = 2, seed = 123); beep(2)
Fit_hbcu_Stan_Summ <- print(Fit_hbcu_Stan, pars=c("theta", "sigma", "z_init"),
                            probs=c(0.1, 0.5, 0.9), digits = 3)

#### Give 5 parms ####

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
Fit_rOhbu_Stan <- sampling(Model, data = Data, 
                           chains = 4, iter = 10000, warmup = 5000, thin = 1,
                           cores = 2, seed = 123); beep(2)
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
Fit_rOhcu_Stan <- sampling(Model, data = Data, 
                           chains = 4, iter = 10000, warmup = 5000, thin = 1,
                           cores = 2, seed = 123); beep(2)
Fit_rOhcu_Stan_Summ <- print(Fit_rOhcu_Stan, pars=c("theta", "sigma", "z_init"),
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
Fit_rhbcu_Stan <- sampling(Model, data = Data, 
                           chains = 4, iter = 10000, warmup = 5000, thin = 1,
                           cores = 2, seed = 123); beep(2)
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
Fit_Ohbcu_Stan <- sampling(Model, data = Data, 
                           chains = 4, iter = 10000, warmup = 5000, thin = 1,
                           cores = 2, seed = 123); beep(2)
Fit_Ohbcu_Stan_Summ <- print(Fit_Ohbcu_Stan, pars=c("theta", "sigma", "z_init"),
                             probs=c(0.1, 0.5, 0.9), digits = 3)
