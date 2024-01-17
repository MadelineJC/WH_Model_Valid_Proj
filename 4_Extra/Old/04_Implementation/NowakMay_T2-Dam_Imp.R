library(rstan)
library(deSolve)
library(GillespieSSA)
library(gdata)
library(bayesplot)
library(tidyverse)
library(brms)
library(beepr)
library(here)

#### Estimate all parms ===================================
## Wrangle data ##
StanData <- read.csv("02_Data/NowakMay_StochSim_T2-Dam.csv")
N <- length(StanData$Time) - 1
ts <- 1:N
y_init <- c(StanData$Parasites[1], StanData$ImmuneCells[1])
y <- as.matrix(StanData[2:(N + 1), 2:3])
y <- cbind(y[ , 1], y[ , 2]); 
StanData <- list(N = N, ts = ts, y_init = y_init, y = y)

## Compile model ##
model <- stan_model("03_Fit/NowakMay_T2-Dam_Fit.stan")

# Fitting
# inits <- list(
#   r = 2.5,
#   O = 0.008,
#   h = 0.06,
#   b = 35,
#   c = 0.2,
#   u = 0.2
# ); list_inits <- list(inits, inits, inits, inits)
Fit_NM_T2_Dam <- sampling(model, data = StanData, chains = 2, iter = 1000, cores = 2, seed = 1)
system("say The run is done!")

## Summarizing the fit
Fit_NM_T2_Dam_Summ <- print(Fit_NM_T2_Dam, pars=c("r", "O", "h", "b", "c", "u", "sigma", "z_init"),
                            probs=c(0.1, 0.5, 0.9), digits = 3)

## Model checks ##
Fit_NM_T2_Dam_Array <- as.array(Fit_NM_T2_Dam)[,,-4] 
# dimnames(Fit_NM_T2_Dam_Array)[[2]] <- c("Chain 1", "Chain 2", "Chain 3", "Chain 4")
dimnames(Fit_NM_T2_Dam_Array)[[2]] <- c("Chain 1", "Chain 2")
dimnames(Fit_NM_T2_Dam_Array)[[3]][1:6] <- c("r = 2.5", "O = 0.008", "h = 0.06", "u = 0.2", "c = 0.2", "b = 35")
parms <- c("r = 2.5", "O = 0.008", "h = 0.06", "u = 0.2", "c = 0.2", "b = 35")
color_scheme_set("viridisA")
mcmc_trace(Fit_NM_T2_Dam_Array, parms)
mcmc_pairs(Fit_NM_T2_Dam_Array, parms)
mcmc_dens(Fit_NM_T2_Dam_Array, parms)

#### Give O and h =========================================
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
    // real O = theta[2];
    // real h = theta[3];
    real b = theta[2];
    real c = theta[3];
    real u = theta[4];
    
    real dP_dt = P*r - H*(0.008*P/(1 + 0.008*0.06*P));
    real dH_dt = b + H*(c*(0.008*P/(1 + 0.008*0.06*P))-u);
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
  real<lower = 0> r;
  // real<lower = 0> O;
  // real<lower = 0> h;
  real<lower = 0> b;
  real<lower = 0> c;
  real<lower = 0> u;
  real<lower = 0> z_init[2];  
  real<lower = 0> sigma[2];   
}
transformed parameters {
  real z[N, 2]
  = integrate_ode_bdf(dz_dt, z_init, 0, ts, {r, b, c, u},
                       rep_array(0.0, 0), rep_array(0, 0));
}
model {
  r ~ normal(1, 3); // r = 2.5
  // O ~ normal(0, 1); // O = 0.008
  // h ~ normal(0, 1); // h = 0.06
  b ~ uniform(0, 100); // b = 35
  c ~ normal(0, 1); // c = 0.2
  u ~ normal(0, 1); // u = 0.2
  sigma ~ lognormal(-1, 1);
  z_init[1] ~ normal(80, 1); // 80
  z_init[2] ~ normal(200, 1); // 200
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
"03_Fit/NowakMay_T2-Dam_Fit_Oh.stan")

## Compile model ##
model <- stan_model("03_Fit/NowakMay_T2-Dam_Fit_Oh.stan")

Fit_NM_T2_Dam <- sampling(model, data = StanData, chains = 1, iter = 1000, cores = 2, 
                          control = list(stepsize = 0.001),
                                         seed = 1) # This step size works well, but the Stan-chosen step size works just as well
Fit_NM_T2_Dam <- sampling(model, data = StanData, chains = 4, iter = 1000, cores = 2, control = list(stepsize = 0.001), seed = 1)

system("say The run is done!")

## Summarizing the fit
Fit_NM_T2_Dam_Summ <- print(Fit_NM_T2_Dam, pars=c("r", "b", "c", "u", "sigma", "z_init"),
                            probs=c(0.1, 0.5, 0.9), digits = 3)

## Model checks ##
parms <- c("r", "b", "c", "u")
Output <- rstan::extract(Fit_NM_T2_Dam, permuted=TRUE, include=TRUE)
fit_Trace <- stan_trace(Fit_NM_T2_Dam, parms); fit_Trace
fit_Pairs <- mcmc_pairs(Fit_NM_T2_Dam, parms); fit_Pairs
fit_Dens <- mcmc_dens(Fit_NM_T2_Dam, parms); fit_Dens


Fit_NM_T2_Dam_Array <- as.array(Fit_NM_T2_Dam)[,,-4] 
dimnames(Fit_NM_T2_Dam_Array)[[2]] <- c("Chain 1", "Chain 2", "Chain 3", "Chain 4")
# dimnames(Fit_NM_T2_Dam_Array)[[2]] <- c("Chain 1", "Chain 2")
dimnames(Fit_NM_T2_Dam_Array)[[3]][1:6] <- c("r = 2.5", "O = 0.008", "h = 0.06", "u = 0.2", "c = 0.2", "b = 35")
parms <- c("r = 2.5", "O = 0.008", "h = 0.06", "u = 0.2", "c = 0.2", "b = 35")
color_scheme_set("viridisA")
mcmc_trace(Fit_NM_T2_Dam_Array, parms)
mcmc_pairs(Fit_NM_T2_Dam_Array, parms)
mcmc_dens(Fit_NM_T2_Dam_Array, parms)


