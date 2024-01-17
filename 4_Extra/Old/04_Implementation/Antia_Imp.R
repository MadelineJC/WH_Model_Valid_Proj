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
StanData <- read_csv(here('./02_Data/Antia_StochSim.csv'))
N <- length(StanData$Time) - 1
ts <- 1:N
y_init <- c(StanData$Parasites[1], StanData$ImmuneCells[1])
y <- as.matrix(StanData[2:(N + 1), 2:3])
y <- cbind(y[ , 1], y[ , 2]); 
StanData <- list(N = N, ts = ts, y_init = y_init, y = y)

## Compile model ##
model <- stan_model(here('./03_Fit/Antia_Fit.stan'))

## Fitting ##
inits <- list(
  r = 0.2,
  k = 0.01,
  p = 1,
  o = 1000
); list_inits <- list(inits, inits, inits, inits)

Fit_Antia <- sampling(model, data = StanData, chains = 4, iter = 1000, cores = 4, seed = 123)

system("say The run is done!")

## Summarizing the fit ##
Fit_Antia_Summ <- print(Fit_Antia, pars=c("r", "k", "p", "o", "sigma", "z_init"),
                            probs=c(0.1, 0.5, 0.9), digits = 3)

## Model checks ##
parms <- c("r", "k", "p", "o")
Output <- rstan::extract(Fit_Antia, permuted=TRUE, include=TRUE)
fit_Trace <- stan_trace(Fit_Antia, parms); fit_Trace
fit_Pairs <- mcmc_pairs(Fit_Antia, parms); fit_Pairs
fit_Dens <- mcmc_dens(Fit_Antia, parms); fit_Dens 
