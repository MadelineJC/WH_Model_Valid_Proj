library(rstan)
library(deSolve)
library(GillespieSSA)
library(gdata)
library(bayesplot)
library(tidyverse)
library(brms)
library(beepr)
library(here)

#### Estimate all parms ####
## Import model
Model <- stan_model("03_Fit/NowakMay_T3_Fit_Oh.stan")
# https://rpubs.com/kaz_yos/stan_jacobian for warning message
# Also, here: https://mc-stan.org/docs/2_18/stan-users-guide/changes-of-variables.html 

## Clean data
InfecData <- read.csv('02_Data/NowakMay_StochSim_T3.csv')
N <- length(InfecData$Time) - 1
ts <- 1:N
y_init <- c(InfecData$Parasites[1], InfecData$ImmuneCells[1])
y <- as.matrix(InfecData[2:(N + 1), 2:3])
y <- cbind(y[ , 1], y[ , 2]); 
Data <- list(N = N, ts = ts, y_init = y_init, y = y)

# Fitting
Fit_NM_T3 <- sampling(Model, data = Data, chains = 4, iter = 1000, cores = 2, seed = 1); beep(3)

## Cleaning the fit
Fit_NM_T3_Array <- as.array(Fit_NM_T3)[,,-4] 
dimnames(Fit_NM_T3_Array)[[2]] <- c("Chain 1", "Chain 2", "Chain 3", "Chain 4")
# dimnames(Fit_NM_T3_Array)[[3]][1:6] <- c("r = 3", "O = 0.0015", "h = 0.2", "b = 35", "c = 0.2", "u = 0.5")
# parms <- c("r = 3", "O = 0.0015", "h = 0.2", "b = 35", "c = 0.2", "u = 0.5")
dimnames(Fit_NM_T3_Array)[[3]][1:4] <- c("r = 3", "b = 35", "c = 0.2", "u = 0.5")
parms <- c("r = 3", "b = 35", "c = 0.2", "u = 0.5")

## Summarising the fit
Fit_NM_T3_Summ <- print(Fit_NM_T3, pars=c("theta", "sigma", "z_init"),
                            probs=c(0.1, 0.5, 0.9), digits = 3)
