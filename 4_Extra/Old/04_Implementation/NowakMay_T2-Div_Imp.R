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
Model <- stan_model(here("./03_Fit/NowakMay_T2-Div_Fit_TRIAL.stan"))
# https://rpubs.com/kaz_yos/stan_jacobian for warning message
# Also, here: https://mc-stan.org/docs/2_18/stan-users-guide/changes-of-variables.html 

## Clean data
InfecData <- read_csv(here('./02_Data/NowakMay_StochSim_T2-Div.csv'))
N <- length(InfecData$Time) - 1
ts <- 1:N
y_init <- c(InfecData$Parasites[1], InfecData$ImmuneCells[1])
y <- as.matrix(InfecData[2:(N + 1), 2:3])
y <- cbind(y[ , 1], y[ , 2]); 
Data <- list(N = N, ts = ts, y_init = y_init, y = y)

# Fitting
Fit_NM_T2_Div <- sampling(Model, data = Data, chains = 4, iter = 1000, cores = 2, seed = 1)
system("say The run is done!")

## Summarising the fit
Fit_NM_T2_Div_Summ <- print(Fit_NM_T2_Div, pars=c("theta", "sigma", "z_init"),
                            probs=c(0.1, 0.5, 0.9), digits = 3)

## Cleaning the fit
Fit_NM_T2_Div_Array <- as.array(Fit_NM_T2_Div)[,,-4] 
dimnames(Fit_NM_T2_Div_Array)[[2]] <- c("Chain 1", "Chain 2", "Chain 3", "Chain 4")
# dimnames(Fit_NM_T2_Osc_Array)[[2]] <- c("Chain 1", "Chain 2")
dimnames(Fit_NM_T2_Div_Array)[[3]][1:5] <- c("r = 1", "O = 0.5", "h = 0.1", "c = 0.8", "u = 0.2")
parms <- c("r = 1", "O = 0.5", "h = 0.1", "c = 0.8", "u = 0.2")
color_scheme_set("viridisA")
mcmc_trace(Fit_NM_T2_Div_Array,parms)
