## Packages ##
library(rstan)
library(bayesplot)
library(tidyverse)
library(deSolve)
library(GillespieSSA)

#### Ex. Marty sent ===========================================================
# Source: https://medewitt.github.io/resources/stan_ar_models.html 

## Data generation ##
fake_data <- arima.sim(n = 200, model = list(ar = c(0.2, 0.5, 0.05)))
ts.plot(fake_data, main = "Time Series Plot of Our Fake Data")

## The model ##
write("
data {
  int<lower=0> K;  // Order of Autoregression
  int<lower=0> N; // number of observations
  real y[N];      // Outcome
}
parameters {
  real alpha;
  real beta[K];
  real sigma;
}
model {
  for (n in (K+1):N) {
    real mu = alpha;
    for (k in 1:K)
      mu += beta[k] * y[n-k];
      y[n] ~ normal(mu, sigma);
  }
}", 
"3_Extras/AutoregressiveModel/GitExample.stan")

## Compile model ##
model <- stan_model("3_Extras/AutoregressiveModel/GitExample.stan")

## Format data ##
stan_dat <- list(
  N = length(fake_data),
  K = 3,
  y = as.vector(fake_data)
)

## Fitting ##
fit <- sampling(model, stan_dat,
                cores = 4, iter = 2000,
                refresh = 0, chains = 4)

## Model checks ##
fit_summ <- print(fit, pars = c("alpha", "beta", "sigma"),
                  probs = c(0.1, 0.5, 0.9), digits = 3)
parms <- c("alpha", "beta[1]", "beta[2]", "beta[3]", "sigma")
Output <- rstan::extract(fit, permuted = TRUE, include = TRUE)
fit_Trace <- stan_trace(fit,parms); fit_Trace
fit_Pairs <- mcmc_pairs(fit,parms); fit_Pairs
fit_Dens <- mcmc_dens(fit,parms); fit_Dens

# Inference
print(fit, pars = "beta")

## Okay, so while the model seems to think it knows what it's doing, these estimates actually aren't incredible.
## Recall that the true values are 0.2, 0.5, and 0.05.

#### Stan user guide ==========================================================
# write("
# data {
#   int<lower=0> N;
#   vector[N] y;
# }
# parameters {
#   real alpha;
#   real beta;
#   real<lower=0> sigma;
# }
# model {
#   for (n in 2:N)
#     y[n] ~ normal(alpha + beta * y[n-1], sigma);
# }",
# "3_Extras/AutoregressiveModel/StanUserGuide_Model.stan")

#### Single Species: Logistic Growth ==========================================

## Let's generate some fake data ##
## Deterministic ##
LogisticGrowth <- function(t,y,p){
  r <- p[1]; K <- p[2]
  N <- y[1]
  dN = r*N*(1 - (N/K))
  list(c(dN))
}
r <- 0.1; K <- 1000
parms <- c(r, K)
N0 <- 10
SysInit <- c(N0)
TT <- seq(1,200,1) 
results <- lsoda(SysInit, TT, LogisticGrowth, parms)
N <- results[,2]
LogisticGrowth_Det = data.frame(results); colnames(LogisticGrowth_Det) <- c("Time", "Abundance")
LogisticGrowth_Det_Plot <- ggplot(LogisticGrowth_Det, aes(x = Time))+
  geom_line(aes(y = Abundance))+
  theme_minimal()
LogisticGrowth_Det_Plot

## Stochastic ##
x0 <- c(N = 10) 
a <- c("r*N",
       "(r*N^2)/K")
nu <- matrix(c(+1,-1), nrow = 1, byrow = TRUE)
r = 0.1; K = 1000
parms1 <- c(r = r, K = K)
tf = 200
method <- "OTL"
simName <- "LogisticGrowth"
set.seed(123)
LogisticGrowth_StochSim <- suppressWarnings(ssa(x0, a, nu, parms1, tf, method, simName,
                                  verbose = FALSE, 
                                  consoleInterval = 1, 
                                  censusInterval = 1, 
                                  maxWallTime = 30, 
                                  ignoreNegativeState = TRUE)) 
LogisticGrowth_Stoch <- LogisticGrowth_StochSim$data 
LogisticGrowth_Stoch <- as.data.frame(LogisticGrowth_Stoch); colnames(LogisticGrowth_Stoch) <- c("Time", "Abundance")
LogisticGrowth_Stoch_Plot <- ggplot(LogisticGrowth_Stoch, aes(x = Time))+
  geom_line(aes(y = Abundance))+ 
  theme_minimal()
LogisticGrowth_Stoch_Plot

## Wrangle data into form Stan likes ##
N <- length(LogisticGrowth_Stoch$Time) - 1
ts <- 1:N
y_init <- c(LogisticGrowth_Stoch$Abundance[1])
# y <- as.matrix(LogisticGrowth_Stoch[2:(N + 1), 2:2])
# y <- as.vector(LogisticGrowth_Stoch$Abundance)
y <- as.vector(LogisticGrowth_Stoch[2:(N + 1), 2:2])
StanData <- list(N = N, ts = ts, y_init = y_init, y = y)

## Let's compile this model, babey ##
write("
data {
  int<lower=0> N; // Number of observations
  //real ts[N];
  //real y_init[1];             
  //real<lower = 0> y[N, 1]; // Outcome
  vector[N] y;
}
parameters {
  real<lower = 0> r; // Growth rate
  real<lower = 0> K; // Carrying capacity
  real<lower=0> sigma; // Error
}
model {
  r ~ normal(0, 10);
  K ~ normal(1000,10);
  //sigma ~ lognormal(-1, 1);
  for (t in 2:N) {
    y[t] ~ normal(y[t-1] + r*y[t-1]*(1-y[t-1]/K), sigma);
    //real mu = y[t-1] + r*y[t-1]*(1-y[t-1]/K);
      //y[t] ~ normal(mu, sigma);
  }
}",
"3_Extras/AutoregressiveModel/Autoregressive_LogisticGrowth.stan")

model <- stan_model("3_Extras/AutoregressiveModel/Autoregressive_LogisticGrowth.stan")

## Fitting ##
fit <- sampling(model, data = StanData, chains = 4, iter = 2000, cores = 4, seed = 123)

## Model checks ##
fit_summ <- print(fit, pars = c("r", "K"),
                  probs = c(0.1, 0.5, 0.9), digits = 3)
parms <- c("r", "K")
Output <- rstan::extract(fit, permuted = TRUE, include = TRUE)
fit_Trace <- stan_trace(fit,parms); fit_Trace
fit_Pairs <- mcmc_pairs(fit,parms); fit_Pairs
fit_Dens <- mcmc_dens(fit,parms); fit_Dens

posterior <- rstan::extract(fit)
r_med <- median(posterior$r); K_med <- median(posterior$K)
r_low <- quantile(posterior$r, 0.025); K_low <- quantile(posterior$K, 0.025)
r_high <- quantile(posterior$r, 0.975); K_high <- quantile(posterior$K, 0.975)

## Med estimates ##
LogisticGrowth <- function(t,y,p){
  r <- p[1]; K <- p[2]
  N <- y[1]
  dN = r*N*(1 - (N/K))
  list(c(dN))
}
r <- r_med; K <- K_med
parms <- c(r, K)
N0 <- 10
SysInit <- c(N0)
TT <- seq(1,200,1) 
results <- lsoda(SysInit, TT, LogisticGrowth, parms)
N <- results[,2]
LG_Med = data.frame(results); colnames(LG_Med) <- c("Time", "Abundance")

## Low estimates ##
LogisticGrowth <- function(t,y,p){
  r <- p[1]; K <- p[2]
  N <- y[1]
  dN = r*N*(1 - (N/K))
  list(c(dN))
}
r <- r_low; K <- K_low
parms <- c(r, K)
N0 <- 10
SysInit <- c(N0)
TT <- seq(1,200,1) 
results <- lsoda(SysInit, TT, LogisticGrowth, parms)
N <- results[,2]
LG_Low = data.frame(results); colnames(LG_Low) <- c("Time", "Abundance")

## High estimates ##
LogisticGrowth <- function(t,y,p){
  r <- p[1]; K <- p[2]
  N <- y[1]
  dN = r*N*(1 - (N/K))
  list(c(dN))
}
r <- r_high; K <- K_high
parms <- c(r, K)
N0 <- 10
SysInit <- c(N0)
TT <- seq(1,200,1) 
results <- lsoda(SysInit, TT, LogisticGrowth, parms)
N <- results[,2]
LG_High = data.frame(results); colnames(LG_High) <- c("Time", "Abundance")

DetDF <- data.frame(seq(1,200,1), LogisticGrowth_Det$Abundance); colnames(DetDF) <- c("Time", "Det")
DataDF <- data.frame(LogisticGrowth_Stoch$Time, LogisticGrowth_Stoch$Abundance); colnames(DataDF) <- c("Time", "Abundance")
EstDF <- data.frame(LG_Med$Time, LG_Med$Abundance, LG_Low$Abundance, LG_High$Abundance); colnames(EstDF) <- c("Time", "Med_Est", "Low_Est", "High_Est")

## NOTE: I know this isn't a great way to plot the output/I'm not actually plotting the estimates, but here's a rough representation of what the model thinks is happening
Comp_Plot <- ggplot() +
  geom_line(data = DetDF, aes(x = Time, y = Det), linetype = "dotted") +
  geom_line(data = DataDF, aes(x = Time, y = Abundance)) +
  geom_line(data = EstDF, aes(x = Time, y = Med_Est), color = "springgreen4") +
  geom_line(data = EstDF, aes(x = Time, y = Low_Est), color = "royalblue1") +
  geom_line(data = EstDF, aes(x = Time, y = High_Est), color = "firebrick3") +
  labs(y = "Abundance")+
  theme_minimal()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
Comp_Plot

## Okay, so that actually went really well! I'd call this a working model. 

#### Two-Species: Logistic Growth =============================================
## Let's generate some fake data ##
## Stochastic ##
x0 <- c(N = 10) 
a <- c("r*N",
       "(r*N^2)/K")
nu <- matrix(c(+1,-1), nrow = 1, byrow = TRUE)
r = 0.1; K = 1000
parms1 <- c(r = r, K = K)
tf = 200
method <- "OTL"
simName <- "LogisticGrowth"
set.seed(123)
LogisticGrowth_StochSim <- suppressWarnings(ssa(x0, a, nu, parms1, tf, method, simName,
                                                verbose = FALSE, 
                                                consoleInterval = 1, 
                                                censusInterval = 1, 
                                                maxWallTime = 30, 
                                                ignoreNegativeState = TRUE)) 
LogisticGrowth_Stoch <- LogisticGrowth_StochSim$data 
LogisticGrowth_Stoch <- as.data.frame(LogisticGrowth_Stoch); colnames(LogisticGrowth_Stoch) <- c("Time", "Abundance")
LogisticGrowth_Stoch_Plot <- ggplot(LogisticGrowth_Stoch, aes(x = Time))+
  geom_line(aes(y = Abundance))+ 
  theme_minimal()
LogisticGrowth_Stoch_Plot

x0 <- c(N = 10) 
a <- c("r*N",
       "(r*N^2)/K")
nu <- matrix(c(+1,-1), nrow = 1, byrow = TRUE)
r = 0.2; K = 1200
parms1 <- c(r = r, K = K)
tf = 200
method <- "OTL"
simName <- "LogisticGrowth"
set.seed(123)
LogisticGrowth_StochSim <- suppressWarnings(ssa(x0, a, nu, parms1, tf, method, simName,
                                                verbose = FALSE, 
                                                consoleInterval = 1, 
                                                censusInterval = 1, 
                                                maxWallTime = 30, 
                                                ignoreNegativeState = TRUE)) 
LogisticGrowth_Stoch_2 <- LogisticGrowth_StochSim$data 
LogisticGrowth_Stoch_2 <- as.data.frame(LogisticGrowth_Stoch_2); colnames(LogisticGrowth_Stoch_2) <- c("Time", "Abundance")
LogisticGrowth_Stoch_Plot <- ggplot(LogisticGrowth_Stoch, aes(x = Time))+
  geom_line(aes(y = Abundance))+ 
  theme_minimal()
LogisticGrowth_Stoch_Plot

## Wrangle data into form Stan likes ##
N <- length(LogisticGrowth_Stoch$Time) - 1
ts <- 1:N
Y_init <- c(LogisticGrowth_Stoch$Abundance[1])
Y <- as.vector(LogisticGrowth_Stoch[2:(N + 1), 2:2])

y_init <- c(LogisticGrowth_Stoch_2$Abundance[1])
y <- as.vector(LogisticGrowth_Stoch_2[2:(N + 1), 2:2])
StanData <- list(N = N, ts = ts, Y_init = Y_init, Y = Y, y_init = y_init, y = y)

## Let's compile this model, babey ##
write("
data {
  int<lower = 0> N; // Number of observations
  vector [N] Y;
  vector [N] y;
}
parameters {
  real<lower = 0> R; // Growth rate of pop'n 1 (Y)
  real<lower = 0> K; // Carrying capacity of pop'n 1 (Y)
  real<lower = 0> r; // Growth rate of pop'n 2 (y)
  real<lower = 0> k; // Carrying capacity of pop'n 2 (y)
  real<lower = 0> sigma; // Error
}
model {
  R ~ normal(0.1, 1);
  K ~ normal(1000, 1);
  r ~ normal(0.2, 1);
  k ~ normal(1200, 1);
  //sigma ~ lognormal(-1, 1);
  for (t in 2:N) {
    Y[t] ~ normal(Y[t-1] + R*Y[t-1]*(1-Y[t-1]/K), sigma); // Pop'n 1
    y[t] ~ normal(y[t-1] + r*y[t-1]*(1-y[t-1]/k), sigma); // Pop'n 2
  }
}",
"3_Extras/AutoregressiveModel/Autoregressive_TwoSp-LogisticGrowth.stan")

model <- stan_model("3_Extras/AutoregressiveModel/Autoregressive_TwoSp-LogisticGrowth.stan")

## Fitting ##
fit <- sampling(model, data = StanData, chains = 4, iter = 2000, cores = 4, seed = 123)

## Model checks ##
fit_summ <- print(fit, pars = c("R", "K", "r", "k"),
                  probs = c(0.1, 0.5, 0.9), digits = 3)
parms <- c("R", "K", "r", "k")
Output <- rstan::extract(fit, permuted = TRUE, include = TRUE)
fit_Trace <- stan_trace(fit,parms); fit_Trace
fit_Pairs <- mcmc_pairs(fit,parms); fit_Pairs
fit_Dens <- mcmc_dens(fit,parms); fit_Dens

## Let's do some plotting to compare data to estimates ##
posterior <- rstan::extract(fit)
R_med <- median(posterior$R); K_med <- median(posterior$K)
r_med <- median(posterior$r); k_med <- median(posterior$k)
R_low <- quantile(posterior$R, 0.025); K_low <- quantile(posterior$K, 0.025)
r_low <- quantile(posterior$r, 0.025); k_low <- quantile(posterior$k, 0.025)
R_high <- quantile(posterior$R, 0.975); K_high <- quantile(posterior$K, 0.975)
r_high <- quantile(posterior$r, 0.975); k_high <- quantile(posterior$k, 0.975)

## Med estimates ##
LogisticGrowth <- function(t,y,p){
  r <- p[1]; K <- p[2]
  N <- y[1]
  dN = r*N*(1 - (N/K))
  list(c(dN))
}
r <- R_med; K <- K_med
parms <- c(r, K)
N0 <- 10
SysInit <- c(N0)
TT <- seq(1,200,1) 
results <- lsoda(SysInit, TT, LogisticGrowth, parms)
N <- results[,2]
LG_Med = data.frame(results); colnames(LG_Med) <- c("Time", "Abundance")

LogisticGrowth <- function(t,y,p){
  r <- p[1]; K <- p[2]
  N <- y[1]
  dN = r*N*(1 - (N/K))
  list(c(dN))
}
r <- r_med; K <- k_med
parms <- c(r, K)
N0 <- 10
SysInit <- c(N0)
TT <- seq(1,200,1) 
results <- lsoda(SysInit, TT, LogisticGrowth, parms)
N <- results[,2]
LG_Med_2 = data.frame(results); colnames(LG_Med_2) <- c("Time", "Abundance")

## Low estimates ##
LogisticGrowth <- function(t,y,p){
  r <- p[1]; K <- p[2]
  N <- y[1]
  dN = r*N*(1 - (N/K))
  list(c(dN))
}
r <- R_low; K <- K_low
parms <- c(r, K)
N0 <- 10
SysInit <- c(N0)
TT <- seq(1,200,1) 
results <- lsoda(SysInit, TT, LogisticGrowth, parms)
N <- results[,2]
LG_Low = data.frame(results); colnames(LG_Low) <- c("Time", "Abundance")

LogisticGrowth <- function(t,y,p){
  r <- p[1]; K <- p[2]
  N <- y[1]
  dN = r*N*(1 - (N/K))
  list(c(dN))
}
r <- r_low; K <- k_low
parms <- c(r, K)
N0 <- 10
SysInit <- c(N0)
TT <- seq(1,200,1) 
results <- lsoda(SysInit, TT, LogisticGrowth, parms)
N <- results[,2]
LG_Low_2 = data.frame(results); colnames(LG_Low_2) <- c("Time", "Abundance")

## High estimates ##
LogisticGrowth <- function(t,y,p){
  r <- p[1]; K <- p[2]
  N <- y[1]
  dN = r*N*(1 - (N/K))
  list(c(dN))
}
r <- R_high; K <- K_high
parms <- c(r, K)
N0 <- 10
SysInit <- c(N0)
TT <- seq(1,200,1) 
results <- lsoda(SysInit, TT, LogisticGrowth, parms)
N <- results[,2]
LG_High = data.frame(results); colnames(LG_High) <- c("Time", "Abundance")

LogisticGrowth <- function(t,y,p){
  r <- p[1]; K <- p[2]
  N <- y[1]
  dN = r*N*(1 - (N/K))
  list(c(dN))
}
r <- r_high; K <- k_high
parms <- c(r, K)
N0 <- 10
SysInit <- c(N0)
TT <- seq(1,200,1) 
results <- lsoda(SysInit, TT, LogisticGrowth, parms)
N <- results[,2]
LG_High_2 = data.frame(results); colnames(LG_High_2) <- c("Time", "Abundance")

LogisticGrowth_Stoch_2 <- LogisticGrowth_Stoch_2[1:187, ]
DataDF <- data.frame(LogisticGrowth_Stoch$Time, LogisticGrowth_Stoch$Abundance, LogisticGrowth_Stoch_2$Abundance); colnames(DataDF) <- c("Time", "Abundance", "Abundance_2")
EstDF <- data.frame(LG_Med$Time, LG_Med$Abundance, LG_Low$Abundance, LG_High$Abundance, LG_Med_2$Abundance, LG_Low_2$Abundance, LG_High_2$Abundance); colnames(EstDF) <- c("Time", "Med_Est", "Low_Est", "High_Est", "Med_2_Est", "Low_2_Est", "High_2_Est")

## NOTE: I know this isn't a great way to plot the output/I'm not actually plotting the estimates, but here's a rough representation of what the model thinks is happening
Comp_Plot <- ggplot() +
  # geom_line(data = DetDF, aes(x = Time, y = Det), linetype = "dotted") +
  geom_line(data = DataDF, aes(x = Time, y = Abundance)) +
  geom_line(data = DataDF, aes(x = Time, y = Abundance_2), linetype = "dotted") +
  geom_line(data = EstDF, aes(x = Time, y = Med_Est), color = "springgreen4") +
  geom_line(data = EstDF, aes(x = Time, y = Med_2_Est), color = "springgreen4", linetype = "dotted") +
  geom_line(data = EstDF, aes(x = Time, y = Low_Est), color = "royalblue1") +
  geom_line(data = EstDF, aes(x = Time, y = Low_2_Est), color = "royalblue1", linetype = "dotted") +
  geom_line(data = EstDF, aes(x = Time, y = High_Est), color = "firebrick3") +
  geom_line(data = EstDF, aes(x = Time, y = High_2_Est), color = "firebrick3", linetype = "dotted") +
  labs(main = "Solid = Pop'n 1; Dotted = Pop'n 2", y = "Abundance")+
  theme_minimal()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
Comp_Plot

## Okay, so this also worked, but these populations don't interact. So, the bar is low. 

#### Two Interacting Species: Antia et al. ====================================
## Deterministic ##
r = 0.2; k = 0.01; p = 1; o = 1000
P <- c(); I <- c()
P[1] = 1; I[1] = 1
ts = 1; t_stop = 55
for (t in 2:t_stop){
  P[t] = (1 + r*ts - (1-exp(-k*I[t - 1]*ts)))*P[t - 1]
  I[t] = (1 + p*(P[t - 1]/(P[t - 1] + o))*ts)*I[t - 1]
}
Antia_Det <- data.frame(seq(1, t_stop, 1), P, I); colnames(Antia_Det) <- c("Time", "P", "I")

## Stochastic ##
DiscreteModel_Dem <- function(P_Last, I_Last){
  ts = 1
  r = 0.2; k = 0.01; p = 1; o = 1000
  
  P_B = rpois(1, r*P_Last*ts) 
  P_D = rbinom(1, P_Last, (1-exp(-k*I_Last*ts))) 
  I_B = rpois(1, p*I_Last*((P_Last/(P_Last + o)))*ts) 
  
  P_Next = P_Last + P_B - P_D
  I_Next = I_Last + I_B
  
  P_Last <- P_Next
  I_Last <- I_Next
  
  if (P_Last < 0){
    P_Last <- 0
  }
  
  if (P_Last > 10^9){
    P_Last <- 0
    I_Last <- 0
  }
  return(c(P_Last, I_Last))
}

ts <- 1; TT <- 55
P_Last <- 1; I_Last <- 1
Antia_Stoch <- data.frame()
set.seed(6)
for (j in 1:TT){
  Output = DiscreteModel_Dem(P_Last, I_Last)
  P_Last = Output[1]
  I_Last = Output[2]
  Addition <- c(P_Last, I_Last)
  Antia_Stoch <- data.frame(rbind(Antia_Stoch, Addition))
}
Antia_Stoch <- cbind(seq(1, TT, 1), Antia_Stoch)
colnames(Antia_Stoch) <- c("Time", "P", "I")

## Plotting the deterministic and stochastic implementations together:
plot(1:55, Antia_Det$P, type = "l", col = "forestgreen",
     xlab = "Time", ylab = "Abundance", ylim = c(0, 1000))
lines(1:55, Antia_Det$I, col = "cornflowerblue")
lines(1:55, Antia_Stoch$P, col = "forestgreen", lty = 3)
lines(1:55, Antia_Stoch$I, col = "cornflowerblue", lty = 3)
legend("topleft", legend = c("Deterministic", "Stochastic"), lty = c(1, 3))

## Data wrangling for Stan ##
x <- which(Antia_Stoch$P == 0)[1]
Antia_Stoch <- Antia_Stoch[1:x - 1, ]
N <- length(Antia_Stoch$Time) - 1
ts <- 1:N
P_init <- c(Antia_Stoch$P[1])
I_init <- c(Antia_Stoch$I[1])
P <- as.vector(Antia_Stoch[2:(N + 1), 2:2])
I <- as.vector(Antia_Stoch[2:(N + 1), 3:3])
StanData <- list(N = N, ts = ts, P_init = P_init, I_init = I_init, P = P, I = I)

## Write model ##
write("
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
}",
"3_Extras/AutoregressiveModel/Autoregressive_Antia.stan")

## Compile model ##
model <- stan_model("3_Extras/AutoregressiveModel/Autoregressive_Antia.stan")

## Fitting ##
fit <- sampling(model, data = StanData, chains = 4, iter = 10000, cores = 4, seed = 1)

## Model checks ##
fit_summ <- print(fit, pars = c("r", "k", "p", "o"),
                  probs = c(0.1, 0.5, 0.9), digits = 3)
parms <- c("r", "k", "p", "o")
Output <- rstan::extract(fit, permuted = TRUE, include = TRUE)
fit_Trace <- stan_trace(fit,parms); fit_Trace
fit_Pairs <- mcmc_pairs(fit,parms); fit_Pairs
fit_Overlay <- mcmc_dens_overlay(fit,parms); fit_Overlay

r = mean(Output[["r"]]); k = mean(Output[["k"]]); p = mean(Output[["p"]]); o = mean(Output[["o"]])
P <- c(); I <- c()
P[1] = 1; I[1] = 1
ts = 1; t_stop = 55
for (t in 2:t_stop){
  P[t] = (1 + r*ts - (1-exp(-k*I[t - 1]*ts)))*P[t - 1]
  I[t] = (1 + p*(P[t - 1]/(P[t - 1] + o))*ts)*I[t - 1]
}
Antia_Est <- data.frame(seq(1, t_stop, 1), P, I); colnames(Antia_Det) <- c("Time", "P", "I")

plot(1:nrow(Antia_Stoch), Antia_Stoch$P, type = "l", col = "forestgreen",
     xlab = "Time", ylab = "Abundance", xlim = c(0, 55), ylim = c(0, 1000))
lines(1:nrow(Antia_Stoch), Antia_Stoch$I, col = "cornflowerblue")
lines(1:nrow(Antia_Est), Antia_Est$P, col = "forestgreen", lty = 3)
lines(1:nrow(Antia_Est), Antia_Est$I, col = "cornflowerblue", lty = 3)
legend("topleft", legend = c("Data", "Estimate"), lty = c(1, 3))

#### Mod. Fenton and Perkins ==================================================
## Deterministic ##
r = 2.5; B = 0.008; h = 0.06; b = 35; e = 0.2; delta = 0.2
P <- c(); I <- c()
P[1] = 80; I[1] = 200
ts = 0.1; t_stop = 50/ts
for (t in 2:t_stop){
  P[t] = (1 + r*ts - (1-exp((-B/(1 + B*h*P[t - 1]))*I[t - 1]*ts)))*P[t - 1]
  # I[t] = I[t - 1] + b*ts + e*(B*P[t - 1]/(1 + B*P[t - 1]*h))*I[t - 1]*ts - (1-exp(-delta*ts))*I[t - 1]
  I[t] = (1 + b*ts*(1/I[t - 1]) + e*(B*P[t - 1]/(1 + B*P[t - 1]*h))*ts - (1-exp(-delta*ts)))*I[t - 1]
}
ModFentonPerkins_Det <- data.frame(seq(0.1, 50, 0.1), P, I); colnames(ModFentonPerkins_Det) <- c("Time", "P", "I")
plot(ModFentonPerkins_Det$Time, ModFentonPerkins_Det$P, type = "l", col = "forestgreen",
     main = "Deterministic", xlab = "Time", ylab = "Abundance", ylim = c(0, max(ModFentonPerkins_Det[ , 2:3])))
lines(ModFentonPerkins_Det$Time, ModFentonPerkins_Det$I, col = "cornflowerblue")

## Stochastic ##
DiscreteModel_Dem <- function(P_Last, I_Last){
  ts = 0.1
  r = 2.5; B = 0.008; h = 0.06; b = 35; e = 0.2; delta = 0.2
  
  P_B = rpois(1, r*P_Last*ts) 
  P_D = rbinom(1, P_Last, (1-exp((-B/(1 + B*h*P_Last))*I_Last*ts))) ## ASK MARTY ABOUT THIS LINE
  I_B = rpois(1, (b + e*(B*P_Last/(1 + B*P_Last*h))*I_Last)*ts)
  I_D = rbinom(1, I_Last, (1-exp(-delta*ts))) ## ASK MARTY ABOUT THIS LINE
  
  P_Next = P_Last + P_B - P_D
  I_Next = I_Last + I_B - I_D
  
  P_Last <- P_Next
  I_Last <- I_Next
  
  ## To prevent the populations from going negative
  if (P_Last < 0){
    P_Last <- 0
  }
  
  return(c(P_Last, I_Last))
}

ts <- 0.1; TT <- 20/ts
P_Last <- 80; I_Last <- 200
ModFentonPerkins_Stoch <- data.frame()
for (j in 1:TT){
  Output = DiscreteModel_Dem(P_Last, I_Last)
  P_Last = Output[1]
  I_Last = Output[2]
  Addition <- c(P_Last, I_Last)
  ModFentonPerkins_Stoch <- data.frame(rbind(ModFentonPerkins_Stoch, Addition))
}
ModFentonPerkins_Stoch <- cbind(seq(1, TT, 1), ModFentonPerkins_Stoch)
colnames(ModFentonPerkins_Stoch) <- c("Time","P", "I")

## Plotting the deterministic and stochastic implementations together:
plot(1:nrow(ModFentonPerkins_Det), ModFentonPerkins_Det$P, type = "l", col = "forestgreen",
     xlab = "Time", ylab = "Abundance", ylim = c(0, 1000))
lines(1:nrow(ModFentonPerkins_Det), ModFentonPerkins_Det$I, col = "cornflowerblue")
lines(1:nrow(ModFentonPerkins_Stoch), ModFentonPerkins_Stoch$P, col = "forestgreen", lty = 3)
lines(1:nrow(ModFentonPerkins_Stoch), ModFentonPerkins_Stoch$I, col = "cornflowerblue", lty = 3)
legend("topleft", legend = c("Deterministic", "Stochastic"), lty = c(1, 3))

## Data wrangling for Stan ##
N <- length(ModFentonPerkins_Stoch$Time) - 1
ts <- 1:N
P <- as.vector(ModFentonPerkins_Stoch[2:(N + 1), 2:2])
I <- as.vector(ModFentonPerkins_Stoch[2:(N + 1), 3:3])
StanData <- list(N = N, ts = ts, P = P, I = I)

write("
data {
  int<lower=0> N;
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
  r ~ normal(2.5, 1); // 2.5
  B ~ normal(0, 1); // 0.008
  h ~ normal(0, 1); // 0.06
  b ~ normal(35, 1); // 35
  e ~ normal(0, 1); // 0.2
  delta ~ normal(0, 1); // 0.2
  sigma ~ lognormal(-1, 1);
  for (t in 2:N) {
    P[t] ~ normal(P[t-1] + r*P[t-1] - (1-exp((-B/(1 + B*h*P[t-1]))*I[t-1]))*P[t-1], sigma);
    I[t] ~ normal(I[t-1] + b + e*(B*P[t-1]/(1 + B*P[t-1]*h))*I[t-1] - (1-exp(-delta))*I[t-1], sigma);
  }
}",
"Autoregressive_Attempt8.stan")

## Compile model ##
model <- stan_model("Autoregressive_ModFentonPerkins.stan")

## Fitting ##
fit <- sampling(model, data = StanData, chains = 4, iter = 2000, cores = 4, seed = 1)

## Model checks ##
fit_summ <- print(fit, pars = c("r", "B", "h", "b", "e", "delta"),
                  probs = c(0.1, 0.5, 0.9), digits = 3)
parms <- c("r", "B", "h", "b", "e", "delta")
Output <- rstan::extract(fit, permuted = TRUE, include = TRUE)
fit_Trace <- stan_trace(fit,parms); fit_Trace
fit_Pairs <- mcmc_pairs(fit,parms); fit_Pairs
fit_Dens <- mcmc_dens(fit,parms); fit_Dens

## Putting the parameter estimates into the deterministic model:
r = 0.241; B = 0.001; h = 0.07; b = 4.344; e = 0.194; delta = 0.024
P <- c(); I <- c()
P[1] = 80; I[1] = 200
ts = 0.1; t_stop = 100/ts
for (t in 2:t_stop){
  P[t] = (1 + r*ts - (1-exp((-B/(1 + B*h*P[t - 1]))*I[t - 1]*ts)))*P[t - 1]
  # I[t] = I[t - 1] + b*ts + e*(B*P[t - 1]/(1 + B*P[t - 1]*h))*I[t - 1]*ts - (1-exp(-delta*ts))*I[t - 1]
  I[t] = (1 + b*ts*(1/I[t - 1]) + e*(B*P[t - 1]/(1 + B*P[t - 1]*h))*ts - (1-exp(-delta*ts)))*I[t - 1]
}
ModFentonPerkins_Det <- data.frame(seq(1,t_stop,1), P, I); colnames(ModFentonPerkins_Det) <- c("Time", "P", "I")
ModFentonPerkins_Det_Plot <- ggplot(ModFentonPerkins_Det,aes(x = Time))+
  geom_line(aes(y = P))+ 
  geom_line(aes(y = I))+
  theme_minimal()
ModFentonPerkins_Det_Plot

#### Mod. Fenton and Perkins (2) ==============================================
## Deterministic ##
r = 2.5; B = 0.012; h = 0.075; b = 35; e = 0.3; delta = 0.41
P <- c(); I <- c()
P[1] = 80; I[1] = 200
ts = 0.1; t_stop = 60/ts
for (t in 2:t_stop){
  # P[t] = P[t - 1] + r*P[t - 1]*ts - (((B*P[t - 1])/(1 + B*P[t - 1]*h))*I[t - 1]*ts)
  # I[t] = I[t - 1] + b*ts + (e*I[t - 1]*((B*P[t - 1])/(1 + B*P[t - 1]*h))*ts) - delta*I[t - 1]*ts
  P[t] = (1 + r*ts - (1-exp((-B/(1 + B*h*P[t - 1]))*I[t - 1]*ts)))*P[t - 1]
  I[t] = (1 + b*ts*(1/I[t - 1]) + e*(B*P[t - 1]/(1 + B*P[t - 1]*h))*ts - (1-exp(-delta*ts)))*I[t - 1]
}
ModFentonPerkins_Det <- data.frame(seq(0.1, 60, 0.1), P, I); colnames(ModFentonPerkins_Det) <- c("Time", "P", "I")
plot(ModFentonPerkins_Det$Time, ModFentonPerkins_Det$P, type = "l", col = "forestgreen",
     main = "Deterministic", xlab = "Time", ylab = "Abundance", ylim = c(0, max(ModFentonPerkins_Det[ , 2:3])))
lines(ModFentonPerkins_Det$Time, ModFentonPerkins_Det$I, col = "cornflowerblue")

## Stochastic ##
DiscreteModel_Dem <- function(P_Last, I_Last){
  ts = 0.1
  r = 2.5; B = 0.012; h = 0.075; b = 35; e = 0.3; delta = 0.41
  
  P_B = rpois(1, r*P_Last*ts) 
  P_D = rbinom(1, P_Last, (1-exp((-B/(1 + B*h*P_Last))*I_Last*ts))) ## ASK MARTY ABOUT THIS LINE
  I_B = rpois(1, (b + e*(B*P_Last/(1 + B*P_Last*h))*I_Last)*ts)
  I_D = rbinom(1, I_Last, (1-exp(-delta*ts))) ## ASK MARTY ABOUT THIS LINE
  
  P_Next = P_Last + P_B - P_D
  I_Next = I_Last + I_B - I_D
  
  P_Last <- P_Next
  I_Last <- I_Next
  
  ## To prevent the populations from going negative
  if (P_Last < 0){
    P_Last <- 0
  }
  
  return(c(P_Last, I_Last))
}

ts <- 0.1; t_stop <- 60/ts
P_Last <- 80; I_Last <- 200
ModFentonPerkins_Stoch <- data.frame()
set.seed(3)
for (j in 1:t_stop){
  Output = DiscreteModel_Dem(P_Last, I_Last)
  P_Last = Output[1]
  I_Last = Output[2]
  Addition <- c(P_Last, I_Last)
  ModFentonPerkins_Stoch <- data.frame(rbind(ModFentonPerkins_Stoch, Addition))
}
ModFentonPerkins_Stoch <- data.frame(seq(0.1, 60, 0.1), ModFentonPerkins_Stoch); colnames(ModFentonPerkins_Stoch) <- c("Time","P", "I")
plot(ModFentonPerkins_Stoch$Time, ModFentonPerkins_Stoch$P, type = "l", col = "forestgreen",
     main = "Stochastic", xlab = "Time", ylab = "Abundance", ylim = c(0, max(ModFentonPerkins_Stoch[ , 2:3])))
lines(ModFentonPerkins_Stoch$Time, ModFentonPerkins_Stoch$I, col = "cornflowerblue")

## Plotting the deterministic and stochastic implementations together:
plot(1:nrow(ModFentonPerkins_Det), ModFentonPerkins_Det$P, type = "l", col = "forestgreen",
     xlab = "Time", ylab = "Abundance", ylim = c(0, 500))
lines(1:nrow(ModFentonPerkins_Det), ModFentonPerkins_Det$I, col = "cornflowerblue")
lines(1:nrow(ModFentonPerkins_Stoch), ModFentonPerkins_Stoch$P, col = "forestgreen", lty = 3)
lines(1:nrow(ModFentonPerkins_Stoch), ModFentonPerkins_Stoch$I, col = "cornflowerblue", lty = 3)
legend("topleft", legend = c("Deterministic", "Stochastic"), lty = c(1, 3))

## Data wrangling for Stan ##
x <- which(ModFentonPerkins_Stoch$P == 0)[1]
# ModFentonPerkins_Stoch <- ModFentonPerkins_Stoch[1:x - 1, ]
N <- length(ModFentonPerkins_Stoch$Time) - 1
ts <- 1:N
P_init <- c(ModFentonPerkins_Stoch$P[80])
I_init <- c(ModFentonPerkins_Stoch$I[200])
P <- as.vector(ModFentonPerkins_Stoch[2:(N + 1), 2:2])
I <- as.vector(ModFentonPerkins_Stoch[2:(N + 1), 3:3])
StanData <- list(N = N, ts = ts, P_init = P_init, I_init = I_init, P = P, I = I)

## Write model ##
write("
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
  r ~ normal(2.5, 1); // 2.5
  B ~ normal(0, 1); // 0.012
  h ~ normal(0, 1); // 0.075
  b ~ normal(35, 1); // 35
  e ~ normal(0, 1); // 0.3
  delta ~ normal(0, 1); // 0.41
  sigma ~ lognormal(-1, 1);
  for (t in 2:N) {
    // P[t] ~ normal(P[t-1] + r*P[t-1]*0.1 - (1-exp((-B/(1 + B*h*P[t-1]))*I[t-1]*0.1))*P[t-1], sigma);
    // I[t] ~ normal(I[t-1] + b*0.1 + e*(B*P[t-1]/(1 + B*P[t-1]*h))*I[t-1]*0.1 - (1-exp(-delta))*I[t-1], sigma);
    P[t] ~ normal((1 + r*0.1 - (1-exp((-B/(1 + B*h*P[t - 1]))*I[t - 1]*0.1)))*P[t - 1], sigma);
    I[t] ~ normal((1 + b*0.1*(1/I[t - 1]) + e*(B*P[t - 1]/(1 + B*P[t - 1]*h))*0.1 - (1-exp(-delta*0.1)))*I[t - 1], sigma);
  }
}",
"3_Extras/AutoregressiveModel/Autoregressive_ModFentonPerkins.stan")

## Compile model ##
model <- stan_model("3_Extras/AutoregressiveModel/Autoregressive_ModFentonPerkins.stan")

## Fitting ##
fit <- sampling(model, data = StanData, chains = 4, iter = 10000, cores = 4, seed = 1)

## Model checks ##
fit_summ <- print(fit, pars = c("r", "B", "h", "b", "e", "delta"),
                  probs = c(0.1, 0.5, 0.9), digits = 3)
parms <- c("r", "B", "h", "b", "e", "delta")
Output <- rstan::extract(fit, permuted = TRUE, include = TRUE)
fit_Trace <- stan_trace(fit,parms); fit_Trace
fit_Pairs <- mcmc_pairs(fit,parms); fit_Pairs
fit_Overlay <- mcmc_dens_overlay(fit,parms); fit_Overlay

r = mean(Output[["r"]]); B = mean(Output[["B"]]); h = mean(Output[["h"]]) 
b = mean(Output[["b"]]); e = mean(Output[["e"]]); delta = mean(Output[["delta"]])
P <- c(); I <- c()
P[1] = 80; I[1] = 200
ts = 0.1; t_stop = 60/ts
for (t in 2:t_stop){
  # P[t] = P[t - 1] + r*P[t - 1]*ts - (((B*P[t - 1])/(1 + B*P[t - 1]*h))*I[t - 1]*ts)
  # I[t] = I[t - 1] + b*ts + (e*I[t - 1]*((B*P[t - 1])/(1 + B*P[t - 1]*h))*ts) - delta*I[t - 1]*ts
  P[t] = (1 + r*ts - (1-exp((-B/(1 + B*h*P[t - 1]))*I[t - 1]*ts)))*P[t - 1]
  I[t] = (1 + b*ts*(1/I[t - 1]) + e*(B*P[t - 1]/(1 + B*P[t - 1]*h))*ts - (1-exp(-delta*ts)))*I[t - 1]
}
ModFentonPerkins_Est <- data.frame(seq(0.1, 60, 0.1), P, I); colnames(ModFentonPerkins_Est) <- c("Time", "P", "I")

plot(ModFentonPerkins_Stoch$Time, ModFentonPerkins_Stoch$P, type = "l", col = "forestgreen",
     xlab = "Time", ylab = "Abundance", ylim = c(0, 500))
lines(ModFentonPerkins_Stoch$Time, ModFentonPerkins_Stoch$I, col = "cornflowerblue")
lines(ModFentonPerkins_Est$Time, ModFentonPerkins_Est$P, col = "forestgreen", lty = 3)
lines(ModFentonPerkins_Est$Time, ModFentonPerkins_Est$I, col = "cornflowerblue", lty = 3)
legend("topleft", legend = c("Data", "Estimate"), lty = c(1, 3))



