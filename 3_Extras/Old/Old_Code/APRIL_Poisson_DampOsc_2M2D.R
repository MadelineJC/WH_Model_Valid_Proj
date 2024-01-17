library(rstan)
library(gdata)
library(bayesplot)
library(tidyverse)
library(brms)
library(deSolve)
library(GillespieSSA)

#### Deterministic ####
TypeIIFR_Model <- function(t,y,p){
  r <- p[1]; O <- p[2]; h <- p[3]
  b <- p[4]; c <- p[5]; u <- p[6]
  P <- y[1]; H <- y[2]
  dP = P*r - H*(O*P/(1 + O*P*h))
  dH = b + H*(c*(O*P/(1 + O*P*h)) - u)
  list(c(dP,dH))
}
r <- 2.5; O <- 0.008; h <- 0.06
b <- 35; c <- 0.2; u <- 0.2
parms <- c(r,O,h,b,c,u)
P0 <- 10; H0 <- 350 
N0 <- c(P0,H0)
TT <- seq(0.1,400,0.1) 
results <- lsoda(N0,TT,TypeIIFR_Model,parms)
P <- results[,2]; H <- results[,3];
TypeII_Full_DF = data.frame(results)
TypeII_Full_Plot <- ggplot(TypeII_Full_DF,aes(x=time))+
  geom_line(aes(y=X1, color="cornflowerblue"))+ # Parasite
  geom_line(aes(y=X2,color="springgreen4"))+ # Host
  labs(title="Full Time-Series; Type II",
       x="Time",
       y="Population Abundance")+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"))+
  theme_minimal()
TypeII_Full_Plot 

#### Stochastic ####
x0 <- c(P=10,H=350) 
a <- c("P*r",
       "H*(O*P/1 + O*P*h)", 
       "b + H*c*(O*P/1 + O*P*h)",
       "H*u")
nu <- matrix(c(+1,-1,0,0,
               0,0,+1,-1), nrow = 2, byrow = TRUE)
r = 2.5; O = 0.008; h = 0.06; b = 35; c = 0.2; u = 0.2
parms1 <- c(r = r, O = O, h = h, b = b, c = c, u = u)
tf = 400
method <- "OTL"
simName <- "WH_IDD_TypeII_FULL"
set.seed(1905)
SimResults_Stoch_TypeII <- suppressWarnings(ssa(x0, a, nu, parms1, tf, method, simName,
                                                verbose = FALSE, 
                                                consoleInterval = 0,
                                                censusInterval = 0.1, 
                                                maxWallTime = 30, 
                                                ignoreNegativeState = TRUE)) 

Stoch_Data_TypeII <- SimResults_Stoch_TypeII$data 
Stoch_Data_TypeII <- as.data.frame(Stoch_Data_TypeII)
write.csv(Stoch_Data_TypeII,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/DampOsc_Data.csv", row.names = FALSE)

TypeII_Full_Stoch_Plot <- ggplot(Stoch_Data_TypeII,aes(x=t))+
  geom_line(aes(y=P, color="cornflowerblue"))+ # Parasite
  geom_line(aes(y=H,color="springgreen4"))+ # Host
  labs(title="W-H InfecDisDyn w/ Type II FcnalResp and DemStoch",
       x="Time",
       y="Population Abundance")+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"))
TypeII_Full_Stoch_Plot

#### Model Fitting ####

#### Fit 9 ####
write("// Stan Model;
functions{
  real[] dZ_dt(
    real t, // Time
    real[] Z, // System state {Parasite, Host}
    real[] theta, // Parms
    real[] x_r, // Real data; below is integer data
    int[] x_i){
    real P = Z[1]; // System state coded as an array, such that Z = (P,H)
    real H = Z[2];

    real r = theta[1];  // Parameters of the system, in order they appear
    real O = theta[2];
    real h = theta[3];
    real b = theta[4];
    real c = theta[5];
    real u = theta[6];

    real dP_dt = P*r - H*(O*P/(1 + O*P*h)); // Mechanistic model
    real dH_dt = b + H*(c*(O*P/(1 + O*P*h))-u);
    return({dP_dt,dH_dt}); // Return the system state
  }
}

data{
  int<lower=0>N; // Define N as non-negative integer
  real ts[N]; // Assigns time points to N
  real y0[2]; // Initial conditions for ODE
  real t0; // Initial time point
  int<lower=0>y[N,2]; // Define y as real and non-negative
}

parameters{
  real<lower=0>r;
  real<lower=0,upper=1>O;
  real<lower=0>h;
  real<lower=0>b;
  real<lower=0>c;
  real<lower=0,upper=1>u;
}

transformed parameters{
  // Stiff solver for ODE
  real Z[N,2]
  = integrate_ode_bdf(dZ_dt,y0,t0,ts,{r, O, h, b, c, u},rep_array(0.0,2),rep_array(0,2),1e-10,1e-10,2e4);
  for (i in 1:N) {
    Z[i,1] = Z[i,1] + 1e-6;
    Z[i,2] = Z[i,2] + 1e-6;
  }
}

model{ // 
  r~normal(2.5,1); // r
  O~uniform(0,1); // O; bounded between 0 and 1
  h~uniform(0,1); // h
  b~normal(35,1); // b
  c~normal(0.2,1); // c
  u~normal(0.2,1); // u; bounded between 0 and 1
  for (k in 1:2){
    y[ ,k]~poisson(Z[ ,k]);
  }
}

generated quantities {
  int y_rep[N, 2];
  for (k in 1:2) {
    for (n in 1:N)
      y_rep[n,k] = poisson_rng(Z[n,k]);
  }
}",
"Stan_Model_TypeII.stan")

stanc("Stan_Model_TypeII.stan") # To check that we wrote a file

Stan_TypeII <- stan_model("Stan_Model_TypeII.stan")

# Squeezing the data into a form that Stan gets
StochDataTypeII = read.csv('/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/DampOsc_Data.csv')
N <- length(StochDataTypeII$t)-1 # df is 2923; N is 2922
ts <- 1:N
y_init <- c(10,350) # Initial states, P = 10; H = 350
y0 <- array(y_init) # For the ODE solver
t0 <- 0.1
y <- as.matrix(StochDataTypeII[2:(N+1),2:3])
y <- cbind(y[ ,1],y[ ,2]); # This worked, sick; where y[,1] is P, and y[,2] is H
Stan_StochDataTypeII <- list(N=N,ts=ts,y0=y0,t0=t0,y=y)

# Fitting the data to the model
fit9 <- sampling(Stan_TypeII,
                 data = Stan_StochDataTypeII,
                 warmup = 200,
                 iter = 400,
                 chains = 2,
                 cores = 2,
                 thin = 1,
                 control = list(max_treedepth=15,
                                adapt_delta=0.99),
                 seed = 123,
                 check_data = TRUE,
                 diagnostic_file = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Fit9.csv",
                 show_messages = TRUE,
                 verbose = TRUE); beep(1)

FR_Parms <- c("r","O","h","b","c","u")

#### Tracking progress: Fit 9 ####
## 400 iter
## 31 min
## Warnings: 3 div trans; bulk ESS too low; tail ESS too low
## Diagnostics
Fit9_Print <- print(fit9, FR_Parms, 
                    probs=c(0.1, 0.5, 0.9), digits = 3)
# Package: rstan
Fit9_Trace <- stan_trace(fit9,FR_Parms)
# Package: bayesplot
color_scheme_set("pink")
Fit9_Dens <- mcmc_dens(fit9,FR_Parms)
Fit9_Overlay <- mcmc_dens_overlay(fit9,FR_Parms)
Fit9_Violin <- mcmc_violin(fit9, pars = FR_Parms, probs = c(0.1, 0.5, 0.9))
Fit9_Pairs <- mcmc_pairs(fit9,pars = FR_Parms)

Fit9_Trace
Fit9_Dens
Fit9_Overlay
Fit9_Violin
Fit9_Pairs

#### Plotting posterior ####

# Extract output info from fit9
Output9 <- rstan::extract(fit9,permuted=TRUE,include=TRUE)
# Get mean of 400 posterior draws for each timepoint in P and H
library(dplyr)
Z_means <- Z_df %>% summarise_each(funs(mean))
# Parse df by P and H
Z_means_P <- Z_means[-c(2923:5844)]
ncol(Z_means_P) # 2922; good
Z_means_H <- Z_means[-c(1:2922)]  
ncol(Z_means_H) # 2922
# Invert dfs
Z_means_P <- t(Z_means_P)
Z_means_H <- t(Z_means_H)
# Re-name columns 
colnames(Z_means_P) <- c("post_means_P")
colnames(Z_means_H) <- c("post_means_H")
head(y) # Recall that y is the data we gave the model (our stochastic data)
y_df <- data.frame(y)
colnames(y_df) <- c("P","H")
# Merge dfs
Mega_df <- cbind(Z_means_P,Z_means_H,y_df)
# Add time steps
ts <- seq(1:2922)
ts <- data.frame(ts)
# Final df
Mega_df <- cbind(Mega_df,ts)

library(ggplot2)
Post_By_Data_Plot <- ggplot(Mega_df,aes(x=ts))+
  geom_line(aes(y=P,colour="springgreen4"))+ # Data P
  geom_line(aes(y=H,colour="cornflowerblue"))+ # Data H
  geom_line(aes(y=post_means_P,colour="lightgreen"),size=1.5)+ # Post. P
  geom_line(aes(y=post_means_H,colour="lightskyblue1"),size=1.5)+ # Post. H
  scale_x_continuous(breaks=seq(0,3000,by=500))+ # Cleaning up the axes
  scale_y_continuous(breaks=seq(0,400,by=50))+
  scale_color_manual(name="Legend",
                     labels=c("Host immune cell abundance from data",
                              "Estimated parasite abundance (av. over 400 posterior draws)",
                              "Estimated host immune cell abundance (av. over 400 posterior draws)",
                              "Parasite abundance from data"),
                     values=c("cornflowerblue","lightgreen","lightskyblue1","springgreen4"))+
  xlab("Time Step")+
  ylab("Population Abundance")+
  ggtitle("Posterior Estimates Plotted Against Data")+
  theme_minimal()
Post_By_Data_Plot

y_rep <- Output9$y_rep
y_rep <- data.frame(y_rep)
y_rep_P <- y_rep[-c(2923:5844)]
y_rep_P <- as.matrix(y_rep_P)
y_rep_H <- y_rep[-c(1:2922)]
y_rep_H <- as.matrix(y_rep_H)
y_df <- as.data.frame(y)
y_P <- y_df[-c(2)]
y_P <- as.vector(y_P$V1,mode = "numeric")
y_H <- y_df[-c(1)]
y_H <- as.vector(y_H$V2,mode = "numeric")

color_scheme_set("green")
P_plot <- ppc_dens_overlay(y_P, y_rep_P[1:100,])+
  theme_minimal()+
  xlab("Abundance")+
  ylab("Density")+
  labs(title="Empirical vs. Estimated Distribution of Parasite Abundance Data",
       subtitle="Showing 100 Draws from the Posterior")
P_plot

color_scheme_set("brightblue")
H_plot <- ppc_dens_overlay(y_H, y_rep_H[1:100,])+
  theme_minimal()+
  xlab("Abundance")+
  ylab("Density")+
  labs(title="Empirical vs. Estimated Distribution of Host Immune Cell Abundance Data",
       subtitle="Showing 100 Draws from the Posterior")
H_plot