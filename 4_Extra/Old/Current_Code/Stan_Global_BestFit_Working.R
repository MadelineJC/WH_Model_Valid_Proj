#### First week of July, 2020 ####
# Here, I tried using data generated using an ODE solver that evaluates every time step instead of every 0.1 time step
# First, I tried just estimating r, b, c, and u, and that worked really well
# Then, I added O, which worked, and tried the same anaylsis with a logged prior on O

## NEXT: Get it to estimate h too

# Models
library(rstan)
library(deSolve)
library(GillespieSSA)
# Vis
library(gdata)
library(bayesplot)
library(tidyverse)
library(brms)
# Bullshit
library(beepr)

#### Deterministic ####
July1 <- function(t,y,p){
  r <- p[1]; O <- p[2]; h <- p[3]
  b <- p[4]; c <- p[5]; u <- p[6]
  P <- y[1]; H <- y[2]
  dP = P*r - H*(O*P/(1 + O*P*h))
  dH = b + H*(c*(O*P/(1 + O*P*h)) - u)
  list(c(dP,dH))
}
r <- 2.5; O <- 0.012; h <- 0.075
b <- 35; c <- 0.3; u <- 0.41
parms <- c(r,O,h,b,c,u)
P0 <- 80; H0 <- 200 
N0 <- c(P0,H0)
TT <- seq(1,20,1) 
results <- lsoda(N0,TT,July1,parms)
P <- results[,2]; H <- results[,3];
July1_DF = data.frame(results)
July1_Plot <- ggplot(July1_DF,aes(x=time))+
  geom_line(aes(y=X1, color="cornflowerblue"))+ # Parasite
  geom_line(aes(y=X2,color="springgreen4")) # Host
July1_Plot

#### Stochastic data gen ####
x0 <- c(P=80,H=200) 
a <- c("P*r",
       "H*(O*P/1 + O*P*h)", 
       "b + H*c*(O*P/1 + O*P*h)",
       "H*u")
nu <- matrix(c(+1,-1,0,0,
               0,0,+1,-1), nrow = 2, byrow = TRUE)
r = 2.5; O = 0.012; h = 0.075; b = 35; c = 0.3; u = 0.41
parms1 <- c(r = r, O = O, h = h, b = b, c = c, u = u)
tf = 20
method <- "OTL"
simName <- "July1"
set.seed(1508)
July1_SLC <- suppressWarnings(ssa(x0, a, nu, parms1, tf, method, simName,
                                  verbose = FALSE, 
                                  consoleInterval = 1, 
                                  censusInterval = 1, 
                                  maxWallTime = 30, 
                                  ignoreNegativeState = TRUE)) 
NewData <- July1_SLC$data 
NewData <- as.data.frame(NewData)
write.csv(NewData,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Data/NewData.csv",row.names = FALSE)
July1_Plot <- ggplot(NewData,aes(x=t))+
  geom_line(aes(y=P, color="cornflowerblue"))+ 
  geom_line(aes(y=H,color="springgreen4"))
July1_Plot

#### fit: Estimate O, and give the model h ####

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
  theta[{2}] ~ normal(0.5,0.5);
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
"GlobalFit_Stan_Working_1.stan")

model <- stan_model("GlobalFit_Stan_Working_1.stan")

N <- length(NewData$t) - 1
ts <- 1:N
y_init <- c(NewData$P[1], NewData$H[1])
y <- as.matrix(NewData[2:(N + 1), 2:3])
y <- cbind(y[ , 1], y[ , 2]); 
Data <- list(N = N, ts = ts, y_init = y_init, y = y)

# Fitting
fit <- sampling(model, data = Data, seed = 12); beep(3)

fit_summ <- print(fit, pars=c("theta", "sigma", "z_init"),
                  probs=c(0.1, 0.5, 0.9), digits = 3)
parms <- c("theta[1]","theta[2]","theta[3]","theta[4]","theta[5]")

#### fit: Model checks ####
Output1 <- rstan::extract(fit,permuted=TRUE,include=TRUE)
# Parms
fit_Trace <- stan_trace(fit,parms)
fit_Dens <- mcmc_dens(fit,parms)
fit_Overlay <- mcmc_dens_overlay(fit,parms)
fit_Violin <- mcmc_violin(fit,parms,probs = c(0.1, 0.5, 0.9))
fit_Pairs <- mcmc_pairs(fit,parms)
fit_Trace
fit_Dens
fit_Overlay
fit_Violin
fit_Pairs
# Data; density
y_rep <- Output1$y_rep
y_rep <- data.frame(y_rep)
length(y_rep)
y_rep_P <- y_rep[-c(21:40)]
y_rep_P <- as.matrix(y_rep_P)
y_rep_H <- y_rep[-c(1:20)]
y_rep_H <- as.matrix(y_rep_H)
y_df <- as.data.frame(y)
y_P <- y_df[-c(2)]
y_P <- as.vector(y_P$V1,mode = "numeric")
y_H <- y_df[-c(1)]
y_H <- as.vector(y_H$V2,mode = "numeric")
color_scheme_set("green")
P_plot <- ppc_dens_overlay(y_P, y_rep_P[1:200,])+
  theme_minimal()+
  xlab("Abundance")+
  ylab("Density")+
  labs(title="Empirical vs. Estimated Distribution of Parasite Abundance Data",
       subtitle="Showing 200 Draws from the Posterior")
P_plot # How estimates line up with data for P
color_scheme_set("brightblue")
H_plot <- ppc_dens_overlay(y_H, y_rep_H[1:200,])+
  theme_minimal()+
  xlab("Abundance")+
  ylab("Density")+
  labs(title="Empirical vs. Estimated Distribution of Host Immune Cell Abundance Data",
       subtitle="Showing 200 Draws from the Posterior")
H_plot # How estimates line up with data for H
# Data; time series
Z_df <- data.frame(y_rep)
Z_means <- Z_df %>% summarise_each(funs(mean))
# Parse df by P and H
Z_means_P <- Z_means[-c(21:40)]
ncol(Z_means_P) # 2922; good
Z_means_H <- Z_means[-c(1:20)]  
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
ts <- seq(1:20)
ts <- data.frame(ts)
# Final df
Mega_df <- cbind(Mega_df,ts)
Post_By_Data_Plot <- ggplot(Mega_df,aes(x=ts))+
  geom_line(aes(y=P,colour="springgreen4"))+ # Data P
  geom_line(aes(y=H,colour="cornflowerblue"))+ # Data H
  geom_line(aes(y=post_means_P,colour="lightgreen"),size=1.5)+ # Post. P
  geom_line(aes(y=post_means_H,colour="lightskyblue1"),size=1.5)+ # Post. H
  scale_x_continuous(breaks=seq(0,3000,by=500))+ # Cleaning up the axes
  scale_y_continuous(breaks=seq(0,400,by=50))+
  scale_color_manual(name="Legend",
                     labels=c("Host immune cell abundance from data",
                              "Estimated parasite abundance (av. over 2000 posterior draws)",
                              "Estimated host immune cell abundance (av. over 2000 posterior draws)",
                              "Parasite abundance from data"),
                     values=c("cornflowerblue","lightgreen","lightskyblue1","springgreen4"))+
  xlab("Time Step")+
  ylab("Population Abundance")+
  ggtitle("Posterior Estimates Plotted Against Data")+
  theme_minimal()
Post_By_Data_Plot
# Okay so this worked. It captured the oscillations

#### fit2: Now, with logged prior on O, and still giving the model h ####

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
  //target+=-log(theta[{2}]);
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
"GlobalFit_Stan_Working_2.stan")

model2 <- stan_model("GlobalFit_Stan_Working_2.stan")

N <- length(NewData$t) - 1
ts <- 1:N
y_init <- c(NewData$P[1], NewData$H[1])
y <- as.matrix(NewData[2:(N + 1), 2:3])
y <- cbind(y[ , 1], y[ , 2]); 
Data <- list(N = N, ts = ts, y_init = y_init, y = y)

# Fitting
fit2 <- sampling(model2, data = Data, seed = 12); beep(3)

fit_summ <- print(fit2, pars=c("theta", "sigma", "z_init"),
                  probs=c(0.1, 0.5, 0.9), digits = 3)
parms <- c("theta[1]","theta[2]","theta[3]","theta[4]","theta[5]")

#### fit2: Model checks ####

Output2 <- rstan::extract(fit2,permuted=TRUE,include=TRUE)
# Parms
fit2_Trace <- stan_trace(fit2,pars)
fit2_Dens <- mcmc_dens(fit2,parms)
fit2_Overlay <- mcmc_dens_overlay(fit2,parms)
fit2_Violin <- mcmc_violin(fit2,parms,probs = c(0.1, 0.5, 0.9))
fit2_Pairs <- mcmc_pairs(fit2,parms)
fit2_Trace
fit2_Dens
fit2_Overlay
fit2_Violin
fit2_Pairs
# Data; density
y_rep <- Output2$y_rep
y_rep <- data.frame(y_rep)
length(y_rep)
y_rep_P <- y_rep[-c(21:40)]
y_rep_P <- as.matrix(y_rep_P)
y_rep_H <- y_rep[-c(1:20)]
y_rep_H <- as.matrix(y_rep_H)
y_df <- as.data.frame(y)
y_P <- y_df[-c(2)]
y_P <- as.vector(y_P$V1,mode = "numeric")
y_H <- y_df[-c(1)]
y_H <- as.vector(y_H$V2,mode = "numeric")
color_scheme_set("green")
P_plot <- ppc_dens_overlay(y_P, y_rep_P[1:200,])+
  theme_minimal()+
  xlab("Abundance")+
  ylab("Density")+
  labs(title="Empirical vs. Estimated Distribution of Parasite Abundance Data",
       subtitle="Showing 200 Draws from the Posterior")
P_plot # How estimates line up with data for P
color_scheme_set("brightblue")
H_plot <- ppc_dens_overlay(y_H, y_rep_H[1:200,])+
  theme_minimal()+
  xlab("Abundance")+
  ylab("Density")+
  labs(title="Empirical vs. Estimated Distribution of Host Immune Cell Abundance Data",
       subtitle="Showing 200 Draws from the Posterior")
H_plot # How estimates line up with data for H
# Data; time series
Z_df <- data.frame(y_rep)
Z_means <- Z_df %>% summarise_each(funs(mean))
# Parse df by P and H
Z_means_P <- Z_means[-c(21:40)]
ncol(Z_means_P) # 2922; good
Z_means_H <- Z_means[-c(1:20)]  
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
ts <- seq(1:20)
ts <- data.frame(ts)
# Final df
Mega_df <- cbind(Mega_df,ts)
Post_By_Data_Plot <- ggplot(Mega_df,aes(x=ts))+
  geom_line(aes(y=P,colour="springgreen4"))+ # Data P
  geom_line(aes(y=H,colour="cornflowerblue"))+ # Data H
  geom_line(aes(y=post_means_P,colour="lightgreen"),size=1.5)+ # Post. P
  geom_line(aes(y=post_means_H,colour="lightskyblue1"),size=1.5)+ # Post. H
  scale_x_continuous(breaks=seq(0,3000,by=500))+ # Cleaning up the axes
  scale_y_continuous(breaks=seq(0,400,by=50))+
  scale_color_manual(name="Legend",
                     labels=c("Host immune cell abundance from data",
                              "Estimated parasite abundance (av. over 2000 posterior draws)",
                              "Estimated host immune cell abundance (av. over 2000 posterior draws)",
                              "Parasite abundance from data"),
                     values=c("cornflowerblue","lightgreen","lightskyblue1","springgreen4"))+
  xlab("Time Step")+
  ylab("Population Abundance")+
  ggtitle("Posterior Estimates Plotted Against Data")+
  theme_minimal()
Post_By_Data_Plot

#### fit3: Estimate O and h ####
### Things I've tried here that haven't worked:
  # (1) Estimating O and h separately, both with logged priors
  # (2) Estimating O and h together, with a logged prior
  # (3) Giving the model O and estimating h, with a logged prior
  # (4) Estimating O and h separately, with a logged prior on h ONLY
  # (5) Estimating O and h together, with a logged prior on Oh, but not O
  # (6) Estimating O and h together, with a logged prior, extremely tight on h
### Yet to try:
  

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
    = integrate_ode_bdf(dz_dt, z_init, 0, ts, theta,
                         rep_array(0.0, 0), rep_array(0, 0),
                         1e-5, 1e-3, 5e2); 
}
model {
  theta[{1}] ~ normal(2.5, 0.1);
  log(theta[{2}]) ~ normal(0.012,0.1);
  log(theta[{3}]) ~ normal(0.075,0.1);
  theta[{4}] ~ normal(35,0.1);
  theta[{5, 6}] ~ normal(0.35,0.5);
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
"GlobalFit_Stan_Working_3.stan")

model3 <- stan_model("GlobalFit_Stan_Working_3.stan")

N <- length(NewData$t) - 1
ts <- 1:N
y_init <- c(NewData$P[1], NewData$H[1])
y <- as.matrix(NewData[2:(N + 1), 2:3])
y <- cbind(y[ , 1], y[ , 2]); 
Data <- list(N = N, ts = ts, y_init = y_init, y = y)

# Fitting
fit3 <- sampling(model3, data = Data, 
                 chains = 1, iter = 4000, warmup = 2000, cores = 1,
                 control = list(max_treedepth=15,adapt_delta=0.99)); beep(3)

fit_summ <- print(fit3, pars=c("theta", "sigma", "z_init"),
                  probs=c(0.1, 0.5, 0.9), digits = 3)
parms <- c("theta[1]","theta[2]","theta[3]","theta[4]","theta[5]")

#### fit3: Model checks ####

Output3 <- rstan::extract(fit3,permuted=TRUE,include=TRUE)
# Parms
fit3_Trace <- stan_trace(fit3,parms)
fit3_Dens <- mcmc_dens(fit3,parms)
fit3_Overlay <- mcmc_dens_overlay(fit3,parms)
fit3_Violin <- mcmc_violin(fit3,parms,probs = c(0.1, 0.5, 0.9))
fit3_Pairs <- mcmc_pairs(fit3,parms)
fit3_Trace
fit3_Dens
fit3_Overlay
fit3_Violin
fit3_Pairs
# Data; density
y_rep <- Output3$y_rep
y_rep <- data.frame(y_rep)
length(y_rep)
y_rep_P <- y_rep[-c(21:40)]
y_rep_P <- as.matrix(y_rep_P)
y_rep_H <- y_rep[-c(1:20)]
y_rep_H <- as.matrix(y_rep_H)
y_df <- as.data.frame(y)
y_P <- y_df[-c(2)]
y_P <- as.vector(y_P$V1,mode = "numeric")
y_H <- y_df[-c(1)]
y_H <- as.vector(y_H$V2,mode = "numeric")
color_scheme_set("green")
P_plot <- ppc_dens_overlay(y_P, y_rep_P[1:200,])+
  theme_minimal()+
  xlab("Abundance")+
  ylab("Density")+
  labs(title="Empirical vs. Estimated Distribution of Parasite Abundance Data",
       subtitle="Showing 200 Draws from the Posterior")
P_plot # How estimates line up with data for P
color_scheme_set("brightblue")
H_plot <- ppc_dens_overlay(y_H, y_rep_H[1:200,])+
  theme_minimal()+
  xlab("Abundance")+
  ylab("Density")+
  labs(title="Empirical vs. Estimated Distribution of Host Immune Cell Abundance Data",
       subtitle="Showing 200 Draws from the Posterior")
H_plot # How estimates line up with data for H
# Data; time series
Z_df <- data.frame(y_rep)
Z_means <- Z_df %>% summarise_each(funs(mean))
# Parse df by P and H
Z_means_P <- Z_means[-c(21:40)]
ncol(Z_means_P) # 2922; good
Z_means_H <- Z_means[-c(1:20)]  
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
ts <- seq(1:20)
ts <- data.frame(ts)
# Final df
Mega_df <- cbind(Mega_df,ts)
Post_By_Data_Plot <- ggplot(Mega_df,aes(x=ts))+
  geom_line(aes(y=P,colour="springgreen4"))+ # Data P
  geom_line(aes(y=H,colour="cornflowerblue"))+ # Data H
  geom_line(aes(y=post_means_P,colour="lightgreen"),size=1.5)+ # Post. P
  geom_line(aes(y=post_means_H,colour="lightskyblue1"),size=1.5)+ # Post. H
  scale_x_continuous(breaks=seq(0,3000,by=500))+ # Cleaning up the axes
  scale_y_continuous(breaks=seq(0,400,by=50))+
  scale_color_manual(name="Legend",
                     labels=c("Host immune cell abundance from data",
                              "Estimated parasite abundance (av. over 2000 posterior draws)",
                              "Estimated host immune cell abundance (av. over 2000 posterior draws)",
                              "Parasite abundance from data"),
                     values=c("cornflowerblue","lightgreen","lightskyblue1","springgreen4"))+
  xlab("Time Step")+
  ylab("Population Abundance")+
  ggtitle("Posterior Estimates Plotted Against Data")+
  theme_minimal()
Post_By_Data_Plot