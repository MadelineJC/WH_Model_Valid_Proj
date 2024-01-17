#### DocString ####
# This code contains parts of the analysis regarding identification of functional
# responses in host-disease systems

# AUTHOR: Madeline Jarvis-Cross & Cole B. Brookson
# DATE OF CREATION: 2020-06-30
#### END ####

library(rstan)
library(deSolve)
library(GillespieSSA)
library(gdata)
library(bayesplot)
library(tidyverse)
library(brms)
library(beepr)
library(here)

#### GIVE: r ####

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
  log(theta[{1,2}]) ~ normal(0.5,0.5);
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
"Stan_Global_Model_r.stan")

model <- stan_model("Stan_Global_Model_r.stan")
NewData = read_csv('Data/NewData.csv')
N <- length(NewData$t) - 1
ts <- 1:N
y_init <- c(NewData$P[1], NewData$H[1])
y <- as.matrix(NewData[2:(N + 1), 2:3])
y <- cbind(y[ , 1], y[ , 2]); 
Data <- list(N = N, ts = ts, y_init = y_init, y = y)

# Fitting
fit_r <- sampling(model, data = Data, chains = 4, iter = 2000, cores = 2, seed = 12); beep(3)

fit_r_array <- as.array(fit_r)[,,-4] 
dimnames(fit_r_array)[[2]] <- c("Chain 1", "Chain 2", "Chain 3", "Chain 4")
dimnames(fit_r_array)[[3]][1:5] <- c("O = 0.012", "h = 0.075", "b = 35", "c = 0.3", "u = 0.41")
parms <- c("O = 0.012", "h = 0.075", "b = 35", "c = 0.3", "u = 0.41")

fit_r_summ <- print(fit_r, pars=c("theta", "sigma", "z_init"),
                    probs=c(0.1, 0.5, 0.9), digits = 3)

#### ... Model Checks ####
Output <- rstan::extract(fit_r,permuted=TRUE,include=TRUE)
# Parms
color_scheme_set("viridisA")
fit_r_Trace <- mcmc_trace(fit_r_array,parms)
# fit_r_Dens <- mcmc_dens(fit_r_array,parms)
fit_r_Overlay <- mcmc_dens_overlay(fit_r_array,parms)
# fit_r_Violin <- mcmc_violin(fit_r_array,parms,probs = c(0.1, 0.5, 0.9))
fit_r_Pairs <- mcmc_pairs(fit_r_array,parms)
fit_r_Trace
# fit_r_Dens
fit_r_Overlay
# fit_r_Violin
fit_r_Pairs
# Data; density
y_rep <- Output$y_rep
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
Z_means = Z_df %>% summarise(across(.cols = everything(),
                                    .fns = mean))##NOTE: summarisze_each_ is deprecated
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

#### GIVE: O ####

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
"Stan_Global_Model_O.stan")

model <- stan_model("Stan_Global_Model_O.stan")
NewData = read_csv('Data/NewData.csv')
N <- length(NewData$t) - 1
ts <- 1:N
y_init <- c(NewData$P[1], NewData$H[1])
y <- as.matrix(NewData[2:(N + 1), 2:3])
y <- cbind(y[ , 1], y[ , 2]); 
Data <- list(N = N, ts = ts, y_init = y_init, y = y)

# Fitting
fit_O <- sampling(model, data = Data, chains = 4, iter = 2000, cores = 2, seed = 12); beep(3)

fit_O_array <- as.array(fit_O)[,,-4] 
dimnames(fit_O_array)[[2]] <- c("Chain 1", "Chain 2", "Chain 3", "Chain 4")
dimnames(fit_O_array)[[3]][1:5] <- c("r = 2.5", "h = 0.075", "b = 35", "c = 0.3", "u = 0.41")
parms <- c("r = 2.5", "h = 0.075", "b = 35", "c = 0.3", "u = 0.41")

fit_O_summ <- print(fit_O, pars=c("theta", "sigma", "z_init"),
                    probs=c(0.1, 0.5, 0.9), digits = 3)

#### ... Model Checks ####
Output <- rstan::extract(fit_O,permuted=TRUE,include=TRUE)
# Parms
color_scheme_set("viridisA")
fit_O_Trace <- mcmc_trace(fit_O_array,parms)
# fit_O_Dens <- mcmc_dens(fit_O_array,parms)
fit_O_Overlay <- mcmc_dens_overlay(fit_O_array,parms)
# fit_O_Violin <- mcmc_violin(fit_O_array,parms,probs = c(0.1, 0.5, 0.9))
fit_O_Pairs <- mcmc_pairs(fit_O_array,parms)
fit_O_Trace
# fit_O_Dens
fit_O_Overlay
# fit_O_Violin
fit_O_Pairs
# Data; density
y_rep <- Output$y_rep
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
  xlab(" ")+
  ylab(" ")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# labs(title="Empirical vs. Estimated Distribution of Parasite Abundance Data",
#      subtitle="Showing 200 Draws from the Posterior")
P_plot # How estimates line up with data for P
color_scheme_set("brightblue")
H_plot <- ppc_dens_overlay(y_H, y_rep_H[1:200,])+
  theme_minimal()+
  xlab(" ")+
  ylab(" ")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# labs(title="Empirical vs. Estimated Distribution of Host Immune Cell Abundance Data",
#      subtitle="Showing 200 Draws from the Posterior")
H_plot # How estimates line up with data for H
# Data; time series
Z_df <- data.frame(y_rep)
Z_means = Z_df %>% summarise(across(.cols = everything(),
                                    .fns = mean))##NOTE: summarisze_each_ is deprecated
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
  geom_line(aes(y=P,colour="springgreen4"),size=1)+ # Data P
  geom_line(aes(y=H,colour="cornflowerblue"),size=1)+ # Data H
  geom_line(aes(y=post_means_P,colour="lightgreen"),size=1)+ # Post. P
  geom_line(aes(y=post_means_H,colour="lightskyblue1"),size=1)+ # Post. H
  scale_x_continuous(breaks=seq(0,3000,by=500))+ # Cleaning up the axes
  scale_y_continuous(breaks=seq(0,400,by=50))+
  scale_color_manual(name="Legend",
                     labels=c("Host immune cell abundance from data",
                              "Est. parasite abundance",
                              "Est. host immune cell abundance",
                              "Parasite abundance from data"),
                     values=c("cornflowerblue","lightgreen","lightskyblue1","springgreen4"),
                     guide = F)+
  xlab(" ")+
  ylab(" ")+
  # ggtitle("Posterior Estimates Plotted Against Data")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Post_By_Data_Plot

#### GIVE: h ####

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
"Stan_Global_Model_h.stan")

model <- stan_model("Stan_Global_Model_h.stan")
NewData = read_csv('Data/NewData.csv')
N <- length(NewData$t) - 1
ts <- 1:N
y_init <- c(NewData$P[1], NewData$H[1])
y <- as.matrix(NewData[2:(N + 1), 2:3])
y <- cbind(y[ , 1], y[ , 2]); 
Data <- list(N = N, ts = ts, y_init = y_init, y = y)

# Fitting
fit_h <- sampling(model, data = Data, chains = 4, iter = 2000, cores = 2, seed = 12); beep(3)

fit_h_array <- as.array(fit_h)[,,-4] 
dimnames(fit_h_array)[[2]] <- c("Chain 1", "Chain 2", "Chain 3", "Chain 4")
dimnames(fit_h_array)[[3]][1:5] <- c("r = 2.5", "O = 0.012", "b = 35", "c = 0.3", "u = 0.41")
parms <- c("r = 2.5", "O = 0.012", "b = 35", "c = 0.3", "u = 0.41")

fit_h_summ <- print(fit_h, pars=c("theta", "sigma", "z_init"),
                    probs=c(0.1, 0.5, 0.9), digits = 3)

#### ... Model Checks ####
Output <- rstan::extract(fit_h,permuted=TRUE,include=TRUE)
# Parms
color_scheme_set("viridisA")
fit_h_Trace <- mcmc_trace(fit_h_array,parms)
# fit_h_Dens <- mcmc_dens(fit_h_array,parms)
fit_h_Overlay <- mcmc_dens_overlay(fit_h_array,parms)
# fit_h_Violin <- mcmc_violin(fit_h_array,parms,probs = c(0.1, 0.5, 0.9))
fit_h_Pairs <- mcmc_pairs(fit_h_array,parms)
fit_h_Trace
# fit_h_Dens
fit_h_Overlay
# fit_h_Violin
fit_h_Pairs
# Data; density
y_rep <- Output$y_rep
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
  xlab(" ")+
  ylab(" ")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# labs(title="Empirical vs. Estimated Distribution of Parasite Abundance Data",
#      subtitle="Showing 200 Draws from the Posterior")
P_plot # How estimates line up with data for P
color_scheme_set("brightblue")
H_plot <- ppc_dens_overlay(y_H, y_rep_H[1:200,])+
  theme_minimal()+
  xlab(" ")+
  ylab(" ")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# labs(title="Empirical vs. Estimated Distribution of Host Immune Cell Abundance Data",
#      subtitle="Showing 200 Draws from the Posterior")
H_plot # How estimates line up with data for H
# Data; time series
Z_df <- data.frame(y_rep)
Z_means = Z_df %>% summarise(across(.cols = everything(),
                                    .fns = mean))##NOTE: summarisze_each_ is deprecated
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
  geom_line(aes(y=P,colour="springgreen4"),size=1)+ # Data P
  geom_line(aes(y=H,colour="cornflowerblue"),size=1)+ # Data H
  geom_line(aes(y=post_means_P,colour="lightgreen"),size=1)+ # Post. P
  geom_line(aes(y=post_means_H,colour="lightskyblue1"),size=1)+ # Post. H
  scale_x_continuous(breaks=seq(0,3000,by=500))+ # Cleaning up the axes
  scale_y_continuous(breaks=seq(0,400,by=50))+
  scale_color_manual(name="Legend",
                     labels=c("Host immune cell abundance from data",
                              "Est. parasite abundance",
                              "Est. host immune cell abundance",
                              "Parasite abundance from data"),
                     values=c("cornflowerblue","lightgreen","lightskyblue1","springgreen4"),
                     guide = F)+
  xlab(" ")+
  ylab(" ")+
  # ggtitle("Posterior Estimates Plotted Against Data")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Post_By_Data_Plot

#### GIVE: b ####

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
  log(theta[{2,3}]) ~ normal(0.5,0.5);
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
"Stan_Global_Model_b.stan")

model <- stan_model("Stan_Global_Model_b.stan")
NewData = read_csv('Data/NewData.csv')
N <- length(NewData$t) - 1
ts <- 1:N
y_init <- c(NewData$P[1], NewData$H[1])
y <- as.matrix(NewData[2:(N + 1), 2:3])
y <- cbind(y[ , 1], y[ , 2]); 
Data <- list(N = N, ts = ts, y_init = y_init, y = y)

# Fitting
fit_b <- sampling(model, data = Data, chains = 4, iter = 2000, cores = 2, seed = 12); beep(3)

fit_b_array <- as.array(fit_b)[,,-4] 
dimnames(fit_b_array)[[2]] <- c("Chain 1", "Chain 2", "Chain 3", "Chain 4")
dimnames(fit_b_array)[[3]][1:5] <- c("r = 2.5", "O = 0.012", "h = 0.075", "c = 0.3", "u = 0.41")
parms <- c("r = 2.5", "O = 0.012", "h = 0.075", "c = 0.3", "u = 0.41")

fit_b_summ <- print(fit_b, pars=c("theta", "sigma", "z_init"),
                    probs=c(0.1, 0.5, 0.9), digits = 3)

#### ... Model Checks ####
Output <- rstan::extract(fit_b,permuted=TRUE,include=TRUE)
# Parms
color_scheme_set("viridisA")
fit_b_Trace <- mcmc_trace(fit_b_array,parms)
# fit_b_Dens <- mcmc_dens(fit_b_array,parms)
fit_b_Overlay <- mcmc_dens_overlay(fit_b_array,parms)
# fit_b_Violin <- mcmc_violin(fit_b_array,parms,probs = c(0.1, 0.5, 0.9))
fit_b_Pairs <- mcmc_pairs(fit_b_array,parms)
fit_b_Trace
# fit_b_Dens
fit_b_Overlay
# fit_b_Violin
fit_b_Pairs
# Data; density
y_rep <- Output$y_rep
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
Z_means = Z_df %>% summarise(across(.cols = everything(),
                                    .fns = mean))##NOTE: summarisze_each_ is deprecated
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

#### GIVE: c ####

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
  log(theta[{2,3}]) ~ normal(0.5,0.5);
  theta[{4}] ~ normal(35,1);
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
"Stan_Global_Model_c.stan")

model <- stan_model("Stan_Global_Model_c.stan")
NewData = read_csv('Data/NewData.csv')
N <- length(NewData$t) - 1
ts <- 1:N
y_init <- c(NewData$P[1], NewData$H[1])
y <- as.matrix(NewData[2:(N + 1), 2:3])
y <- cbind(y[ , 1], y[ , 2]); 
Data <- list(N = N, ts = ts, y_init = y_init, y = y)

# Fitting
fit_c <- sampling(model, data = Data, chains = 4, iter = 2000, cores = 2, seed = 12); beep(3)

fit_c_array <- as.array(fit_c)[,,-4] 
dimnames(fit_c_array)[[2]] <- c("Chain 1", "Chain 2", "Chain 3", "Chain 4")
dimnames(fit_c_array)[[3]][1:5] <- c("r = 2.5", "O = 0.012", "h = 0.075", "b = 35", "u = 0.41")
parms <- c("r = 2.5", "O = 0.012", "h = 0.075", "b = 35", "u = 0.41")

fit_c_summ <- print(fit_c, pars=c("theta", "sigma", "z_init"),
                    probs=c(0.1, 0.5, 0.9), digits = 3)

#### ... Model Checks ####
Output <- rstan::extract(fit_c,permuted=TRUE,include=TRUE)
# Parms
color_scheme_set("viridisA")
fit_c_Trace <- mcmc_trace(fit_c_array,parms)
# fit_c_Dens <- mcmc_dens(fit_c_array,parms)
fit_c_Overlay <- mcmc_dens_overlay(fit_c_array,parms)
# fit_c_Violin <- mcmc_violin(fit_c_array,parms,probs = c(0.1, 0.5, 0.9))
fit_c_Pairs <- mcmc_pairs(fit_c_array,parms)
fit_c_Trace
# fit_c_Dens
fit_c_Overlay
# fit_c_Violin
fit_c_Pairs
# Data; density
y_rep <- Output$y_rep
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
Z_means = Z_df %>% summarise(across(.cols = everything(),
                                    .fns = mean))##NOTE: summarisze_each_ is deprecated
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

#### GIVE: u ####

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
  log(theta[{2,3}]) ~ normal(0.5,0.5);
  theta[{4}] ~ normal(35,1);
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
"Stan_Global_Model_u.stan")

model <- stan_model("Stan_Global_Model_u.stan")
NewData = read_csv('Data/NewData.csv')
N <- length(NewData$t) - 1
ts <- 1:N
y_init <- c(NewData$P[1], NewData$H[1])
y <- as.matrix(NewData[2:(N + 1), 2:3])
y <- cbind(y[ , 1], y[ , 2]); 
Data <- list(N = N, ts = ts, y_init = y_init, y = y)

# Fitting
fit_u <- sampling(model, data = Data, chains = 4, iter = 2000, cores = 2, seed = 12); beep(3)

fit_u_array <- as.array(fit_u)[,,-4] 
dimnames(fit_u_array)[[2]] <- c("Chain 1", "Chain 2", "Chain 3", "Chain 4")
dimnames(fit_u_array)[[3]][1:5] <- c("r = 2.5", "O = 0.012", "h = 0.075", "b = 35", "c = 0.3")
parms <- c("r = 2.5", "O = 0.012", "h = 0.075", "b = 35", "c = 0.3")

fit_u_summ <- print(fit_u, pars=c("theta", "sigma", "z_init"),
                    probs=c(0.1, 0.5, 0.9), digits = 3)

#### ... Model Checks ####
Output <- rstan::extract(fit_u,permuted=TRUE,include=TRUE)
# Parms
color_scheme_set("viridisA")
fit_u_Trace <- mcmc_trace(fit_u_array,parms)
# fit_u_Dens <- mcmc_dens(fit_u_array,parms)
fit_u_Overlay <- mcmc_dens_overlay(fit_u_array,parms)
# fit_u_Violin <- mcmc_violin(fit_u_array,parms,probs = c(0.1, 0.5, 0.9))
fit_u_Pairs <- mcmc_pairs(fit_u_array,parms)
fit_u_Trace
# fit_u_Dens
fit_u_Overlay
# fit_u_Violin
fit_u_Pairs
# Data; density
y_rep <- Output$y_rep
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
Z_means = Z_df %>% summarise(across(.cols = everything(),
                                    .fns = mean))##NOTE: summarisze_each_ is deprecated
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

#### GIVE: none ####

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
  log(theta[{2,3}]) ~ normal(0.5,0.5);
  theta[{4}] ~ normal(35,1);
  theta[{5}] ~ normal(0.5,0.5);
  theta[{6}] ~ normal(0.5,0.5);
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
"Stan_Global_Model.stan")

model <- stan_model("Stan_Global_Model.stan")
NewData = read_csv('Data/NewData.csv')
N <- length(NewData$t) - 1
ts <- 1:N
y_init <- c(NewData$P[1], NewData$H[1])
y <- as.matrix(NewData[2:(N + 1), 2:3])
y <- cbind(y[ , 1], y[ , 2]); 
Data <- list(N = N, ts = ts, y_init = y_init, y = y)

# Fitting
fit <- sampling(model, data = Data, chains = 4, iter = 2000, cores = 2, seed = 12); beep(3)

fit_array <- as.array(fit)[,,-4] 
dimnames(fit_array)[[2]] <- c("Chain 1", "Chain 2", "Chain 3", "Chain 4")
dimnames(fit_array)[[3]][1:6] <- c("r = 2.5", "O = 0.012", "h = 0.075", "b = 35", "c = 0.3", "u = 0.41")
parms <- c("r = 2.5", "O = 0.012", "h = 0.075", "b = 35", "c = 0.3", "u = 0.41")

fit_summ <- print(fit, pars=c("theta", "sigma", "z_init"),
                  probs=c(0.1, 0.5, 0.9), digits = 3)

#### ... Model Checks ####
Output <- rstan::extract(fit,permuted=TRUE,include=TRUE)
# Parms
color_scheme_set("viridisA")
fit_Trace <- mcmc_trace(fit_array,parms)
# fit_Dens <- mcmc_dens(fit_array,parms)
fit_Overlay <- mcmc_dens_overlay(fit_array,parms)
# fit_Violin <- mcmc_violin(fit_array,parms,probs = c(0.1, 0.5, 0.9))
fit_Pairs <- mcmc_pairs(fit_array,parms)
fit_Trace
# fit_Dens
fit_Overlay
# fit_Violin
fit_Pairs
# Data; density
y_rep <- Output$y_rep
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
Z_means = Z_df %>% summarise(across(.cols = everything(),
                                    .fns = mean))##NOTE: summarisze_each_ is deprecated
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

#### GIVE: O and h ####

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
  theta[{2}] ~ normal(35,1);
  theta[{3}] ~ normal(0.5,0.5);
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
"Stan_Global_Model_Oh.stan")

model <- stan_model("Stan_Global_Model_Oh.stan")
NewData = read_csv('Data/NewData.csv')
N <- length(NewData$t) - 1
ts <- 1:N
y_init <- c(NewData$P[1], NewData$H[1])
y <- as.matrix(NewData[2:(N + 1), 2:3])
y <- cbind(y[ , 1], y[ , 2]); 
Data <- list(N = N, ts = ts, y_init = y_init, y = y)

# Fitting
fit_Oh <- sampling(model, data = Data, chains = 4, iter = 2000, cores = 2, seed = 12); beep(3)

fit_Oh_array <- as.array(fit_Oh)[,,-4] 
dimnames(fit_Oh_array)[[2]] <- c("Chain 1", "Chain 2", "Chain 3", "Chain 4")
dimnames(fit_Oh_array)[[3]][1:4] <- c("r = 2.5", "b = 35", "c = 0.3", "u = 0.41")
parms <- c("r = 2.5", "b = 35", "c = 0.3", "u = 0.41")

fit_Oh_summ <- print(fit_Oh, pars=c("theta", "sigma", "z_init"),
                     probs=c(0.1, 0.5, 0.9), digits = 3)

#### ... Model Checks ####
Output <- rstan::extract(fit_Oh,permuted=TRUE,include=TRUE)
# Parms
color_scheme_set("viridisA")
fit_Oh_Trace <- mcmc_trace(fit_Oh_array,parms)
# fit_Oh_Dens <- mcmc_dens(fit_Oh_array,parms)
fit_Oh_Overlay <- mcmc_dens_overlay(fit_Oh_array,parms)
# fit_Oh_Violin <- mcmc_violin(fit_Oh_array,parms,probs = c(0.1, 0.5, 0.9))
fit_Oh_Pairs <- mcmc_pairs(fit_Oh_array,parms)
fit_Oh_Trace
# fit_Oh_Dens
fit_Oh_Overlay
# fit_Oh_Violin
fit_Oh_Pairs
# Data; density
y_rep <- Output$y_rep
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
  xlab(" ")+
  ylab(" ")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# labs(title="Empirical vs. Estimated Distribution of Parasite Abundance Data",
#      subtitle="Showing 200 Draws from the Posterior")
P_plot # How estimates line up with data for P
color_scheme_set("brightblue")
H_plot <- ppc_dens_overlay(y_H, y_rep_H[1:200,])+
  theme_minimal()+
  xlab(" ")+
  ylab(" ")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# labs(title="Empirical vs. Estimated Distribution of Host Immune Cell Abundance Data",
#      subtitle="Showing 200 Draws from the Posterior")
H_plot # How estimates line up with data for H
# Data; time series
Z_df <- data.frame(y_rep)
Z_means = Z_df %>% summarise(across(.cols = everything(),
                                    .fns = mean))##NOTE: summarisze_each_ is deprecated
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
  geom_line(aes(y=P,colour="springgreen4"),size=1)+ # Data P
  geom_line(aes(y=H,colour="cornflowerblue"),size=1)+ # Data H
  geom_line(aes(y=post_means_P,colour="lightgreen"),size=1)+ # Post. P
  geom_line(aes(y=post_means_H,colour="lightskyblue1"),size=1)+ # Post. H
  scale_x_continuous(breaks=seq(0,3000,by=500))+ # Cleaning up the axes
  scale_y_continuous(breaks=seq(0,400,by=50))+
  scale_color_manual(name="Legend",
                     labels=c("Host immune cell abundance from data",
                              "Est. parasite abundance",
                              "Est. host immune cell abundance",
                              "Parasite abundance from data"),
                     values=c("cornflowerblue","lightgreen","lightskyblue1","springgreen4"),
                     guide = F)+
  xlab(" ")+
  ylab(" ")+
  # ggtitle("Posterior Estimates Plotted Against Data")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Post_By_Data_Plot

#### Plotting priors ####
# theta[{1}] ~ normal(2.5, 1); r
# log(theta[{1,2}]) ~ normal(0.5,0.5); O, h
# theta[{3}] ~ normal(35,1); b
# theta[{4, 5}] ~ normal(0.5,0.5); c, u

r <- rnorm(1000, 2.5, 1); dr = density(r)
O <- rnorm(2000, 0.5, 0.5) 
O <- O[O >= 0]
O <- O[O <= 1]
dO = density(O)
h <- rnorm(1000, 0.5, 0.5); dh = density(h)
b <- rnorm(1000, 35, 1); db = density(b)
c <- rnorm(1000, 0.5, 0.5); dc = density(c)
u <- rnorm(2000, 0.5, 0.5)
u <- u[u >= 0]
u <- u[u <= 1]
du = density(u)

plot(dr, type = "l", xlim = c(0,50), ylim = c(0,1))
lines(dO)
lines(dh)
lines(db)
lines(dc)
lines(du)

