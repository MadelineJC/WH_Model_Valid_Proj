library(rstan)
library(gdata)
library(bayesplot)
library(tidyverse)
library(brms)
library(deSolve)
library(GillespieSSA)

#### Deterministic ####
BC <- function(t,y,p){
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
TT <- seq(400)
results <- lsoda(N0,TT,BC,parms)
P <- results[,2]; H <- results[,3];
BC_Det_DF = data.frame(results)
BC_Det_Plot <- ggplot(BC_Det_DF,aes(x=time))+
  geom_line(aes(y=X1, color="cornflowerblue"))+ # Parasite
  geom_line(aes(y=X2,color="springgreen4"))+ # Host
  labs(title="Solution 1: Type II Deterministic Time Series",
       x="Time",
       y="Population Abundance")+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"))+
  theme_minimal()
BC_Det_Plot 

#### Stochastic ####
x0 <- c(P=80,H=200) 
a <- c("P*r",
       "H*(O*P/1 + O*P*h)", 
       "b + H*c*(O*P/1 + O*P*h)",
       "H*u")
nu <- matrix(c(+1,-1,0,0,
               0,0,+1,-1), nrow = 2, byrow = TRUE)
r = 2.5; O = 0.012; h = 0.075; b = 35; c = 0.3; u = 0.41
parms1 <- c(r = r, O = O, h = h, b = b, c = c, u = u)
tf = 400
method <- "OTL"
simName <- "BC"
set.seed(1508)
BC_Results <- suppressWarnings(ssa(x0, a, nu, parms1, tf, method, simName,
                                   verbose = FALSE, 
                                   consoleInterval = 0, 
                                   censusInterval = 0.1, 
                                   maxWallTime = 30, 
                                   ignoreNegativeState = TRUE)) 
BC_Stoch <- BC_Results$data 
BC_Stoch <- as.data.frame(BC_Stoch)
BC_Stoch <- BC_Stoch[-c(2476:3493),]
write.csv(BC_Stoch,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/SLC_Data.csv",row.names = FALSE)
BC_Stoch_Plot <- ggplot(BC_Stoch,aes(x=t))+
  geom_line(aes(y=P, color="cornflowerblue"))+ 
  geom_line(aes(y=H,color="springgreen4"))+ 
  labs(title="Solution 1: Type II Stochastic Time Series",
       x="Time",
       y="Population Abundance")+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"))+
  theme_minimal()
BC_Stoch_Plot

#### Stan ####
write("
functions{
  real[] dZ_dt(
    real t, 
    real[] Z, 
    real[] theta, 
    real[] x_r, 
    int[] x_i){
    real P = Z[1]; 
    real H = Z[2];

    real r = theta[1]; 
    real O = theta[2];
    real h = theta[3];
    real b = theta[4];
    real c = theta[5];
    real u = theta[6];

    real dP_dt = P*r - H*(O*P/(1 + O*P*h));
    real dH_dt = b + H*(c*(O*P/(1 + O*P*h))-u);
    return({dP_dt,dH_dt});
  }
}
data {
  int<lower=0> N;           
  real ts[N];                 
  real y_init[2];             
  real<lower=0> y[N, 2];    
}
parameters {
  real<lower=0>r;
  real<lower=0,upper=1>O;
  real<lower=0>h;
  real<lower=0>b;
  real<lower=0,upper=1>c;
  real<lower=0>u;
  real<lower=0>Z_init[2];  
  real<lower=0>sigma[2];   
}
transformed parameters {
  real Z[N, 2]
    = integrate_ode_rk45(dZ_dt, Z_init, 0, ts, {r, O, h, b, c, u},
                         rep_array(0.0, 0), rep_array(0, 0),
                         1e-10, 1e-10, 5e2);
}
model {
  r~normal(2.5,1);
  O~normal(0.015,0.5); 
  h~normal(0.075,0.5); 
  b~normal(35,1); 
  c~normal(0.3,0.5); 
  u~normal(0.4,0.5);
  sigma~lognormal(-1, 1);
  Z_init~lognormal(log(3), 1);
  for (k in 1:2) {
    y_init[k] ~ lognormal(log(Z_init[k]), sigma[k]);
    y[ , k] ~ lognormal(log(Z[ , k]), sigma[k]);
  }
}
generated quantities {
  real y_init_rep[2];
  real y_rep[N, 2];
  for (k in 1:2) {
    y_init_rep[k] = lognormal_rng(log(Z_init[k]), sigma[k]);
    for (n in 1:N)
      y_rep[n, k] = lognormal_rng(log(Z[n, k]), sigma[k]);
  }
}",
"BC_2.stan")

stanc("BC_2.stan") 

BC_2 <- stan_model("BC_2.stan")

BC_Stoch_Data = read.csv('/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Data/June30.csv')
N <- length(BC_Stoch_Data$t)-1 # df is 871; N is 870
ts <- 1:N
y_init <- c(80,200) 
y <- as.matrix(BC_Stoch_Data[2:(N+1),2:3])
y <- cbind(y[ ,1],y[ ,2]); 
BC_Stan_Data <- list(N=N,ts=ts,y_init=y_init,y=y)

fitx <- sampling(BC_2,
                  data = BC_Stan_Data,
                  warmup = 200,
                  iter = 400,
                  chains = 1,
                  cores = 1,
                  thin = 1,
                  control = list(max_treedepth=15,
                                 adapt_delta=0.99),
                  seed = 123,
                  check_data = TRUE,
                  #diagnostic_file = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Fit_LN.SLC.csv",
                  show_messages = TRUE,
                  verbose = TRUE); beep(3)

# Visualization
FR_Parms <- c("r","O","h","b","c","u")
fit15_Print <- print(fitx, 
                     probs=c(0.1, 0.5, 0.9), digits = 3)
fit14_Trace <- stan_trace(fit14,FR_Parms)
color_scheme_set("pink")
fit14_Dens <- mcmc_dens(fit14,FR_Parms)
fit14_Overlay <- mcmc_dens_overlay(fit14,FR_Parms)
fit14_Violin <- mcmc_violin(fit14, pars = FR_Parms, probs = c(0.1, 0.5, 0.9))
fit14_Pairs <- mcmc_pairs(fit14,pars = FR_Parms)
fit14_Trace
fit14_Dens
fit14_Overlay
fit14_Violin
fit14_Pairs

# Extract output info from fit9
Output15 <- rstan::extract(fit15,permuted=TRUE,include=TRUE)
y_rep <- Output15$y_rep
Post <- Output15$Z
Post_df <- data.frame(Post)
# Get mean of 400 posterior draws for each timepoint in P and H
Post_means <- Post_df %>% summarise_each(funs(mean))
# Parse df by P and H
Post_P <- Post_means[-c(871:1740)]
ncol(Post_P) # 2922; good
Post_H <- Post_means[-c(1:870)]  
# Invert dfs
Post_P <- t(Post_P)
Post_H <- t(Post_H)
# Re-name columns 
colnames(Post_P) <- c("post_means_P")
colnames(Post_H) <- c("post_means_H")
head(y) # Recall that y is the data we gave the model (our stochastic data)
y_df <- data.frame(y)
colnames(y_df) <- c("P","H")
# Merge dfs
Mega_df <- cbind(Post_P,Post_H,y_df)
# Add time steps
ts <- seq(1:1740)
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
                              "Estimated parasite abundance (av. over 400 posterior draws)",
                              "Estimated host immune cell abundance (av. over 400 posterior draws)",
                              "Parasite abundance from data"),
                     values=c("cornflowerblue","lightgreen","lightskyblue1","springgreen4"))+
  xlab("Time Step")+
  ylab("Population Abundance")+
  ggtitle("Posterior Estimates Plotted Against Data")+
  theme_minimal()
Post_By_Data_Plot

y_rep <- data.frame(y_rep)
y_rep_P <- y_rep[-c(871:1740)]
y_rep_P <- as.matrix(y_rep_P)
y_rep_H <- y_rep[-c(1:870)]
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
