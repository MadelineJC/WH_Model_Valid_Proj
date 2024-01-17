# Type III Model; Type II Data

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

    real dP_dt = P*r - H*(O*(P^2)/(1 + O*(P^2)*h)); // Mechanistic model; TypeIII
    real dH_dt = b + H*(c*(O*(P^2)/(1 + O*(P^2)*h))-u);
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

model{ // Ignore the means/variances for now...
  r~normal(2.5,1); // r
  O~beta(2,2); // O; bounded between 0 and 1
  h~normal(0.06,1); // h
  b~normal(35,1); // b
  c~normal(0.2,1); // c
  u~beta(2,2); // u; bounded between 0 and 1
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
StochDataTypeII = read.csv('/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/20Apr_StochData.csv')
N <- length(StochDataTypeII$t)-1 # df is 2923; N is 2922
ts <- 1:N
y_init <- c(10,350) # Initial states, P = 10; H = 350
y0 <- array(y_init) # For the ODE solver
t0 <- 0.1
y <- as.matrix(StochDataTypeII[2:(N+1),2:3])
y <- cbind(y[ ,1],y[ ,2]); # This worked, sick; where y[,1] is P, and y[,2] is H
Stan_StochDataTypeII <- list(N=N,ts=ts,y0=y0,t0=t0,y=y)

# Fitting the data to the model
Stan_Model_TypeII <- stan_model("Stan_Model_TypeII.stan")
fit.3M2D_2 <- sampling(Stan_TypeII,
                 data = Stan_StochDataTypeII,
                 warmup = 200,
                 iter = 400,
                 chains = 2,
                 cores = 2,
                 thin = 1,
                 # control = list(max_treedepth=15,
                 # adapt_delta=0.99),
                 seed = 1234,
                 check_data = TRUE,
                 diagnostic_file = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Fit.3M2D_2.csv",
                 show_messages = TRUE,
                 verbose = TRUE); beep(3)

FR_Parms <- c("r","O","h","b","c","u")
Fit.3M2D_Print <- print(fit.3M2D, FR_Parms, 
                    probs=c(0.1, 0.5, 0.9), digits = 3)
Output.3M2D <- rstan::extract(fit.3M2D,permuted=TRUE,include=TRUE)
y_rep <- Output.3M2D$y_rep
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

color_scheme_set("pink")
fit.3M2D_2_Trace <- stan_trace(fit.3M2D_2,FR_Parms)
fit.3M2D_2_Dens <- mcmc_dens(fit.3M2D_2,FR_Parms)
fit.3M2D_2_Overlay <- mcmc_dens_overlay(fit.3M2D_2,FR_Parms)
fit.3M2D_2_Violin <- mcmc_violin(fit.3M2D_2, pars = FR_Parms, probs = c(0.1, 0.5, 0.9))
fit.3M2D_2_Pairs <- mcmc_pairs(fit.3M2D_2,pars = FR_Parms)

fit.3M2D_2_Trace
fit.3M2D_2_Dens
fit.3M2D_2_Overlay
fit.3M2D_2_Violin
fit.3M2D_2_Pairs




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

    real dP_dt = P*r - H*(O*(P^2)/(1 + O*(P^2)*h)); // Mechanistic model; TypeIII
    real dH_dt = b + H*(c*(O*(P^2)/(1 + O*(P^2)*h))-u);
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
  real<lower=0,upper=1>u;
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

BC_Stoch_Data = read.csv('/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/BC_Stoch_Data.csv')
N <- length(BC_Stoch_Data$t)-1 # df is 871; N is 870
ts <- 1:N
y_init <- c(80,200) 
y <- as.matrix(BC_Stoch_Data[2:(N+1),2:3])
y <- cbind(y[ ,1],y[ ,2]); 
BC_Stan_Data <- list(N=N,ts=ts,y_init=y_init,y=y)

fit.3M2D_3_2 <- sampling(BC_2,
                  data = BC_Stan_Data,
                  warmup = 200,
                  iter = 400,
                  chains = 2,
                  cores = 2,
                  thin = 1,
                  # control = list(max_treedepth=15,
                                 # adapt_delta=0.99),
                  seed = 123,
                  check_data = TRUE,
                  diagnostic_file = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Fit.3M2D_3_2.csv",
                  show_messages = TRUE,
                  verbose = TRUE); beep(3)

Fit.3M2D_3_Print <- print(fit.3M2D_3, FR_Parms, 
                        probs=c(0.1, 0.5, 0.9), digits = 3)
Output.3M2D_3 <- rstan::extract(fit.3M2D_3,permuted=TRUE,include=TRUE)
y_rep <- Output.3M2D_3$y_rep
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

color_scheme_set("pink")
fit.3M2D_3_2_Trace <- stan_trace(fit.3M2D_3_2,FR_Parms)
fit.3M2D_3_2_Dens <- mcmc_dens(fit.3M2D_3_2,FR_Parms)
fit.3M2D_3_2_Overlay <- mcmc_dens_overlay(fit.3M2D_3_2,FR_Parms)
fit.3M2D_3_2_Violin <- mcmc_violin(fit.3M2D_3_2, pars = FR_Parms, probs = c(0.1, 0.5, 0.9))
fit.3M2D_3_2_Pairs <- mcmc_pairs(fit.3M2D_3_2,pars = FR_Parms)

fit.3M2D_3_2_Trace
fit.3M2D_3_2_Dens
fit.3M2D_3_2_Overlay
fit.3M2D_3_2_Violin
fit.3M2D_3_2_Pairs
