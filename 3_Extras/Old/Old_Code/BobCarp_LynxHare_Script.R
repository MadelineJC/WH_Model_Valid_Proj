#### Bob's Carpenter's work ####

# Data
Data <- read.csv('/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Data/BC_LH_Data.csv')
lynx_hare_df <- as.data.frame(Data)
summary(lynx_hare_df)
BC_Plot <- ggplot(lynx_hare_df,aes(x=Year))+
  geom_line(aes(y=Lynx, color="cornflowerblue"))+ # Parasite
  geom_line(aes(y=Hare,color="springgreen4")) # Host
BC_Plot

# Stat model
write("
functions {
  real[] dz_dt(real t,       // time
               real[] z,     // system state {prey, predator}
               real[] theta, // parameters
               real[] x_r,   // unused data
               int[] x_i) {
    real u = z[1];
    real v = z[2];

    real alpha = theta[1];
    real beta = theta[2];
    real gamma = theta[3];
    real delta = theta[4];

    real du_dt = (alpha - beta * v) * u;
    real dv_dt = (-gamma + delta * u) * v;
    return { du_dt, dv_dt };
  }
}
data {
  int<lower = 0> N;           // number of measurement times
  real ts[N];                 // measurement times > 0
  real y_init[2];             // initial measured populations
  real<lower = 0> y[N, 2];    // measured populations
}
parameters {
  real<lower = 0> theta[4];   // { alpha, beta, gamma, delta }
  real<lower = 0> z_init[2];  // initial population
  real<lower = 0> sigma[2];   // measurement errors
}
transformed parameters {
  real z[N, 2]
    = integrate_ode_rk45(dz_dt, z_init, 0, ts, theta,
                         rep_array(0.0, 0), rep_array(0, 0),
                         1e-5, 1e-3, 5e2);
}
model {
  theta[{1, 3}] ~ normal(1, 0.5);
  theta[{2, 4}] ~ normal(0.05, 0.05);
  sigma ~ lognormal(-1, 1);
  z_init ~ lognormal(log(10), 1);
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
"BC_LH_Model.stan")

model <- stan_model("BC_LH_Model.stan")

N <- length(lynx_hare_df$Year) - 1
ts <- 1:N
y_init <- c(lynx_hare_df$Hare[1], lynx_hare_df$Lynx[1])
y <- as.matrix(lynx_hare_df[2:(N + 1), 2:3])
y <- cbind(y[ , 2], y[ , 1]); # hare, lynx order
lynx_hare_data <- list(N = N, ts = ts, y_init = y_init, y = y)

# Fitting
fit <- sampling(model, data = lynx_hare_data, seed = 123); beep(3)

# Checking
print(fit, pars=c("theta", "sigma", "z_init"),
      probs=c(0.1, 0.5, 0.9), digits = 3)
## Output
#            mean se_mean    sd    10%    50%    90% n_eff  Rhat
#theta[1]   0.545   0.002 0.064  0.465  0.542  0.630  1076 1.002
#theta[2]   0.028   0.000 0.004  0.022  0.027  0.033  1195 1.001
#theta[3]   0.803   0.003 0.092  0.692  0.797  0.926   993 1.002
#theta[4]   0.024   0.000 0.004  0.020  0.024  0.029  1062 1.001
#sigma[1]   0.250   0.001 0.045  0.200  0.243  0.307  2537 1.001
#sigma[2]   0.252   0.001 0.044  0.200  0.245  0.309  2692 1.000
#z_init[1] 33.956   0.057 2.856 30.415 33.908 37.630  2474 1.000
#z_init[2]  5.933   0.012 0.535  5.273  5.912  6.614  2095 1.002

# Where true values for parms are:
# theta[1]: 0.55
# theta[2]: 0.028
# theta[3]: 0.84
# theta[4]: 0.026

Output <- rstan::extract(fit,permuted=TRUE,include=TRUE)
parms <- c("theta[1]","theta[2]","theta[3]","theta[4]","sigma[1]","sigma[2]","z_init[1]","z_init[2]")
# Parms
color_scheme_set("mix-blue-red")
fit_Trace <- stan_trace(fit,parms)
fit_Dens <- mcmc_dens(fit,parms)
fit_Overlay <- mcmc_dens_overlay(fit,parms)
fit_Violin <- mcmc_violin(fit,parms,probs = c(0.1, 0.5, 0.9))
fit_Pairs <- mcmc_pairs(fit,parms)
fit_Highlight_1 <- mcmc_trace_highlight(fit,pars=parms,highlight=1)
fit_Highlight_2 <- mcmc_trace_highlight(fit,pars=parms,highlight=2)
fit_Highlight_3 <- mcmc_trace_highlight(fit,pars=parms,highlight=3)
fit_Highlight_4 <- mcmc_trace_highlight(fit,pars=parms,highlight=4)
fit_Trace
fit_Dens
fit_Overlay
fit_Violin
fit_Pairs
fit_Highlight_1
fit_Highlight_2
fit_Highlight_3
fit_Highlight_4
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
  labs(title="Empirical vs. Estimated Distribution of Hare Abundance Data",
       subtitle="Showing 200 Draws from the Posterior")
P_plot # How estimates line up with data for P
color_scheme_set("brightblue")
H_plot <- ppc_dens_overlay(y_H, y_rep_H[1:200,])+
  theme_minimal()+
  xlab("Abundance")+
  ylab("Density")+
  labs(title="Empirical vs. Estimated Distribution of Lynx Abundance Data",
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