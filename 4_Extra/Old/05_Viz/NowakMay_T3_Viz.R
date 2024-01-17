## Packages ##
library(rstan)
library(gdata)
library(bayesplot)
library(tidyverse)
library(brms)

## Fit object ##
Fit_NM_T3

## Basic model checks ##
Fit_NM_T3_Summ <- print(Fit_NM_T3, pars=c("r", "b", "c", "u", "sigma"),
                            probs=c(0.1, 0.5, 0.9), digits = 3)
parms <- c("r", "b", "c", "u")
orange_scheme <- c("#ffebcc", "#ffcc80",
                            "#ffad33", "#e68a00",
                            "#995c00", "#663d00")
color_scheme_set(orange_scheme)
Output <- rstan::extract(Fit_NM_T3, permuted=TRUE, include=TRUE)
fit_Trace <- stan_trace(Fit_NM_T3, parms); fit_Trace
fit_Pairs <- mcmc_pairs(Fit_NM_T3, parms); fit_Pairs
fit_Dens <- mcmc_dens(Fit_NM_T3, parms); fit_Dens
fit_Overlay <- mcmc_dens_overlay(Fit_NM_T3, parms); fit_Overlay

## Compare data to estimates of abundance time series ##
Z_df <- data.frame(y_rep)
Z_means = Z_df %>% summarise(across(.cols = everything(),
                                    .fns = mean)) ##NOTE: summarise_each_ is deprecated
# Parse df by P and H
ncol(Z_means); x <- ncol(Z_means)/2
Z_means_P <- Z_means[-c(x + 1:ncol(Z_means))]
ncol(Z_means_P) 
Z_means_H <- Z_means[-c(1:x)]  
ncol(Z_means_H)
# Invert dfs
Z_means_P <- t(Z_means_P)
Z_means_H <- t(Z_means_H)
# Re-name columns 
colnames(Z_means_P) <- c("post_means_P")
colnames(Z_means_H) <- c("post_means_H")

## "Real data" ##
StanData <- read.csv('02_Data/NowakMay_StochSim_T3.csv')
N <- length(StanData$Time) - 1
ts <- 1:N
y_init <- c(StanData$Parasites[1], StanData$ImmuneCells[1])
y <- as.matrix(StanData[2:(N + 1), 2:3])
head(y) # Recall that y is the data we gave the model (our stochastic data)
y_df <- data.frame(y)
colnames(y_df) <- c("P","H")

# Merge dfs
Mega_df <- cbind(Z_means_P,Z_means_H,y_df)
# Add time steps
ts <- seq(nrow(Mega_df))
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
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Post_By_Data_Plot

## And now without labels for clarity ##
Post_By_Data_Plot <- ggplot(Mega_df,aes(x=ts))+
  geom_line(aes(y=P,colour="springgreen4"))+ # Data P
  geom_line(aes(y=H,colour="cornflowerblue"))+ # Data H
  geom_line(aes(y=post_means_P,colour="lightgreen"),size=1.5)+ # Post. P
  geom_line(aes(y=post_means_H,colour="lightskyblue1"),size=1.5)+ # Post. H
  scale_x_continuous(breaks = seq(0, 100, by = 25)) + # Cleaning up the axes
  scale_y_continuous(breaks = seq(0, 500, by = 100)) +
  scale_color_manual(name="Legend",
                     labels=c("Host immune cell abundance from data",
                              "Estimated parasite abundance (av. over 2000 posterior draws)",
                              "Estimated host immune cell abundance (av. over 2000 posterior draws)",
                              "Parasite abundance from data"),
                     values=c("cornflowerblue","lightgreen","lightskyblue1","springgreen4"))+
  xlab(" ")+
  ylab(" ")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
Post_By_Data_Plot

#### AR ===================================================
#### ... Parms ============================================
## Fit object ##
Fit_NM_T3

## Basic model checks ##
Fit_NM_T3_Summ <- print(Fit_NM_T3, pars=c("r", "b", "c", "u", "sigma"),
                            probs=c(0.1, 0.5, 0.9), digits = 3)
parms <- c("r", "b", "c", "u")
color_scheme_set("purple")
Output <- rstan::extract(Fit_NM_T3, permuted=TRUE, include=TRUE)
fit_Trace <- stan_trace(Fit_NM_T3, parms); fit_Trace
fit_Pairs <- mcmc_pairs(Fit_NM_T3, parms); fit_Pairs
fit_Dens <- mcmc_dens(Fit_NM_T3, parms); fit_Dens
fit_Overlay <- mcmc_dens_overlay(Fit_NM_T3, parms); fit_Overlay

#### ... Data comparison ==================================
## Let's do some plotting to compare data to estimates ##
posterior <- rstan::extract(Fit_Antia)
r_med <- median(posterior$r); K_med <- median(posterior$k); p_med <- median(posterior$p); o_med <- median(posterior$o)
r_low <- quantile(posterior$r, 0.025); K_low <- quantile(posterior$k, 0.025); p_low <- quantile(posterior$p, 0.025); o_low <- quantile(posterior$o, 0.025)
r_high <- quantile(posterior$r, 0.975); K_high <- quantile(posterior$k, 0.975); p_high <- quantile(posterior$p, 0.975); o_high <- quantile(posterior$o, 0.975)

## Med estimates ##
Antia <- function(t,y,p){
  r <- p[1]; K <- p[2]; b <- p[3]; o <- p[4]
  P <- y[1]; I <- y[2]
  dP = r*P - K*P*I
  dI = b*I*(P/(P+o))
  list(c(dP, dI))
}
r <- r_med; b <- b_med; c <- c_med; u <- u_med
parms <- c(r, b, c, u)
P0 <- 1; I0 <- 1
SysInit <- c(P0, I0)
TT <- seq(1,20,0.1) 
results <- lsoda(SysInit, TT, Antia, parms)
N <- results[,2]
Antia_Med = data.frame(results); colnames(Antia_Med) <- c("Time", "P", "I")

## Low estimates ##
Antia <- function(t,y,p){
  r <- p[1]; K <- p[2]; b <- p[3]; o <- p[4]
  P <- y[1]; I <- y[2]
  dP = r*P - K*P*I
  dI = b*I*(P/(P+o))
  list(c(dP, dI))
}
r <- r_low; K <- K_low; b <- p_low; o <- o_low
parms <- c(r, K, b, o)
P0 <- 1; I0 <- 1
SysInit <- c(P0, I0)
TT <- seq(1,20,0.1) 
results <- lsoda(SysInit, TT, Antia, parms)
N <- results[,2]
Antia_Low = data.frame(results); colnames(Antia_Low) <- c("Time", "P", "I")

## High estimates ##
Antia <- function(t,y,p){
  r <- p[1]; K <- p[2]; b <- p[3]; o <- p[4]
  P <- y[1]; I <- y[2]
  dP = r*P - K*P*I
  dI = b*I*(P/(P+o))
  list(c(dP, dI))
}
r <- r_high; K <- K_high; b <- p_high; o <- o_high
parms <- c(r, K, b, o)
P0 <- 1; I0 <- 1
SysInit <- c(P0, I0)
TT <- seq(1,20,0.1) 
results <- lsoda(SysInit, TT, Antia, parms)
N <- results[,2]
Antia_High = data.frame(results); colnames(Antia_High) <- c("Time", "P", "I")

DetDF <- data.frame(Antia_Det); colnames(DetDF) <- c("Time", "P_Det",  "I_Det")
DataDF <- data.frame(Antia_Stoch); colnames(DataDF) <- c("Time", "P_Stoch",  "I_Stoch")
EstDF <- data.frame(Antia_Med$Time, Antia_Med$P, Antia_Med$I, Antia_Low$P, Antia_Low$I, Antia_High$P, Antia_High$I); colnames(EstDF) <- c("Time", "Med_P_Est", "Med_I_Est", "Low_P_Est", "Low_I_Est", "High_P_Est", "High_I_Est")

## NOTE: I know this isn't a great way to plot the output/I'm not actually plotting the estimates, but here's a rough representation of what the model thinks is happening
Comp_Plot <- ggplot() +
  geom_line(data = DetDF, aes(x = Time, y = P_Det), linetype = "dotted") +
  geom_line(data = DetDF, aes(x = Time, y = I_Det), linetype = "dotted") +
  geom_line(data = DataDF, aes(x = Time, y = P_Stoch)) +
  geom_line(data = DataDF, aes(x = Time, y = I_Stoch)) +
  geom_line(data = EstDF, aes(x = Time, y = Med_P_Est), color = "springgreen4") +
  geom_line(data = EstDF, aes(x = Time, y = Med_I_Est), color = "springgreen4") +
  geom_line(data = EstDF, aes(x = Time, y = Low_P_Est), color = "royalblue1") +
  geom_line(data = EstDF, aes(x = Time, y = Low_I_Est), color = "royalblue1") +
  geom_line(data = EstDF, aes(x = Time, y = High_P_Est), color = "firebrick3") +
  geom_line(data = EstDF, aes(x = Time, y = High_I_Est), color = "firebrick3") +
  labs(x = " ", y = " ")+
  theme_minimal()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
Comp_Plot





