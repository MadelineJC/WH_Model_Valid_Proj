## Packages ##
library(rstan)
library(gdata)
library(bayesplot)
library(tidyverse)
library(brms)

## Fit object ##
Fit_NM_T2_Osc

## Basic model checks ##
Fit_NM_T2_Osc_Summ <- print(Fit_NM_T2_Osc, pars=c("r", "b", "c", "u", "sigma", "z_init"),
                            probs=c(0.1, 0.5, 0.9), digits = 3)
parms <- c("r", "b", "c", "u")
orange_scheme <- c("#ffebcc", "#ffcc80",
                   "#ffad33", "#e68a00",
                   "#995c00", "#663d00")
color_scheme_set(orange_scheme)
Output <- rstan::extract(Fit_NM_T2_Osc, permuted=TRUE, include=TRUE)
fit_Trace <- stan_trace(Fit_NM_T2_Osc, parms); fit_Trace
fit_Pairs <- mcmc_pairs(Fit_NM_T2_Osc, parms); fit_Pairs
fit_Dens <- mcmc_dens(Fit_NM_T2_Osc, parms); fit_Dens
fit_Overlay <- mcmc_dens_overlay(Fit_NM_T2_Osc, parms); fit_Overlay

## Compare data to estimates of abundance time series ##
y_rep <- Output$y_rep
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
StanData <- read_csv(here('./02_Data/NowakMay_StochSim_T2-Osc.csv'))
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
  scale_x_continuous(breaks = seq(0, 20, by = 5)) + # Cleaning up the axes
  scale_y_continuous(breaks = seq(0, 250, by = 50)) +
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
