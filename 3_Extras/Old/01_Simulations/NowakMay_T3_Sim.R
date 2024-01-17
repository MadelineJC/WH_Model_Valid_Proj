#### DocString ####
# This code contains the data simulation and model fitting for the Nowak and May (2000) model
# With the Type III

# AUTHOR: Madeline Jarvis-Cross 
# DATE OF CREATION: 2022-09-07
#### END ####

#### Packages ####
library(rstan)
library(deSolve)
library(GillespieSSA)
library(gdata)
library(bayesplot)
library(tidyverse)
library(brms)
library(beepr)
library(here)

#### Deterministic simulation ####
NowakMay_T3 <- function(t,y,p){
  r <- p[1]; O <- p[2]; h <- p[3]
  b <- p[4]; c <- p[5]; u <- p[6]
  P <- y[1]; H <- y[2]
  dP = P*r - H*(O*(P^2)/(1 + O*(P^2)*h))
  dH = b + H*(c*(O*(P^2)/(1 + O*(P^2)*h)) - u)
  list(c(dP,dH))
}
r <- 3; O <- 0.0015; h <- 0.2
b <- 35; c <- 0.2; u <- 0.5
parms <- c(r, O, h, b, c, u)
P0 <- 1; H0 <- 1
N0 <- c(P0, H0)
TT <- seq(0, 25, 0.1) 
results <- lsoda(N0,TT,NowakMay_T3,parms)
P <- results[,2]; H <- results[,3];
NowakMay_T3_Sim = data.frame(results)
colnames(NowakMay_T3_Sim) <- c("Time", "Parasites", "ImmuneCells")
write.csv(NowakMay_T3_Sim, '02_Data/NowakMay_DetSim_T3.csv')
NowakMay_T3_Plot <- ggplot(NowakMay_T3_Sim,aes(x = Time))+
  geom_line(aes(y = Parasites, color="cornflowerblue"), size = 1)+ # Parasite
  geom_line(aes(y = ImmuneCells,color="springgreen4"), size = 1)+ # Host
  theme_minimal()+
  xlab(" ")+
  ylab(" ")+
  ggtitle("")+
  # ylim(0, 600)+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"),
                     guide = F)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
NowakMay_T3_Plot

#### Stochastic simulation ####
x0 <- c(P = 1, H = 1) 
a <- c("P*r",
       "H*(O*(P^2)/1 + O*(P^2)*h)", 
       "b + H*c*(O*(P^2)/1 + O*(P^2)*h)",
       "H*u")
nu <- matrix(c(+1,-1,0,0,
               0,0,+1,-1), nrow = 2, byrow = TRUE)
r = 3; O = 0.0015; h = 0.2; b = 35; c = 0.2; u = 0.5
parms1 <- c(r = r, O = O, h = h, b = b, c = c, u = u)
tf = 100
method <- "OTL"
simName <- "NowakMay_T3"
set.seed(5)
NowakMay_T3_SLC <- suppressWarnings(ssa(x0, a, nu, parms1, tf, method, simName,
                                            verbose = FALSE, 
                                            consoleInterval = 1, 
                                            censusInterval = 1, 
                                            maxWallTime = 30, 
                                            ignoreNegativeState = TRUE)) 
NowakMay_T3_Sim <- NowakMay_T3_SLC$data 
NowakMay_T3_Sim <- as.data.frame(NowakMay_T3_Sim)
colnames(NowakMay_T3_Sim) <- c("Time", "Parasites", "ImmuneCells")
write.csv(NowakMay_T3_Sim, '02_Data/NowakMay_StochSim_T3.csv')
NowakMay_T3_Plot <- ggplot(NowakMay_T3_Sim,aes(x = Time))+
  geom_line(aes(y = Parasites, color="cornflowerblue"), size = 1)+ 
  geom_line(aes(y = ImmuneCells,color="springgreen4"), size = 1)+
  theme_minimal()+
  xlab(" ")+
  ylab(" ")+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"),
                     guide = F)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
NowakMay_T3_Plot
