#### DocString ####
# This code contains the data simulation and model fitting for the Nowak and May (2000) model
# With the Type II, unstable limit cycle parameterization

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
NowakMay_T2_Div <- function(t,y,p){
  r <- p[1]; O <- p[2]; h <- p[3]
  c <- p[4]; u <- p[5]
  P <- y[1]; H <- y[2]
  dP = P*r - H*(O*P/(1 + O*P*h))
  dH = H*(c*(O*P/(1 + O*P*h)) - u)
  list(c(dP,dH))
}
r <- 1; O <- 0.5; h <- 0.1
c <- 0.8; u <- 0.2
parms <- c(r, O, h, c, u)
P0 <- 1; H0 <- 1 
N0 <- c(P0, H0)
TT <- seq(0, 300, 0.1) 
results <- lsoda(N0,TT,NowakMay_T2_Div,parms)
P <- results[,2]; H <- results[,3];
NowakMay_T2_Div_Sim = data.frame(results)
colnames(NowakMay_T2_Div_Sim) <- c("Time", "Parasites", "ImmuneCells")
write_csv(NowakMay_T2_Div_Sim, here('./02_Data/NowakMay_DetSim_T2-Div.csv'))
NowakMay_T2_Div_Plot <- ggplot(NowakMay_T2_Div_Sim,aes(x = Time))+
  geom_line(aes(y = Parasites, color="cornflowerblue"), size = 1)+ # Parasite
  geom_line(aes(y = ImmuneCells,color="springgreen4"), size = 1)+ # Host
  theme_minimal()+
  xlab(" ")+
  ylab(" ")+
  ggtitle("")+
  #ylim(0, 300)+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"),
                     guide = F)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
NowakMay_T2_Div_Plot

#### Stochastic simulation ####
x0 <- c(P = 1, H = 1) 
a <- c("P*r",
       "H*(O*P/1 + O*P*h)", 
       "H*c*(O*P/1 + O*P*h)",
       "H*u")
nu <- matrix(c(+1,-1,0,0,
               0,0,+1,-1), nrow = 2, byrow = TRUE)
r = 1; O = 0.5; h = 0.1; c = 0.8; u = 0.2
parms1 <- c(r = r, O = O, h = h, c = c, u = u)
tf = 100
method <- "OTL"
simName <- "NowakMay_T2_Div"
set.seed(9)
NowakMay_T2_Div_SLC <- suppressWarnings(ssa(x0, a, nu, parms1, tf, method, simName,
                                            verbose = FALSE, 
                                            consoleInterval = 1, 
                                            censusInterval = 1, 
                                            maxWallTime = 30, 
                                            ignoreNegativeState = TRUE)) 
NowakMay_T2_Div_Sim <- NowakMay_T2_Div_SLC$data 
NowakMay_T2_Div_Sim <- as.data.frame(NowakMay_T2_Div_Sim)
colnames(NowakMay_T2_Div_Sim) <- c("Time", "Parasites", "ImmuneCells")
write_csv(NowakMay_T2_Div_Sim, here('./02_Data/NowakMay_StochSim_T2-Div.csv'))
NowakMay_T2_Div_Plot <- ggplot(NowakMay_T2_Div_Sim,aes(x = Time))+
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
NowakMay_T2_Div_Plot