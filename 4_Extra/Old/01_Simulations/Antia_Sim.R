#### DocString ####
# This code contains the data simulation and model fitting for the Antia et al. model

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

#### Mapping parms ####
# In Nowak and May (2000) model:
# r = Replication rate of parasite
  # Correlates with "r" in Antia model
# O = Recognition rate of parasite by host immune system
  # Correlates with "k" in Antia model
# h = Handling time of parasite by immune system
  # Correlates with "k" in Antia model
# b = Immigration rate of immune cells in absence of infection
  # NO CORRELATE
# c = Activation/proliferation rate of host immune system
  # Correlates with "p" in Antia model
# u = Natural mortality rate of host immune cells
  # NO CORRELATE

#### Deterministic Simulation ####
Antia_Model <- function(t, y, p1){
  r <- p1[1]; k <- p1[2]; p <- p1[3]; o <- p1[4] 
  P <- y[1]; I <- y[2]
  # r = Replication rate of parasite (0.1-10)
  # p = Max. growth rate of immune system (1.0)
  # k = Rate of destruction of parasite by immune system (10^(-3))
  # o = Parasite density at which growth rate of IS is half its max. (10^3)
  dP = r*P - k*P*I
  dI = p*I*(P/(P + o))
  list(c(dP, dI))
}
r <- 0.2; k <- 0.01; p <- 1; o <- 1000 
parms <- c(r, k, p, o)
P0 <- 1; I0 <- 1
N0 <- c(P0, I0)
TT <- seq(0, 50, 0.1)
results <- lsoda(N0, TT, Antia_Model, parms, verbose = TRUE)
P <- results[,2]; I <- results[,3];
plot(results)
Antia_Sim = data.frame(results)
colnames(Antia_Sim) <- c("Time", "Parasites", "ImmuneCells")
write.csv(Antia_Sim, file = here('./02_Data/Antia_DetSim.csv'), row.names=FALSE)
Antia_Plot <- ggplot(Antia_Sim,aes(x = Time))+
  geom_line(aes(y = Parasites, color="cornflowerblue"), size = 1)+ # Parasite
  geom_line(aes(y = ImmuneCells,color="springgreen4"), size = 1)+ # Host
  labs(title=" ",
       x=" ",
       y=" ")+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"),
                     guide = F)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Antia_Plot

#### Stochastic Simulation ####
x0 <- c(P = 1, I = 1) 
a <- c("P*r",
       "k*P*I", 
       "p*I*(P/(P + o))")
nu <- matrix(c(+1, -1, 0,
               0, 0, +1), nrow = 2, byrow = TRUE)

r <- 0.2; k <- 0.01; p <- 1; o <- 1000
parms1 <- c(r = r, k = k, p = p, o = o)
tf = 50
method <- "OTL"
simName <- "Antia"
set.seed(5)
Antia_SLC <- suppressWarnings(ssa(x0, a, nu, parms1, tf, method, simName,
                                  verbose = FALSE, 
                                  consoleInterval = 1, 
                                  censusInterval = 1, 
                                  maxWallTime = 30, 
                                  ignoreNegativeState = TRUE)) 
Antia_Sim <- Antia_SLC$data 
Antia_Sim <- as.data.frame(Antia_Sim)
colnames(Antia_Sim) <- c("Time", "Parasites", "ImmuneCells")

FirstZero <- which(Antia_Sim$Parasites == 0)[1]
Antia_Sim <- Antia_Sim[1:FirstZero - 1, ]

write_csv(Antia_Sim, file = here('./02_Data/Antia_StochSim.csv'))
Antia_Plot <- ggplot(Antia_Sim,aes(x = Time))+
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
Antia_Plot

