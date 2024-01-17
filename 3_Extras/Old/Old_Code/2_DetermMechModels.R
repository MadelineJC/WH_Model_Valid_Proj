library(tidyverse)
library(deSolve)

# Modeling the basic fcnal response curves
O <- 0.1 # Where O is theta, denoting the regognition rate of the parasite by the IS
h <- 0.2 # Where h is the handling time the IS requires to subdue a pathogen

TypeIModel <- function(P1){
  FR1 <- O*P1
  return(FR1)
}

TypeIIModel <- function(P2){
  FR2 <- O*P2/(1 + O*h*P2)
  return(FR2)
}

TypeIIIModel <- function(P3){
  FR3 <- O*P3^2/(1 + O*h*P3^2)
  return(FR3)
}

# Making deterministic mechanistic models with pairs of ODEs; analyzing using deSolve
  # Parameters should stay the same through different fcnal responses, but the one's I've set here are garbage


#### Type I ####

TypeIFR_Model <- function(t,y,p){
  # Parameter values
  r <- p[1]; O <- p[2]; # Where r is intrinsic growth rate, and O is theta, which is the recognition rate of the parasite by the IS
  b <- p[3]; c <- p[4]; u <- p[5] # Where b is the immigration rate of immune cells, c is the proliferation rate of the immune system in repsonse to phagocytosis, and u is mu, which is the natural mortality rate of immune cells
  # State variables
  P <- y[1]; H <- y[2]
  # ODE model
  dP = P*r - H*(O*P)
  dH = b + H*(c*(O*P)-u)
  list(c(dP,dH))
}
# Parameters
r <- 3; O <- 0.002;
b <- 25; c <- 0.05; u <- 0.2

# Initial conditions
P0 <- 1; H0 <- 1 # P0 set to 1 to indicate this is beginning of infection
N0 <- c(P0,H0)

# Time points for sol'ns
TT <- seq(0.1,100,0.1)

# Simulate the model
results <- lsoda(N0,TT,TypeIFR_Model,parms)

# Extract results
P <- results[,2]; H <- results[,3];

plot(results)

df = data.frame(results)
TypeIFR_Model <- ggplot(df,aes(x=time))+
  geom_line(aes(y=X1, color="cornflowerblue"))+ # Parasite
  geom_line(aes(y=X2,color="springgreen4"))+ # Host
  labs(title="Within-Host Infection Status Over Time",
       x="Time",
       y="Population Abundance")+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"))
TypeIFR_Model

#### Type II ####

TypeIIFR_Model <- function(t,y,p){
  # Parameter values
  r <- p[1]; O <- p[2]; h <- p[3]
  b <- p[4]; c <- p[5]; u <- p[6]
  # State variables
  P <- y[1]; H <- y[2]
  # ODE model
  dP = P*r - H*(O*P/1 + O*P*h)
  dH = b + H*(c*(O*P/1 + O*P*h)-u)
  list(c(dP,dH))
}
# Parameters
r <- 10; O <- 0.002; h <- 0.3
b <- 2; c <- 0.05; u <- 0.8

# Initial conditions
P0 <- 1; H0 <- 1 # P0 set to 1 to indicate this is beginning of infection
N0 <- c(P0,H0)

# Time points for sol'ns
TT <- seq(0.1,100,0.1)

# Simulate the model
results <- lsoda(N0,TT,TypeIIFR_Model,parms)

# Extract results
P <- results[,2]; H <- results[,3];

plot(results)

df = data.frame(results)
TypeIIFR_Model <- ggplot(df,aes(x=time))+
  geom_line(aes(y=X1, color="cornflowerblue"))+ # Parasite
  geom_line(aes(y=X2,color="springgreen4"))+ # Host
  labs(title="Within-Host Infection Status Over Time",
       x="Time",
       y="Population Abundance")+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"))
TypeIIFR_Model 

#### Type III ####

TypeIIIFR_Model <- function(t,y,p){
  # Parameter values
  r <- p[1]; O <- p[2]; h <- p[3]
  b <- p[4]; c <- p[5]; u <- p[6]
  # State variables
  P <- y[1]; H <- y[2]
  # ODE model
  dP = P*r - H*(O*P^2/1 + O*P*h)
  dH = b + H*(c*(O*P^2/1 + O*P^2*h)-u)
  list(c(dP,dH))
}
# Parameters
r <- 10; O <- 0.002; h <- 0.3
b <- 2; c <- 0.05; u <- 0.8

# Initial conditions
P0 <- 1; H0 <- 1 # P0 set to 1 to indicate this is beginning of infection
N0 <- c(P0,H0)

# Time points for sol'ns
TT <- seq(0.1,500,0.1)

# Simulate the model
results <- lsoda(N0,TT,TypeIIIFR_Model,parms)

# Extract results
P <- results[,2]; H <- results[,3];

plot(results)

df = data.frame(results)
TypeIIIFR_Model <- ggplot(df,aes(x=time))+
  geom_line(aes(y=X1, color="cornflowerblue"))+ # Parasite
  geom_line(aes(y=X2,color="springgreen4"))+ # Host
  labs(title="Within-Host Infection Status Over Time",
       x="Time",
       y="Population Abundance")+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"))
TypeIIIFR_Model


