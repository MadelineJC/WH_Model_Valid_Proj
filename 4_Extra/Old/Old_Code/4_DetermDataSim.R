#### Let's simulate some data, my dudes ####

library(deSolve)
library(dplyr)
library(ggplot2)

# We want to simulate 8 sets of data for each model
  # (1) Full dataset
  # (2) First quarter of time from onset to equilibrium
  # (3) Second quarter of time from onset to equilibrium
  # (4) Third quarter of time from onset to equilibrium
  # (5) Fourth quarter of time from onset to equilibrium
  # (6) First half of time from onset to equilibrium
  # (7) Second half of time from onset to equilibrium
  # (8) Middle half of time from onset to equilibrium
  # (9) Random subset of data (40% kept)
  # (10) Combed data (every 4th row kept)

#### Type II ####

#### (1) FULL DATASET, FROM ONSET OF INFECTION TO EQUILIBRIUM ####
TypeIIFR_Model <- function(t,y,p){
  # Parameter values
  r <- p[1]; O <- p[2]; h <- p[3]
  b <- p[4]; c <- p[5]; u <- p[6]
  # State variables
  P <- y[1]; H <- y[2]
  # ODE model
  dP = P*r - H*(O*P/(1 + O*P*h))
  dH = b + H*(c*(O*P/(1 + O*P*h)) - u)
  list(c(dP,dH))
}
# Parameters
r <- 2.5; O <- 0.008; h <- 0.06
b <- 35; c <- 0.2; u <- 0.2
parms <- c(r,O,h,b,c,u)

# Initial conditions
P0 <- 10; H0 <- 350 # P0 set to 1 to indicate this is beginning of infection; H0 high because you have an immune system in absence of infection, that grows in response to infection (c)
N0 <- c(P0,H0)

# Time points for sol'ns
TT <- seq(0.1,400,0.1) 

# Simulate the model
results <- lsoda(N0,TT,TypeIIFR_Model,parms)

# Extract results
P <- results[,2]; H <- results[,3];

plot(results) # Quick check to make sure these run to equilibrium, which seems to happen super early anyway (around 120 (originally set to 1000))

# Re-name results and make into data frame
TypeII_Full_DF = data.frame(results)
write.csv(TypeII_Full_DF,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/20Apr_DetermData", row.names = FALSE)

# Plot
TypeII_Full_Plot <- ggplot(TypeII_Full_DF,aes(x=time))+
  geom_line(aes(y=X1, color="cornflowerblue"))+ # Parasite
  geom_line(aes(y=X2,color="springgreen4"))+ # Host
  labs(title="Full Time-Series; Type II",
       x="Time",
       y="Population Abundance")+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"))
TypeII_Full_Plot 

#### SEPARATING DATA INTO SMALLER GROUPS ####

## (2) FIRST QUARTER OF DATA
TypeII_DF_1.30 <- TypeII_Full_DF[-c(301:1200),] # 1-30

TypeII_DF_1.30_Plot <- ggplot(TypeII_DF_1.30,aes(x=time))+
  geom_line(aes(y=X1, color="cornflowerblue"))+ # Parasite
  geom_line(aes(y=X2,color="springgreen4"))+ # Host
  labs(title="First Quarter of TS (1-30)",
       x="Time",
       y="Population Abundance")+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"))
TypeII_DF_1.30_Plot

## (3) SECOND QUARTER OF DATA 
TypeII_DF_30.60 <- TypeII_Full_DF[-c(1:299,601:1200),] # 30-60

TypeII_DF_30.60_Plot <- ggplot(TypeII_DF_30.60,aes(x=time))+
  geom_line(aes(y=X1, color="cornflowerblue"))+ # Parasite
  geom_line(aes(y=X2,color="springgreen4"))+ # Host
  labs(title="Second Quarter of TS (30-60)",
       x="Time",
       y="Population Abundance")+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"))
TypeII_DF_30.60_Plot

## (4) THIRD QUARTER OF DATA 
TypeII_DF_60.90 <- TypeII_Full_DF[-c(1:599,901:1200),] # 60-90

TypeII_DF_60.90_Plot <- ggplot(TypeII_DF_60.90,aes(x=time))+
  geom_line(aes(y=X1, color="cornflowerblue"))+ # Parasite
  geom_line(aes(y=X2,color="springgreen4"))+ # Host
  labs(title="Third Quarter of TS (60-90)",
       x="Time",
       y="Population Abundance")+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"))
TypeII_DF_60.90_Plot

## (5) FOURTH QUARTER OF DATA 
TypeII_DF_90.120 <- TypeII_Full_DF[-c(1:899),] # 90-120

TypeII_DF_90.120_Plot <- ggplot(TypeII_DF_90.120,aes(x=time))+
  geom_line(aes(y=X1, color="cornflowerblue"))+ # Parasite
  geom_line(aes(y=X2,color="springgreen4"))+ # Host
  labs(title="Fourth Quarter of TS (90-120)",
       x="Time",
       y="Population Abundance")+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"))
TypeII_DF_90.120_Plot

## (6) FIRST HALF OF DATA 
TypeII_DF_1.60 <- TypeII_Full_DF[-c(601:1200),] # 1-60

TypeII_DF_1.60_Plot <- ggplot(TypeII_DF_1.60,aes(x=time))+
  geom_line(aes(y=X1, color="cornflowerblue"))+ # Parasite
  geom_line(aes(y=X2,color="springgreen4"))+ # Host
  labs(title="First Half of TS (1-60)",
       x="Time",
       y="Population Abundance")+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"))
TypeII_DF_1.60_Plot

## (7) SECOND HALF OF DATA 
TypeII_DF_60.120 <- TypeII_Full_DF[-c(1:599),] # 60-120

TypeII_DF_60.120_Plot <- ggplot(TypeII_DF_60.120,aes(x=time))+
  geom_line(aes(y=X1, color="cornflowerblue"))+ # Parasite
  geom_line(aes(y=X2,color="springgreen4"))+ # Host
  labs(title="Second Half of TS (60-120)",
       x="Time",
       y="Population Abundance")+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"))
TypeII_DF_60.120_Plot

## (8) MIDDLE HALF OF DATA 
TypeII_DF_30.90 <- TypeII_Full_DF[-c(1:299,901:1200),] # 30-90

TypeII_DF_30.90_Plot <- ggplot(TypeII_DF_30.90,aes(x=time))+
  geom_line(aes(y=X1, color="cornflowerblue"))+ # Parasite
  geom_line(aes(y=X2,color="springgreen4"))+ # Host
  labs(title="Middle Half of TS (30-90)",
       x="Time",
       y="Population Abundance")+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"))
TypeII_DF_30.90_Plot

## (9) RANDOM SUBSET OF DATA (40% KEPT)
nrow(TypeII_Full_DF) # 1200
TypeII_Subset_DF <- TypeII_Full_DF[sample(nrow(TypeII_Full_DF),480),] # Randomly select 40% of data
TypeII_Subset_DF <- TypeII_Subset_DF %>% arrange(time) # Arrange in chronological order

TypeII_Subset_DF_Plot <- ggplot(TypeII_Subset_DF,aes(x=time))+
  geom_line(aes(y=X1, color="cornflowerblue"))+ # Parasite
  geom_line(aes(y=X2,color="springgreen4"))+ # Host
  labs(title="Randomly Sampled Data (40%)",
       x="Time",
       y="Population Abundance")+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"))
TypeII_Subset_DF_Plot

## (10) COMBED DATA (EVERY 4TH ROW KEPT)
toDelete <- seq(1, nrow(TypeII_Full_DF), 2)
TypeII_Combed_DF <- TypeII_Full_DF[toDelete,]
# Do it again
toDelete <- seq(1, nrow(TypeII_Combed_DF), 2)
TypeII_Combed_DF <- TypeII_Combed_DF[toDelete,] # 300 rows

TypeII_Combed_DF_Plot <- ggplot(TypeII_Combed_DF,aes(x=time))+
  geom_line(aes(y=X1, color="cornflowerblue"))+ # Parasite
  geom_line(aes(y=X2,color="springgreen4"))+ # Host
  labs(title="Combed Data (Every 4th Row Present)",
       x="Time",
       y="Population Abundance")+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"))
TypeII_Combed_DF_Plot

#### Type III ####

#### (1) FULL DATASET, FROM ONSET OF INFECTION TO EQUILIBRIUM ####
TypeIIIFR_Model <- function(t,y,p){
  # Parameter values
  r <- p[1]; O <- p[2]; h <- p[3]
  b <- p[4]; c <- p[5]; u <- p[6]
  # State variables
  P <- y[1]; H <- y[2]
  # ODE model
  dP = P*r - H*(O*(P^2)/(1 + O*(P^2)*h))
  dH = b + H*(c*(O*(P^2)/(1 + O*(P^2)*h)) - u)
  list(c(dP,dH))
}
# Parameters
r <- 2.7; O <- 0.0012; h <- 10
b <- 18; c <- 0.02; u <- 0.6

# Initial conditions
P0 <- 1; H0 <- 18 # P0 set to 1 to indicate this is beginning of infection; H0 high because you have an immune system in absence of infection, that grows in response to infection (c)
N0 <- c(P0,H0)

# Time points for sol'ns
TT <- seq(0.1,200,0.1) # THIS SYSTEM TECHNICALLY EQUILIBRIATES AT 2288 (~230 TIME STEPS)!

# Simulate the model
results <- lsoda(N0,TT,TypeIIIFR_Model,parms)

# Extract results
P <- results[,2]; H <- results[,3];

plot(results) # Quick check to make sure these run to equilibrium, which seems to happen super early anyway (around 120 (originally set to 1000))

# Re-name results and make into data frame
TypeIII_Full_DF = data.frame(results)

# Plot
TypeIII_Full_Plot <- ggplot(TypeIII_Full_DF,aes(x=time))+
  geom_line(aes(y=X1, color="cornflowerblue"))+ # Parasite
  geom_line(aes(y=X2,color="springgreen4"))+ # Host
  labs(title="Full Time-Series (1-120)",
       x="Time",
       y="Population Abundance")+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"))
TypeIII_Full_Plot

####

TypeIIIFR_Model <- function(t,y,p){
  # Parameter values
  r <- p[1]; O <- p[2]; h <- p[3]
  b <- p[4]; c <- p[5]; u <- p[6]
  # State variables
  P <- y[1]; H <- y[2]
  # ODE model
  dP = P*r - H*(O*(P^2)/(1 + O*(P^2)*h))
  dH = b + H*(c*(O*(P^2)/(1 + O*(P^2)*h)) - u)
  list(c(dP,dH))
}
# Parameters
r <- 2.7; O <- 0.0012; h <- 10
b <- 18; c <- 0.02; u <- 0.6

# Initial conditions
P0 <- 100; H0 <- 350 # P0 set to 1 to indicate this is beginning of infection; H0 high because you have an immune system in absence of infection, that grows in response to infection (c)
N0 <- c(P0,H0)

# Time points for sol'ns
TT <- seq(0.1,200,0.1) # THIS SYSTEM TECHNICALLY EQUILIBRIATES AT 2288 (~230 TIME STEPS)!

# Simulate the model
results <- lsoda(N0,TT,TypeIIIFR_Model,parms)

# Extract results
P <- results[,2]; H <- results[,3];

plot(results) # Quick check to make sure these run to equilibrium, which seems to happen super early anyway (around 120 (originally set to 1000))

# Re-name results and make into data frame
TypeIII_Full_DF = data.frame(results)

# Plot
TypeIII_Full_Plot <- ggplot(TypeIII_Full_DF,aes(x=time))+
  geom_line(aes(y=X1, color="cornflowerblue"))+ # Parasite
  geom_line(aes(y=X2,color="springgreen4"))+ # Host
  labs(title="Full Time-Series (1-120)",
       x="Time",
       y="Population Abundance")+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"))
TypeIII_Full_Plot
