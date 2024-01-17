#### Let's add some demographic stochasticity to my simulated data, my dudes ####

library(tidyverse)
library(GillespieSSA)

#### Type II ####

# Numerical vector of initial states
x0 <- c(P=10,H=350) # H0 = b because b is a constant
# Make the individual propensity function
a <- c("P*r",
       "H*(O*P/1 + O*P*h)", 
       "b + H*c*(O*P/1 + O*P*h)",
       "H*u")
# Make state change matrix
nu <- matrix(c(+1,-1,0,0,
               0,0,+1,-1), nrow = 2, byrow = TRUE)
# Set parameters
r = 2.5; O = 0.008; h = 0.06; b = 35; c = 0.2; u = 0.2
parms1 <- c(r = r, O = O, h = h, b = b, c = c, u = u)

# Time frame
tf = 400
# Method
method <- "OTL"
# Sim name
simName <- "WH_IDD_TypeII_FULL"

# Run simulation
set.seed(1905)
SimResults_Stoch_TypeII <- suppressWarnings(ssa(x0, a, nu, parms1, tf, method, simName,
                                                verbose = FALSE, # Indicate status of simulation (slows down sim)
                                                consoleInterval = 0, # Interval at which SSA produces simulation status on console
                                                censusInterval = 0.1, # Interval btw recording the state of the system
                                                maxWallTime = 30, # Cut off simulation at prescribed time
                                                ignoreNegativeState = TRUE)) # Prevent abundances from going negative

Stoch_Data_TypeII <- SimResults_Stoch_TypeII$data # Take just the data part of the simulation and put it in a new variable
Stoch_Data_TypeII <- as.data.frame(Stoch_Data_TypeII)
write.csv(Stoch_Data_TypeII,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/20Apr_StochData.csv", row.names = FALSE)

# Plot
TypeII_Full_Stoch_Plot <- ggplot(Stoch_Data_TypeII,aes(x=t))+
  geom_line(aes(y=P, color="cornflowerblue"))+ # Parasite
  geom_line(aes(y=H,color="springgreen4"))+ # Host
  labs(title="W-H InfecDisDyn w/ Type II FcnalResp and DemStoch",
       x="Time",
       y="Population Abundance")+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"))
TypeII_Full_Stoch_Plot

#### Separating the stochastic simulation into parts ####

# (1) FULL DATASET
Stoch_Data_TypeII
nrow(Stoch_Data_TypeII) # 1953

# (2) FIRST QUARTER 
Stoch_TypeII_DF_1.500 <- Stoch_Data_TypeII[-c(501:1942),] # 1-500

# (3) SECOND QUARTER
Stoch_TypeII_DF_500.1000 <- Stoch_Data_TypeII[-c(1:499,1001:1942),] # 500-1000

# (4) THIRD QUARTER
Stoch_TypeII_DF_1000.1500 <- Stoch_Data_TypeII[-c(1:999,1501:1942),] # 1000-1500

# (5) FOURTH QUARTER
Stoch_TypeII_DF_1500.2000 <- Stoch_Data_TypeII[-c(1:1499),] # 1500-2000

# (6) FIRST HALF
Stoch_TypeII_DF_1.1000 <- Stoch_Data_TypeII[-c(1001:1942),] # 1-1000

# (7) SECOND HALF
Stoch_TypeII_DF_1000.2000 <- Stoch_Data_TypeII[-c(1:999),] # 1000-2000

# (8) MIDDLE HALF
Stoch_TypeII_DF_500.1500 <- Stoch_Data_TypeII[-c(1:499,1501:1942),] # 500-1500

# (9) RANDOM SUBSET OF DATA (20% KEPT)
nrow(Stoch_Data_TypeII) # 1953
Stoch_TypeII_Subset <- Stoch_Data_TypeII[sample(nrow(Stoch_Data_TypeII),(nrow(Stoch_Data_TypeII)*0.2)),] # Randomly select 20% of data
Stoch_TypeII_Subset <- Stoch_TypeII_Subset %>% arrange(t) # Arrange in chronological order

# (10) COMBED DATA (EVERY 4TH ROW KEPT)
toDelete <- seq(1, nrow(Stoch_Data_TypeII), 2)
Stoch_Data_TypeII_Combed <- Stoch_Data_TypeII[toDelete,]
# Do it again
toDelete <- seq(1, nrow(Stoch_Data_TypeII_Combed), 2)
Stoch_Data_TypeII_Combed <- Stoch_Data_TypeII_Combed[toDelete,] # 489 rows

#### Type III ####

x0 <- c(P=1,H=18)
a <- c("P*r",
       "H*(O*(P^2)/1 + O*(P^2)*h)", 
       "b + H*c*(O*(P^2)/1 + O*(P^2)*h)",
       "H*u")
nu <- matrix(c(+1,-1,0,0,
               0,0,+1,-1), nrow = 2, byrow = TRUE)
r = 2.7; O = 0.0012; h = 10; b = 18; c = 0.02; u = 0.6
parms1 <- c(r = r, O = O, h = h, b = b, c = c, u = u)
tf = 200
method <- "OTL"
simName <- "WH_IDD_TypeIII_FULL"

set.seed(2022)
SimResults_Stoch_TypeIII <- suppressWarnings(ssa(x0, a, nu, parms1, tf, method, simName,
                                                 verbose = FALSE,
                                                 consoleInterval = 0,
                                                 censusInterval = 0.1,
                                                 maxWallTime = 30,
                                                 ignoreNegativeState = TRUE))

Stoch_Data_TypeIII <- SimResults_Stoch_TypeIII$data 
Stoch_Data_TypeIII <- as.data.frame(Stoch_Data_TypeIII)

# Plot
TypeIII_Full_Stoch_Plot <- ggplot(Stoch_Data_TypeIII,aes(x=t))+
  geom_line(aes(y=P, color="violetred"))+ # Parasite
  geom_line(aes(y=H,color="darkgoldenrod1"))+ # Host
  labs(title="W-H InfecDisDyn w/ Type III FcnalResp and DemStoch",
       x="Time",
       y="Population Abundance")+
  scale_color_manual(name="Legend",
                     labels=c("Host Immune Cells",
                              "Parasite"),
                     values=c("darkgoldenrod1","violetred"))+
  annotate('text', label = paste('r =',r), x = 20, y = 70)+
  annotate('text', label = paste('O =',O), x = 60, y = 70)+
  annotate('text', label = paste('h =',h), x = 100, y = 70)+
  annotate('text', label = paste('b =',b), x = 140, y = 70)+
  annotate('text', label = paste('c =',c), x = 20, y = 60)+
  annotate('text', label = paste('u =',u), x = 60, y = 60)+
  annotate('text', label = paste('P0 =',"1"), x = 100, y = 60)+
  annotate('text', label = paste('H0 =',"18"), x = 140, y = 60)
TypeIII_Full_Stoch_Plot

#### Separating the stochastic simulation into parts ####

# (1) FULL DATASET
Stoch_Data_TypeII
nrow(Stoch_Data_TypeIII) # 1730

# (2) FIRST QUARTER 
Stoch_TypeIII_DF_1.500 <- Stoch_Data_TypeIII[-c(501:1942),] # 1-500

# (3) SECOND QUARTER
Stoch_TypeIII_DF_500.1000 <- Stoch_Data_TypeIII[-c(1:499,1001:1942),] # 500-1000

# (4) THIRD QUARTER
Stoch_TypeIII_DF_1000.1500 <- Stoch_Data_TypeIII[-c(1:999,1501:1942),] # 1000-1500

# (5) FOURTH QUARTER
Stoch_TypeIII_DF_1500.2000 <- Stoch_Data_TypeIII[-c(1:1499),] # 1500-2000

# (6) FIRST HALF
Stoch_TypeIII_DF_1.1000 <- Stoch_Data_TypeIII[-c(1001:1942),] # 1-1000

# (7) SECOND HALF
Stoch_TypeIII_DF_1000.2000 <- Stoch_Data_TypeIII[-c(1:999),] # 1000-2000

# (8) MIDDLE HALF
Stoch_TypeIII_DF_500.1500 <- Stoch_Data_TypeIII[-c(1:499,1501:1942),] # 500-1500

# (9) RANDOM SUBSET OF DATA (20% KEPT)
nrow(Stoch_Data_TypeIII) # 1730
Stoch_TypeIII_Subset <- Stoch_Data_TypeIII[sample(nrow(Stoch_Data_TypeIII),(nrow(Stoch_Data_TypeIII)*0.2)),] # Randomly select 20% of data
Stoch_TypeIII_Subset <- Stoch_TypeIII_Subset %>% arrange(t) # Arrange in chronological order

# (10) COMBED DATA (EVERY 4TH ROW KEPT)
toDelete <- seq(1, nrow(Stoch_Data_TypeIII), 2)
Stoch_Data_TypeIII_Combed <- Stoch_Data_TypeIII[toDelete,]
# Do it again
toDelete <- seq(1, nrow(Stoch_Data_TypeIII_Combed), 2)
Stoch_Data_TypeIII_Combed <- Stoch_Data_TypeIII_Combed[toDelete,] # 433 rows
