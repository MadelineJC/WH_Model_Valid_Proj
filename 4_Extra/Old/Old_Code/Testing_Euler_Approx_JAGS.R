library(rjags)
library(runjags)
library(modeest) 
library(mc2d)
library(readr)
library(here)
library(deSolve)
library(GillespieSSA)
library(tidyverse)

#### Low res data ####
#data = read_csv(here('./Data/NewData.csv'))
data = read_csv("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Data/NewData.csv")

ts = length(data$t)
C_P = as.integer(data$P); C_H = as.integer(data$H)

# Initial parameter values as a list of lists (one list/chain)
IC=list(r=2.5, O=0.012, h=0.075, b=35, c=0.3, u=0.41); inits=list(IC,IC,IC,IC) 

Fit = run.jags(#model = here('./Current_Code/OSA_JAGS_Para.R'),
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/OSA_JAGS.R",
  data = list('ts' = ts, 'C_P' = C_P, 'C_H' = C_H),
  monitor = c('r','O','h','b','c','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 10000,
  adapt = 10000,
  sample = 1000,
  inits = inits,
  summarise = TRUE)

failed.jags(c('model','data','inits')) # To see where stuff goes funky

# Euler approx testing for ts = 1
ts[21]
P[1] <- 80
H[1] <- 200
for (t in 2:ts){
  P[t] <- P[t-1]*r - H[t-1]*(O*P[t-1]/(1 + O*P[t-1]*h))
  H[t] <- b + H[t-1]*(c*(O*P[t-1]/(1 + O*P[t-1]*h)) - u)}

#### High res data ####

# Data gen
SmallTS_D <- function(t,y,p){
  r <- p[1]; O <- p[2]; h <- p[3]
  b <- p[4]; c <- p[5]; u <- p[6]
  P <- y[1]; H <- y[2]
  dP = P*r - H*(O*P/(1 + O*P*h))
  dH = b + H*(c*(O*P/(1 + O*P*h)) - u)
  list(c(dP,dH))
}
r <- 2.5; O <- 0.012; h <- 0.075
b <- 35; c <- 0.3; u <- 0.41
parms <- c(r,O,h,b,c,u)
P0 <- 80; H0 <- 200 
N0 <- c(P0,H0)
TT <- seq(0,17,0.1) 
results <- lsoda(N0,TT,SmallTS_D,parms)
P <- results[,2]; H <- results[,3];
SmallTS_D_DF = data.frame(results)
SmallTS_D_Plot <- ggplot(SmallTS_D_DF,aes(x=time))+
  geom_line(aes(y=X1, color="cornflowerblue"))+ # Parasite
  geom_line(aes(y=X2,color="springgreen4"))+ # Host
  theme_minimal()+
  xlab("Time")+
  ylab("Abundance")+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"))
SmallTS_D_Plot

x0 <- c(P=80,H=200) 
a <- c("P*r",
       "H*(O*P/1 + O*P*h)", 
       "b + H*c*(O*P/1 + O*P*h)",
       "H*u")
nu <- matrix(c(+1,-1,0,0,
               0,0,+1,-1), nrow = 2, byrow = TRUE)
r = 2.5; O = 0.012; h = 0.075; b = 35; c = 0.3; u = 0.41
parms1 <- c(r = r, O = O, h = h, b = b, c = c, u = u)
tf = 20
method <- "OTL"
simName <- "SmallTS"
set.seed(1508)
SmallTS <- suppressWarnings(ssa(x0, a, nu, parms1, tf, method, simName,
                                verbose = FALSE, 
                                consoleInterval = 0.1, 
                                censusInterval = 0.1, 
                                maxWallTime = 30, 
                                ignoreNegativeState = TRUE)) 
SmallTS_Data <- SmallTS$data 
SmallTS_Data <- as.data.frame(SmallTS_Data) # Length = 174 (so like 17 time steps)
SmallTS_Data <- SmallTS_Data[-c(172:174),] # Make 171 to match approx.
SmallTS_Plot <- ggplot(SmallTS_Data,aes(x=t))+
  geom_line(aes(y=P, color="cornflowerblue"))+ 
  geom_line(aes(y=H,color="springgreen4"))+
  theme_minimal()+
  xlab("Time")+
  ylab("Abundance")+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"))
SmallTS_Plot

#### Euler approx; ts = 0.1
r=2.5; O=0.012; h=0.075; b=35; c=0.3; u=0.41
P = vector(mode = 'numeric', length = 171)
H = vector(mode = 'numeric', length = 171)
P[1] <- 80
H[1] <- 200 # Dat
dt <- 0.1 # Eval. interval
F1 <- function(i) {
  F1 = P[i-1]*r - H[i-1]*(O*P[i-1]/(1 + O*P[i-1]*h))
  return(F1)
}
F2 = function(i) {
  F2 <- b + H[i-1]*(c*(O*P[i-1]/(1 + O*P[i-1]*h)) - u)
  return(F2)
}
T_end <- 17 # For 17 time steps, evaluated every 0.1; ASK MARTY ABOUT THIS
TT <- seq(0,T_end,dt) # Time points where evaluation takes place
L_T <- length(TT) # 171
for (i in 2:L_T){
  P[i] <- P[i-1] + F1(i)*dt
  H[i] <- H[i-1] + F2(i)*dt
}

# Run JAGS
data = SmallTS_Data
L_T = length(data$t) #171
C_P = as.integer(data$P); C_H = as.integer(data$H)
dt = 0.1
Fit = run.jags(#model = here('./Current_Code/OSA_JAGS_Para.R'),
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/OSA_JAGS_SmallTS.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','O','h','b','c','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 10000,
  adapt = 10000,
  sample = 1000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE)
summary(Fit) # Shitty Rhats and neffs; estimation not awful?
plot(Fit) # Quite shit

#### Now, let's try to improve the fit, fuck around, etc. ####
# Give model h
Fit2 = run.jags(#model = here('./Current_Code/OSA_JAGS_Para.R'),
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/OSA_JAGS_h.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','b','c','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 10000,
  adapt = 10000,
  sample = 1000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE)
summary(Fit2) # Sucks much worse than same thing in Stan; giving value of h doesn't seem to improve on previous fit
plot(Fit2) # Also about as shit as last fit

# Give model h and O
Fit3 = run.jags(#model = here('./Current_Code/OSA_JAGS_Para.R'),
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/OSA_JAGS_Oh.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','b','c','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 10000,
  adapt = 10000,
  sample = 1000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE)
summary(Fit3) # Much better; Rhats are between 1.0 and 1.1; neffs are between 60 and 600
plot(Fit3) # Chains converged for all parms (u is less good (Rhat = 1.1))

#### Plot Euler approx against deSolve against GSSA data ####
library(tidyverse)
Melt_MegaDF <- reshape2::melt(MegaDF,id.var="TS")
my.labs <- c("Euler Para",
             "Euler Host",
             "Gillespie Para",
             "Gillespie Host",
             "RK Para",
             "RK Host")
MegaPlot <- ggplot(data = Melt_MegaDF, aes(TS,value,colour=variable,linetype=variable)) +
  geom_line()
MegaPlot <- MegaPlot + scale_linetype_manual(name="Type of Data",values=c(1,1,1,1,4,4),
                                             labels=my.labs)
MegaPlot <- MegaPlot + scale_color_manual(name = "Type of Data",
                                          values=c("deepskyblue1","deepskyblue1","plum1","plum1","darkgoldenrod1","darkgoldenrod1"),
                                          labels=my.labs)
MegaPlot <- MegaPlot + theme(
  legend.title=element_blank())+
  theme_minimal()+
  xlab("Time")+
  ylab("Abundance")
MegaPlot

#### ts = 0.01 ####
# Making data
x0 <- c(P=80,H=200) 
a <- c("P*r",
       "H*(O*P/1 + O*P*h)", 
       "b + H*c*(O*P/1 + O*P*h)",
       "H*u")
nu <- matrix(c(+1,-1,0,0,
               0,0,+1,-1), nrow = 2, byrow = TRUE)
r = 2.5; O = 0.012; h = 0.075; b = 35; c = 0.3; u = 0.41
parms1 <- c(r = r, O = O, h = h, b = b, c = c, u = u)
tf = 20
method <- "OTL"
simName <- "SmallTS"
set.seed(1508)
SmallTS <- suppressWarnings(ssa(x0, a, nu, parms1, tf, method, simName,
                                verbose = FALSE, 
                                consoleInterval = 0.01, 
                                censusInterval = 0.01, 
                                maxWallTime = 30, 
                                ignoreNegativeState = TRUE)) 
SmallTS_Data <- SmallTS$data 
SmallTS_Data <- as.data.frame(SmallTS_Data) 
SmallTS_Plot <- ggplot(SmallTS_Data,aes(x=t))+
  geom_line(aes(y=P, color="cornflowerblue"))+ 
  geom_line(aes(y=H,color="springgreen4"))+
  theme_minimal()+
  xlab("Time")+
  ylab("Abundance")+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"))
SmallTS_Plot
data = SmallTS_Data
L_T = length(data$t) #171
C_P = as.integer(data$P); C_H = as.integer(data$H)

dt <- 0.01 # Same for all subsequnt models
Fit4 = run.jags(#model = here('./Current_Code/OSA_JAGS_Para.R'),
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/OSA_JAGS_SmallTS.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','O','h','b','c','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 10000,
  adapt = 10000,
  sample = 1000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE)
summary(Fit4) 
plot(Fit4)

# Give it h
Fit5 = run.jags(#model = here('./Current_Code/OSA_JAGS_Para.R'),
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/OSA_JAGS_h.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','O','b','c','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 10000,
  adapt = 10000,
  sample = 1000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE)
summary(Fit5) 
plot(Fit5)

# Give it O and h
Fit6 = run.jags(#model = here('./Current_Code/OSA_JAGS_Para.R'),
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/OSA_JAGS_Oh.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','b','c','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 10000,
  adapt = 10000,
  sample = 1000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE)
summary(Fit6) 
plot(Fit6) 

#### ts = 0.001 ####
# Making data
x0 <- c(P=80,H=200) 
a <- c("P*r",
       "H*(O*P/1 + O*P*h)", 
       "b + H*c*(O*P/1 + O*P*h)",
       "H*u")
nu <- matrix(c(+1,-1,0,0,
               0,0,+1,-1), nrow = 2, byrow = TRUE)
r = 2.5; O = 0.012; h = 0.075; b = 35; c = 0.3; u = 0.41
parms1 <- c(r = r, O = O, h = h, b = b, c = c, u = u)
tf = 20
method <- "OTL"
simName <- "SmallTS"
set.seed(1508)
SmallTS <- suppressWarnings(ssa(x0, a, nu, parms1, tf, method, simName,
                                verbose = FALSE, 
                                consoleInterval = 0.001, 
                                censusInterval = 0.001, 
                                maxWallTime = 30, 
                                ignoreNegativeState = TRUE)) 
SmallTS_Data <- SmallTS$data 
SmallTS_Data <- as.data.frame(SmallTS_Data) 
SmallTS_Plot <- ggplot(SmallTS_Data,aes(x=t))+
  geom_line(aes(y=P, color="cornflowerblue"))+ 
  geom_line(aes(y=H,color="springgreen4"))+
  theme_minimal()+
  xlab("Time")+
  ylab("Abundance")+
  scale_color_manual(name="Legend",
                     labels=c("Parasite",
                              "Host Immune Cells"),
                     values=c("springgreen4","cornflowerblue"))
SmallTS_Plot
data = SmallTS_Data
L_T = length(data$t) 
C_P = as.integer(data$P); C_H = as.integer(data$H)
dt <- 0.001 # Same for all subsequnt models

Fit7 = run.jags(#model = here('./Current_Code/OSA_JAGS_Para.R'),
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/OSA_JAGS_SmallTS.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','O','h','b','c','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 10000,
  adapt = 10000,
  sample = 1000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE)
summary(Fit7) 
plot(Fit7)

Fit8 = run.jags(#model = here('./Current_Code/OSA_JAGS_Para.R'),
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/OSA_JAGS_h.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','O','b','c','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 10000,
  adapt = 10000,
  sample = 1000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE)
summary(Fit8) 
plot(Fit8)

Fit9 = run.jags(#model = here('./Current_Code/OSA_JAGS_Para.R'),
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/OSA_JAGS_Oh.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','b','c','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 10000,
  adapt = 10000,
  sample = 1000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE)
summary(Fit9) 
plot(Fit9)

#### Summary ####
# Fit gets better as ts gets smaller
# Fit doesn't get noticably better when you give h, but does when you give O and h