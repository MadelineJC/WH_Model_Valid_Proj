#### INFO ####
# This script implements the Stan models

# Models
library(rstan)
library(deSolve)
library(GillespieSSA)
# Vis
library(gdata)
library(bayesplot)
library(tidyverse)
library(brms)
# Bullshit
library(beepr)

#### Data ####
Data = read.csv("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Data/NewData.csv")
N <- length(Data$t) 
ts <- 1:N
y <- as.matrix(Data[1:N,2:3])
y <- cbind(y[ ,1],y[ ,2]); # Where y[,1] is P, and y[,2] is H
Stan_Data <- list(N=N,ts=ts,y=y)

### Stan script 1 ####
write("
functions{
  real[] dZ_dt(
    real t,
    real[] Z,
    real[] theta,
    real[] x_r,
    int[] x_i){
    real P = Z[1];
    real H = Z[2];

    real r = theta[1];
    real O = theta[2];
    real h = theta[3];
    real b = theta[4];
    real c = theta[5];
    real u = theta[6];

    real dP_dt = P*r - H*(O*P/(1 + O*h*P));
    real dH_dt = b + H*(c*(O*P/(1 + O*h*P))-u);
    return({dP_dt,dH_dt});
  }
}
data {
  int<lower=0>N;
  real ts[N];
  real<lower=0>y[N,2];
}
parameters {
  real<lower=0>r;
  real<lower=0,upper=1>O;
  real<lower=0>h;
  real<lower=0>b;
  real<lower=0>c;
  real<lower=0,upper=1>u;
  real<lower=0>sigma[2];
}
transformed parameters{
real Z[N, 2];
for(i in 2:N){
  Z[i:i,:] = integrate_ode_rk45(dZ_dt, // Function
  y[i-1], // Initial value (empirical data point at previous time step)
  ts[i-1], // Initial time step
  ts[i:i], // Next time step (time step to be solved/estimated)
  {r, O, h, b, c, u},
  rep_array(0.0,2),rep_array(0,2),1e-10,1e-10,2e4);
  }
}
model {
  r~normal(2.5,1);
  O~normal(0.01,2);
  h~normal(0.07,2);
  b~normal(35,1);
  c~normal(0.3,1);
  u~normal(0.4,1);
  sigma~lognormal(-1, 1);
  for (k in 1:2) {
    y[ , k] ~ lognormal(log(Z[ , k]), sigma[k]);
  }
}
generated quantities {
  real y_rep[N, 2];
  for (k in 1:2) {
    for (n in 1:N)
      y_rep[n, k] = lognormal_rng(log(Z[n, k]), sigma[k]);
  }
}",
"Stan_OSA_Model_Working.stan")
Stan_OSA <- stan_model("Stan_OSA_Model_Working.stan")

#init_theta <- list(list(r=2.5,O=0.012,h=0.075,b=35,c=0.3,u=0.41))

Fit <- sampling(Stan_OSA,
                data = Stan_Data,
                warmup = 250,
                iter = 500,
                chains = 2,
                cores = 2,
                thin = 1,
                control = list(max_treedepth = 15,
                               adapt_delta = 0.99),
                #init = init_theta,
                seed = 123,
                check_data = TRUE,
                show_messages = TRUE,
                verbose = TRUE)

#### When next step is based off data ####

n0 = matrix(nrow = 21, ncol = 2)
colnames(n0) = c('P0', 'H0')

tt = matrix(nrow = 21, ncol = 2)
colnames(tt) = c('from', 'to')

results = matrix(nrow = 21, ncol = 3)
colnames(results) = c('time','Para', 'Host')

NewData <- NewData[-c(1)]
NewData <- cbind(c(1:21),NewData)
Para = c(NewData$P)
Host = c(NewData$H)


for(i in 1:21) {
  n0_curr = c(Para[i], Host[i])
  n0[i,] = n0_curr
  tt_curr = seq(i-1,i,1)
  tt[i,] = tt_curr
  results_curr = rk(n0_curr,tt_curr,July1,parms)
  results[i,] = results_curr[2,]
}
OffDataResults <- data.frame(results)
colnames(OffDataResults) <- c('Time','DataP','DataH')

#### When next step is based off previous calculation ####

n02 = matrix(nrow = 21, ncol = 2)
colnames(n02) = c('P0', 'H0')

tt = matrix(nrow = 21, ncol = 2)
colnames(tt) = c('from', 'to')

results2 = matrix(nrow = 21, ncol = 3)
colnames(results2) = c('time','Para', 'Host')

NewData <- NewData[-c(1)]
NewData <- cbind(c(1:21),NewData)
Para = c(NewData$P)
Host = c(NewData$H)

n02_curr = c(Para[1], Host[1])
n02[1,] = n02_curr
tt_curr = seq(0,1,1)
tt[1,] = tt_curr
results2_curr = rk(n02_curr,tt_curr,July1,parms)
results2[1,] = results2_curr[2,]
j = 1

for(i in 2:21) {
  n03_curr = c(as.numeric(results2[j,2]), as.numeric(results2[j,3]))
  n02[i,] = n03_curr
  tt_curr = seq(i-1,i,1)
  tt[i,] = tt_curr
  results2_curr = rk(n03_curr,tt_curr,July1,parms)
  results2[i,] = results2_curr[2,]
  j = j+1
}
OffCalcResults <- data.frame(results2)
colnames(OffCalcResults) <- c('Time','CalcP','CalcH')

#### DF garbage
## 4 dfs
# Original:
July1_DF # X1, X2
# Stoch:
NewData # P, H
NewData <- NewData[-c(1)]
NewData <- NewData[-c(21),]
# OffData:
OffDataResults # DataP, DataH
OffDataResults <- OffDataResults[-c(1)]
OffDataResults <- rbind(c(80,200),OffDataResults)
OffDataResults <- OffDataResults[-c(21:22),]
# OffCalc:
OffCalcResults # CalcP, CalcH
OffCalcResults <- OffCalcResults[-c(1)]
OffCalcResults <- rbind(c(80,200),OffCalcResults)
OffCalcResults <- OffCalcResults[-c(21:22),]
# View
View(July1_DF)
View(NewData)
View(OffDataResults)
View(OffCalcResults)
# Combine
DF <- cbind(NewData,July1_DF)
colnames(DF) <- c('NewDataP','NewDataH','Time','July1P','July1H')
DF2 <- cbind(DF,OffDataResults)
DF3 <- cbind(DF2,OffCalcResults)
DF <- DF3

#### Plot
Plot <- ggplot(DF,aes(x=Time))+
  geom_line(aes(y=NewDataP),linetype="solid",colour="black")+ # Stoch para
  geom_line(aes(y=NewDataH),linetype="solid",colour="black")+ # Stoch host
  geom_line(aes(y=July1P),linetype="twodash",colour="red")+ # Det para
  geom_line(aes(y=July1H),linetype="twodash",colour="red")+ # Det host
  geom_line(aes(y=DataP),linetype="dotdash",colour="blue")+ # Off data para
  geom_line(aes(y=DataH),linetype="dotdash",colour="blue")+ # Off data host
  geom_line(aes(y=CalcP),linetype="dotted",colour="green")+ # Off calc para
  geom_line(aes(y=CalcH),linetype="dotted",colour="green")+ # Off calc host
  scale_x_continuous(breaks=seq(0,20,by=1))+ # Cleaning up the axes
  scale_y_continuous(breaks=seq(0,300,by=50))+
  xlab("Time")+
  ylab("Population Abundance")+
  ggtitle("ODE Solutions using RK")+
  theme_minimal()
Plot

#### Stan script 2 ####
# Based off Conor's script from the Stan forum
write("
functions{
  real[] m_dts(int N, real[] Z, real[] theta){
    real P[N];
    real H[N];
    P[1] = Z[1];
    H[1] = Z[2];

    for (t in 1:(N-1)){
      P[t+1] = P[t]*theta[1] - H[t]*(theta[2]*P[t]/(1 + theta[2]*theta[3]*P[t]));
      H[t+1] = theta[4] + H[t]*(theta[5]*(theta[2]*P[t]/(1 + theta[2]*theta[3]*P[t]))-theta[6]);
    }
    return P[];
    return H[];
  }
}

data {
  int<lower=0>N; 
  real<lower=0>y[N,2]; 
}

transformed data{
}

parameters {
  real<lower=0>r;
  real<lower=0,upper=1>O;
  real<lower=0>h;
  real<lower=0>b;
  real<lower=0>c;
  real<lower=0,upper=1>u;
  real<lower=0>y_init_1[1];
  real<lower=0>y_init_2[1];
  real<lower=0>sigma[2];
}

transformed parameters{
  real Z[N, 2];
  real theta[6];

  theta[1] = r;
  theta[2] = O;
  theta[3] = h;
  theta[4] = b;
  theta[5] = c;
  theta[6] = u;
  
  Z[1] = m_dts(N,y_init_1,theta);
  Z[2] = m_dts(N,y_init_2,theta);
}

model {
  r~normal(2.5,1);
  O~normal(0.01,2);
  h~normal(0.07,2);
  b~normal(35,1);
  c~normal(0.3,1);
  u~normal(0.4,1);
  y_init_1~lognormal(log(140), 1);
  y_init_2~lognormal(log(140), 1);
  sigma~lognormal(-1, 1);
  y[1] ~ lognormal(log(Z[1]), sigma[1]);
  y[2] ~ lognormal(log(Z[2]), sigma[2]);
}

generated quantities {
  real y_rep[N, 2];
  y_rep[1] = lognormal_rng(log(Z[1]), sigma[1]);
  y_rep[2] = lognormal_rng(log(Z[2]), sigma[2]);
}",
"Stan_OSA_Model_Working2.stan")
Stan_OSA_2 <- stan_model("Stan_OSA_Model_Working2.stan")

Fit <- sampling(Stan_OSA_2,
                data = Stan_Data,
                warmup = 250,
                iter = 500,
                chains = 2,
                cores = 2,
                thin = 1,
                control = list(max_treedepth = 15,
                               adapt_delta = 0.99),
                #init = init_theta,
                seed = 123,
                check_data = TRUE,
                show_messages = TRUE,
                verbose = TRUE)