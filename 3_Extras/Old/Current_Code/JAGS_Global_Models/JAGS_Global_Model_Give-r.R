# Functional response project
## One-step-ahead model with JAGS instead of Stan because I'm stupid

# data {
#   for (i in 2:ts) {
#     New[i] <- y[i] - y[i-1]
#   }
# }

model { # Keep priors simple for now
  #r~dunif(0,5) # Intrinsic growth rate of parasite
  O~dunif(0,1) # Rate at which immune system correctly recognizes parasites
  h~dunif(0,1) # Handelling time
  b~dunif(0,100) # Immigration rate of immune cells to infection site
  c~dunif(0,1) # Growth rate of immune system as a response to infection
  u~dunif(0,1) # Natural mortality rate of immune cells
  
  #r~dunif(2,3) # Intrinsic growth rate of parasite
  #O~dunif(0,0.5) # Rate at which immune system correctly recognizes parasites
  #h~dunif(0,0.5) # Handelling time
  #b~dunif(30,40) # Immigration rate of immune cells to infection site
  #c~dunif(0,1) # Growth rate of immune system as a response to infection
  #u~dunif(0,1) # Natural mortality rate of immune cells
  
  # State vars
  P[1] <- 80
  H[1] <- 200 # Dat
  
  y_hat_P[1] <- 80
  y_hat_H[1] <- 200 # Est
  
  for (i in 2:L_T){
    P[i] <- P[i-1] + (P[i-1]*2.5 - H[i-1]*(O*P[i-1]/(1 + O*P[i-1]*h)))*dt
    H[i] <- H[i-1] + (b + H[i-1]*(c*(O*P[i-1]/(1 + O*P[i-1]*h)) - u))*dt
    
    y_hat_P[i] <- P[i];
    y_hat_H[i] <- H[i]; # Expectation
    
    C_P[i] ~ dpois(y_hat_P[i]);
    C_H[i] ~ dpois(y_hat_H[i]);
  }
}
