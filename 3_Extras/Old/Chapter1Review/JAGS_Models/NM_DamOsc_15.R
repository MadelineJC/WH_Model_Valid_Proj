
model { # Keep priors simple for now
  r ~ dunif(0, 10) # Intrinsic growth rate of parasite
  # O ~ dunif(0, 1) # Rate at which immune system correctly recognizes parasites
  # h ~ dunif(0, 1) # Handelling time
  b ~ dunif(0, 1000) # Immigration rate of immune cells to infection site
  c ~ dunif(0, 1) # Growth rate of immune system as a response to infection
  u ~ dunif(0, 1) # Natural mortality rate of immune cells
  
  # State vars
  P[1] <- 80
  H[1] <- 200 # Dat
  
  y_hat_P[1] <- 80
  y_hat_H[1] <- 200 # Est
  
  for (i in 2:L_T){
    P[i] <- P[i-1] + (P[i-1]*r - H[i-1]*(0.008*P[i-1]/(1 + 0.008*P[i-1]*0.06)))*dt
    H[i] <- H[i-1] + (b + H[i-1]*(c*(0.008*P[i-1]/(1 + 0.008*P[i-1]*0.06)) - u))*dt
    
    y_hat_P[i] <- P[i];
    y_hat_H[i] <- H[i]; # Expectation
    
    C_P[i] ~ dpois(y_hat_P[i]);
    C_H[i] ~ dpois(y_hat_H[i]);
  }
}

