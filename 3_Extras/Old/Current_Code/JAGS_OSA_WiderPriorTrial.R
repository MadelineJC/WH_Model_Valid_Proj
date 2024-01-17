# Trying some of the OSA runs with wider priors to see if we still get those weird boundary effects, where the mean estimate of the parm is right-ish, but its CIs are at the bounds/chain isn't concetrated

## Changed the priors to: ##
#r~dunif(0,10) # Intrinsic growth rate of parasite
#O~dunif(0,10) # Rate at which immune system correctly recognizes parasites
#h~dunif(0,10) # Handelling time
#b~dunif(0,1000) # Immigration rate of immune cells to infection site
#c~dunif(0,10) # Growth rate of immune system as a response to infection
#u~dunif(0,10) # Natural mortality rate of immune cells

#### Give NOTHING ####
IC=list(r=2.5, O=0.012, h=0.075, b=35, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
OSA_2_Fit = run.jags(#model = here('./Current_Code/OSA_JAGS_Para.R'),
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Current_Code/OneStepAhead_JAGS_Model_2.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','O','h','b','c','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(3)
summary(OSA_2_Fit) 
plot(OSA_2_Fit)

#### Give r ####
IC=list(O=0.012, h=0.075, b=35, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
OSA_2_Fit_r = run.jags(#model = here('./Current_Code/OSA_JAGS_Para.R'),
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_AddInfoTrials_Models/OneStepAhead_JAGS_Model_2_Give-r.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('O','h','b','c','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(3)
summary(OSA_2_Fit_r) 
plot(OSA_2_Fit_r)

#### Give O ####
IC=list(r=2.5, h=0.075, b=35, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
OSA_2_Fit_O = run.jags(#model = here('./Current_Code/OSA_JAGS_Para.R'),
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_AddInfoTrials_Models/OneStepAhead_JAGS_Model_2_Give-O.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','h','b','c','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(3)
summary(OSA_2_Fit_O) 
plot(OSA_2_Fit_O)

#### Give h ####
IC=list(r=2.5, O=0.012, b=35, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
OSA_2_Fit_h = run.jags(#model = here('./Current_Code/OSA_JAGS_Para.R'),
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_AddInfoTrials_Models/OneStepAhead_JAGS_Model_2_Give-h.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','O','b','c','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(3)
summary(OSA_2_Fit_h) 
plot(OSA_2_Fit_h)

#### Give rO ####
IC=list(h=0.075, b=35, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
OSA_2_Fit_rO = run.jags(#model = here('./Current_Code/OSA_JAGS_Para.R'),
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_AddInfoTrials_Models/OneStepAhead_JAGS_Model_2_Give-rO.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('h','b','c','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(3)
summary(OSA_2_Fit_rO) 
plot(OSA_2_Fit_rO)

#### Give rh ####
IC=list(O=0.012, b=35, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
OSA_2_Fit_rh = run.jags(#model = here('./Current_Code/OSA_JAGS_Para.R'),
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_AddInfoTrials_Models/OneStepAhead_JAGS_Model_2_Give-rh.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('O','b','c','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(3)
summary(OSA_2_Fit_rh) 
plot(OSA_2_Fit_rh)

#### Give rb ####
IC=list(O=0.012, h=0.075, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
OSA_2_Fit_rb = run.jags(#model = here('./Current_Code/OSA_JAGS_Para.R'),
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_AddInfoTrials_Models/OneStepAhead_JAGS_Model_2_Give-rb.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('O','h','c','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(3)
summary(OSA_2_Fit_rb) 
plot(OSA_2_Fit_rb)

#### Give rOh ####
IC=list(b=35, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
OSA_2_Fit_rOh = run.jags(#model = here('./Current_Code/OSA_JAGS_Para.R'),
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_AddInfoTrials_Models/OneStepAhead_JAGS_Model_2_Give-rOh.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('b','c','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(OSA_2_Fit_rOh) 
plot(OSA_2_Fit_rOh)

#### Give rOb ####
IC=list(h=0.075, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
OSA_2_Fit_rOb = run.jags(#model = here('./Current_Code/OSA_JAGS_Para.R'),
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_AddInfoTrials_Models/OneStepAhead_JAGS_Model_2_Give-rOb.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('h','c','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(OSA_2_Fit_rOb) 
plot(OSA_2_Fit_rOb)

#### Give rOc ####
IC=list(h=0.075, b=35, u=0.41); inits=list(IC,IC,IC,IC)
OSA_2_Fit_rOc = run.jags(#model = here('./Current_Code/OSA_JAGS_Para.R'),
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_AddInfoTrials_Models/OneStepAhead_JAGS_Model_2_Give-rOc.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('h','b','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(OSA_2_Fit_rOc) 
plot(OSA_2_Fit_rOc)

#### Give rOhb ####
IC=list(c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
OSA_2_Fit_rOhb = run.jags(#model = here('./Current_Code/OSA_JAGS_Para.R'),
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_AddInfoTrials_Models/OneStepAhead_JAGS_Model_2_Give-rOhb.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('c','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(OSA_2_Fit_rOhb) 
plot(OSA_2_Fit_rOhb)

#### Give rOhc ####
IC=list(b=35, u=0.41); inits=list(IC,IC,IC,IC)
OSA_2_Fit_rOhc = run.jags(#model = here('./Current_Code/OSA_JAGS_Para.R'),
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_AddInfoTrials_Models/OneStepAhead_JAGS_Model_2_Give-rOhc.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('b','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(OSA_2_Fit_rOhc) 
plot(OSA_2_Fit_rOhc)

#### Give rOhu ####
IC=list(b=35, c=0.3); inits=list(IC,IC,IC,IC)
OSA_2_Fit_rOhu = run.jags(#model = here('./Current_Code/OSA_JAGS_Para.R'),
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_AddInfoTrials_Models/OneStepAhead_JAGS_Model_2_Give-rOhu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('b','c'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(OSA_2_Fit_rOhu) 
plot(OSA_2_Fit_rOhu)

#### Give rOhbc ####
IC=list(u=0.41); inits=list(IC,IC,IC,IC)
OSA_2_Fit_rOhbc = run.jags(#model = here('./Current_Code/OSA_JAGS_Para.R'),
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_AddInfoTrials_Models/OneStepAhead_JAGS_Model_2_Give-rOhbc.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(OSA_2_Fit_rOhbc) 
plot(OSA_2_Fit_rOhbc)

#### Outputs ####

OSA_2_Fit_Summ <- data.frame(summary(OSA_2_Fit))
write.csv(OSA_2_Fit_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/OSA_2_Fit.csv",row.names=TRUE)

OSA_2_Fit_r_Summ <- data.frame(summary(OSA_2_Fit_r))
write.csv(OSA_2_Fit_r_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/OSA_2_Fit_r.csv",row.names=TRUE)

OSA_2_Fit_O_Summ <- data.frame(summary(OSA_2_Fit_O))
write.csv(OSA_2_Fit_O_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/OSA_2_Fit_O.csv",row.names=TRUE)

OSA_2_Fit_h_Summ <- data.frame(summary(OSA_2_Fit_h))
write.csv(OSA_2_Fit_h_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/OSA_2_Fit_h.csv",row.names=TRUE)

OSA_2_Fit_rO_Summ <- data.frame(summary(OSA_2_Fit_rO))
write.csv(OSA_2_Fit_rO_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/OSA_2_Fit_rO.csv",row.names=TRUE)

OSA_2_Fit_rh_Summ <- data.frame(summary(OSA_2_Fit_rh))
write.csv(OSA_2_Fit_rh_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/OSA_2_Fit_rh.csv",row.names=TRUE)

OSA_2_Fit_rb_Summ <- data.frame(summary(OSA_2_Fit_rb))
write.csv(OSA_2_Fit_rb_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/OSA_2_Fit_rb.csv",row.names=TRUE)

OSA_2_Fit_rOh_Summ <- data.frame(summary(OSA_2_Fit_rOh))
write.csv(OSA_2_Fit_rOh_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/OSA_2_Fit_rOh.csv",row.names=TRUE)

OSA_2_Fit_rOb_Summ <- data.frame(summary(OSA_2_Fit_rOb))
write.csv(OSA_2_Fit_rOb_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/OSA_2_Fit_rOb.csv",row.names=TRUE)

OSA_2_Fit_rOc_Summ <- data.frame(summary(OSA_2_Fit_rOc))
write.csv(OSA_2_Fit_rOc_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/OSA_2_Fit_rOc.csv",row.names=TRUE)

OSA_2_Fit_rOhb_Summ <- data.frame(summary(OSA_2_Fit_rOhb))
write.csv(OSA_2_Fit_rOhb_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/OSA_2_Fit_rOhb.csv",row.names=TRUE)

OSA_2_Fit_rOhc_Summ <- data.frame(summary(OSA_2_Fit_rOhc))
write.csv(OSA_2_Fit_rOhc_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/OSA_2_Fit_rOhc.csv",row.names=TRUE)

OSA_2_Fit_rOhu_Summ <- data.frame(summary(OSA_2_Fit_rOhu))
write.csv(OSA_2_Fit_rOhu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/OSA_2_Fit_rOhu.csv",row.names=TRUE)

OSA_2_Fit_rOhbc_Summ <- data.frame(summary(OSA_2_Fit_rOhbc))
write.csv(OSA_2_Fit_rOhbc_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/OSA_2_Fit_rOhbc.csv",row.names=TRUE)

#### Plots ####

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/Plots2/OSA_2_Fit.pdf",onefile=TRUE)
plot(OSA_2_Fit)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/Plots2/OSA_2_Fit_r.pdf",onefile=TRUE)
plot(OSA_2_Fit_r)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/Plots2/OSA_2_Fit_O.pdf",onefile=TRUE)
plot(OSA_2_Fit_O)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/Plots2/OSA_2_Fit_h.pdf",onefile=TRUE)
plot(OSA_2_Fit_h)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/Plots2/OSA_2_Fit_rO.pdf",onefile=TRUE)
plot(OSA_2_Fit_rO)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/Plots2/OSA_2_Fit_rh.pdf",onefile=TRUE)
plot(OSA_2_Fit_rh)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/Plots2/OSA_2_Fit_rb.pdf",onefile=TRUE)
plot(OSA_2_Fit_rb)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/Plots2/OSA_2_Fit_rOh.pdf",onefile=TRUE)
plot(OSA_2_Fit_rOh)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/Plots2/OSA_2_Fit_rOb.pdf",onefile=TRUE)
plot(OSA_2_Fit_rOb)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/Plots2/OSA_2_Fit_rOc.pdf",onefile=TRUE)
plot(OSA_2_Fit_rOc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/Plots2/OSA_2_Fit_rOhb.pdf",onefile=TRUE)
plot(OSA_2_Fit_rOhb)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/Plots2/OSA_2_Fit_rOhc.pdf",onefile=TRUE)
plot(OSA_2_Fit_rOhc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/Plots2/OSA_2_Fit_rOhu.pdf",onefile=TRUE)
plot(OSA_2_Fit_rOhu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_3/Output/OSA_JAGS_AddInfoOutputs/Plots2/OSA_2_Fit_rOhbcu.pdf",onefile=TRUE)
plot(OSA_2_Fit_rOhbc)
dev.off()