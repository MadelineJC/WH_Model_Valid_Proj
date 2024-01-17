#### Description ####
# This script initializes model fitting for 63 models
# Each of the 63 models differ in which parameters (out of 6 total) are being estimated
# The model fitting is being implemented using a global fitting framework in JAGS (Just Another Gibbs Sampler)

library(rjags)
library(runjags)
library(modeest) 
library(mc2d)
library(readr)
library(beepr)

#### Things Needed ####

#IC=list(r=2.5, O=0.012, h=0.075, b=35, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
#monitor = c('r','O','h','b','c','u')

SmallTS_Data <- read_csv("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Data/SmallTS_Data.csv")
View(SmallTS_Data)
data = SmallTS_Data
L_T = length(data$t)
C_P = as.integer(data$P); C_H = as.integer(data$H)
dt <- 0.001 # Same for all subsequnt models; small step size for "better" Euler approx.

#### Give NOTHING ####
IC=list(r=2.5, O=0.012, h=0.075, b=35, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model.R",
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
summary(Global_Fit) 
plot(Global_Fit)

#### Give 1 parm ####

#### Give r ####
IC=list(O=0.012, h=0.075, b=35, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_r = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-r.R",
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
summary(Global_Fit_r) 
plot(Global_Fit_r)

#### Give O ####
IC=list(r=2.5, h=0.075, b=35, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_O = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-O.R",
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
summary(Global_Fit_O) 
plot(Global_Fit_O)

#### Give h ####
IC=list(r=2.5, O=0.012, b=35, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_h = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-h.R",
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
summary(Global_Fit_h) 
plot(Global_Fit_h)

#### Give b ####
IC=list(r=2.5, O=0.012, h=0.075, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_b = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-b.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','O','h','c','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(3)
summary(Global_Fit_b) 
plot(Global_Fit_b)

#### Give c ####
IC=list(r=2.5, O=0.012, h=0.075, b=35, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_c = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-c.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','O','h','b','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(3)
summary(Global_Fit_c) 
plot(Global_Fit_c)

#### Give u ####
IC=list(r=2.5, O=0.012, h=0.075, b=35, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_u = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-u.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','O','h','b','c'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(3)
summary(Global_Fit_u) 
plot(Global_Fit_u)

#### Give 2 parms ####

#### Give rO ####
IC=list(h=0.075, b=35, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rO = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rO.R",
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
summary(Global_Fit_rO) 
plot(Global_Fit_rO)

#### Give rh ####
IC=list(O=0.012, b=35, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rh = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rh.R",
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
summary(Global_Fit_rh) 
plot(Global_Fit_rh)

#### Give rb ####
IC=list(O=0.012, h=0.075, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rb = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rb.R",
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
summary(Global_Fit_rb) 
plot(Global_Fit_rb)

#### Give rc ####
IC=list(O=0.012, h=0.075, b=35, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rc = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rc.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('O','h','b','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(3)
summary(Global_Fit_rc) 
plot(Global_Fit_rc)

#### Give ru ####
IC=list(O=0.012, h=0.075, b=35, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_ru = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-ru.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('O','h','b','c'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(3)
summary(Global_Fit_ru) 
plot(Global_Fit_ru)

#### Give Oh ####
IC=list(r=2.5, b=35, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_Oh = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-Oh.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','b','c','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(3)
summary(Global_Fit_Oh) 
plot(Global_Fit_Oh)

#### Give Ob ####
IC=list(r=2.5, h=0.075, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_Ob = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-Ob.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','h','c','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(3)
summary(Global_Fit_Ob) 
plot(Global_Fit_Ob)

#### Give Oc ####
IC=list(r=2.5, h=0.075, b=35, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_Oc = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-Oc.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','h','b','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(3)
summary(Global_Fit_Oc) 
plot(Global_Fit_Oc)

#### Give Ou ####
IC=list(r=2.5, h=0.075, b=35, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_Ou = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-Ou.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','h','b','c'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(3)
summary(Global_Fit_Ou) 
plot(Global_Fit_Ou)

#### Give hb ####
IC=list(r=2.5, O=0.012, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_hb = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-hb.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','O','c','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(3)
summary(Global_Fit_hb) 
plot(Global_Fit_hb)

#### Give hc ####
IC=list(r=2.5, O=0.012, b=35, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_hc = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-hc.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','O','b','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(3)
summary(Global_Fit_hc) 
plot(Global_Fit_hc)

#### Give hu ####
IC=list(r=2.5, O=0.012, b=35, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_hu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-hu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','O','b','c'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(3)
summary(Global_Fit_hu) 
plot(Global_Fit_hu)

#### Give bc ####
IC=list(r=2.5, O=0.012, h=0.075, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_bc = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-bc.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','O','h','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(3)
summary(Global_Fit_bc) 
plot(Global_Fit_bc)

#### Give bu ####
IC=list(r=2.5, O=0.012, h=0.075, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_bu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-bu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','O','h','c'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(3)
summary(Global_Fit_bu) 
plot(Global_Fit_bu)

#### Give cu ####
IC=list(r=2.5, O=0.012, h=0.075, b=35); inits=list(IC,IC,IC,IC)
Global_Fit_cu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-cu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','O','h','b'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(3)
summary(Global_Fit_cu) 
plot(Global_Fit_cu)

#### Give 3 parms ####

#### Give rOh ####
IC=list(b=35, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rOh = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rOh.R",
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
summary(Global_Fit_rOh) 
plot(Global_Fit_rOh)

#### Give rOb ####
IC=list(h=0.075, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rOb = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rOb.R",
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
summary(Global_Fit_rOb) 
plot(Global_Fit_rOb)

#### Give rOc ####
IC=list(h=0.075, b=35, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rOc = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rOc.R",
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
summary(Global_Fit_rOc) 
plot(Global_Fit_rOc)

#### Give rOu ####
IC=list(h=0.075, b=35, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_rOu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rOu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('h','b','c'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_rOu) 
plot(Global_Fit_rOu)

#### Give rhb ####
IC=list(O=0.012, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rhb = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rhb.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('O','c','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_rhb) 
plot(Global_Fit_rhb)

#### Give rhc ####
IC=list(O=0.012, b=35, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rhc = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rhc.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('O','b','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_rhc) 
plot(Global_Fit_rhc)

#### Give rhu ####
IC=list(O=0.012, b=35, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_rhu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rhu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('O','b','c'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_rhu) 
plot(Global_Fit_rhu)

#### Give rbc ####
IC=list(O=0.012, h=0.075, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rbc = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rbc.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('O','h','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_rbc) 
plot(Global_Fit_rbc)

#### Give rbu ####
IC=list(O=0.012, h=0.075, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_rbu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rbu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('O','h','c'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_rbu) 
plot(Global_Fit_rbu)

#### Give rcu ####
IC=list(O=0.012, h=0.075, b=35); inits=list(IC,IC,IC,IC)
Global_Fit_rcu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rcu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('O','h','b'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_rcu) 
plot(Global_Fit_rcu)

#### Give Ohb ####
IC=list(r=2.5, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_Ohb = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-Ohb.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','c','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_Ohb) 
plot(Global_Fit_Ohb)

#### Give Ohc ####
IC=list(r=2.5, b=35, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_Ohc = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-Ohc.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','b','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_Ohc) 
plot(Global_Fit_Ohc)

#### Give Ohu ####
IC=list(r=2.5, b=35, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_Ohu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-Ohu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','b','c'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_Ohu) 
plot(Global_Fit_Ohu)

#### Give Obc ####
IC=list(r=2.5, h=0.075, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_Obc = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-Obc.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','h','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_Obc) 
plot(Global_Fit_Obc)

#### Give Obu ####
IC=list(r=2.5, h=0.075, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_Obu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-Obu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','h','c'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_Obu) 
plot(Global_Fit_Obu)

#### Give Ocu ####
IC=list(r=2.5, h=0.075, b=35); inits=list(IC,IC,IC,IC)
Global_Fit_Ocu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-Ocu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','h','b'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_Ocu) 
plot(Global_Fit_Ocu)

#### Give hbc ####
IC=list(r=2.5, O=0.012, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_hbc = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-hbc.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','O','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_hbc) 
plot(Global_Fit_hbc)

#### Give hbu ####
IC=list(r=2.5, O=0.012, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_hbu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-hbu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','O','c'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_hbu) 
plot(Global_Fit_hbu)

#### Give hcu ####
IC=list(r=2.5, O=0.012, b=35); inits=list(IC,IC,IC,IC)
Global_Fit_hcu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-hcu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','O','b'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_hcu) 
plot(Global_Fit_hcu)

#### Give bcu ####
IC=list(r=2.5, O=0.012, h=0.075); inits=list(IC,IC,IC,IC)
Global_Fit_bcu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-bcu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','O','h'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_bcu) 
plot(Global_Fit_bcu)

#### Give 4 parms ####

#### Give rOhb ####
IC=list(c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rOhb = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rOhb.R",
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
summary(Global_Fit_rOhb) 
plot(Global_Fit_rOhb)

#### Give rOhc ####
IC=list(b=35, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rOhc = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rOhc.R",
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
summary(Global_Fit_rOhc) 
plot(Global_Fit_rOhc)

#### Give rOhu ####
IC=list(b=35, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_rOhu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rOhu.R",
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
summary(Global_Fit_rOhu) 
plot(Global_Fit_rOhu)

#### Give rObc ####
IC=list(h=0.075, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rObc = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rObc.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('h','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_rObc) 
plot(Global_Fit_rObc)

#### Give rObu ####
IC=list(h=0.075, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_rObu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rObu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('h','c'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_rObu) 
plot(Global_Fit_rObu)

#### Give rOcu ####
IC=list(h=0.075, b=35); inits=list(IC,IC,IC,IC)
Global_Fit_rOcu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rOcu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('h','b'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_rOcu) 
plot(Global_Fit_rOcu)

#### Give rhbc ####
IC=list(O=0.012, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rhbc = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rhbc.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('O','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_rhbc) 
plot(Global_Fit_rhbc)

#### Give rhbu ####
IC=list(O=0.012, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_rhbu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rhbu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('O','c'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_rhbu) 
plot(Global_Fit_rhbu)

#### Give rhcu ####
IC=list(O=0.012, b=35); inits=list(IC,IC,IC,IC)
Global_Fit_rhcu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rhcu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('O','b'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_rhcu) 
plot(Global_Fit_rhcu)

#### Give rbcu ####
IC=list(O=0.012, h=0.075); inits=list(IC,IC,IC,IC)
Global_Fit_rbcu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rbcu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('O','h'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_rbcu) 
plot(Global_Fit_rbcu)

#### Give Ohbc ####
IC=list(r=2.5, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_Ohbc = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-Ohbc.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','u'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_Ohbc) 
plot(Global_Fit_Ohbc)

#### Give Ohbu ####
IC=list(r=2.5, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_Ohbu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-Ohbu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','c'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_Ohbu) 
plot(Global_Fit_Ohbu)

#### Give Ohcu ####
IC=list(r=2.5, b=35); inits=list(IC,IC,IC,IC)
Global_Fit_Ohcu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-Ohcu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','b'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_Ohcu) 
plot(Global_Fit_Ohcu)

#### Give Obcu ####
IC=list(r=2.5, h=0.075); inits=list(IC,IC,IC,IC)
Global_Fit_Obcu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-Obcu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','h'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_Obcu) 
plot(Global_Fit_Obcu)

#### Give hbcu ####
IC=list(r=2.5, O=0.012); inits=list(IC,IC,IC,IC)
Global_Fit_hbcu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-hbcu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r','O'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_hbcu) 
plot(Global_Fit_hbcu)

#### Give 5 parms ####

#### Give rOhbc ####
IC=list(u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rOhbc = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rOhbc.R",
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
summary(Global_Fit_rOhbc) 
plot(Global_Fit_rOhbc)

#### Give rOhbu ####
IC=list(c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_rOhbu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rOhbu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('c'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_rOhbu) 
plot(Global_Fit_rOhbu)

#### Give rOhcu ####
IC=list(b=35); inits=list(IC,IC,IC,IC)
Global_Fit_rOhcu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rOhcu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('b'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_rOhcu) 
plot(Global_Fit_rOhcu)

#### Give rObcu ####
IC=list(h=0.075); inits=list(IC,IC,IC,IC)
Global_Fit_rObcu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rObcu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('h'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_rObcu) 
plot(Global_Fit_rObcu)

#### Give rhbcu ####
IC=list(O=0.012); inits=list(IC,IC,IC,IC)
Global_Fit_rhbcu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-rhbcu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('O'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_rhbcu) 
plot(Global_Fit_rhbcu)

#### Give Ohbcu ####
IC=list(r=2.5); inits=list(IC,IC,IC,IC)
Global_Fit_Ohbcu = run.jags(
  model = "/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/JAGS_Global_Models/JAGS_Global_Model_Give-Ohbcu.R",
  data = list('L_T' = L_T, 'C_P' = C_P, 'C_H' = C_H, 'dt' = dt),
  monitor = c('r'),
  n.chains = 4,
  method = 'rjags',
  burnin = 20000,
  adapt = 20000,
  sample = 2000,
  inits = inits,
  summarise = TRUE,
  plots = TRUE); beep(2)
summary(Global_Fit_Ohbcu) 
plot(Global_Fit_Ohbcu)

#### Make outputs in .csv ####

Fit_Summ <- data.frame(summary(Fit))
write.csv(Fit_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit.csv",row.names=TRUE)

Fit_r_Summ <- data.frame(summary(Fit_r))
write.csv(Fit_r_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_r.csv",row.names=TRUE)

Fit_O_Summ <- data.frame(summary(Fit_O))
write.csv(Fit_O_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_O.csv",row.names=TRUE)

Fit_h_Summ <- data.frame(summary(Fit_h))
write.csv(Fit_h_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_h.csv",row.names=TRUE)

Fit_b_Summ <- data.frame(summary(Fit_b))
write.csv(Fit_b_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_b.csv",row.names=TRUE)

Fit_c_Summ <- data.frame(summary(Fit_c))
write.csv(Fit_c_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_c.csv",row.names=TRUE)

Fit_u_Summ <- data.frame(summary(Fit_u))
write.csv(Fit_u_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_u.csv",row.names=TRUE)

Fit_rO_Summ <- data.frame(summary(Fit_rO))
write.csv(Fit_rO_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rO.csv",row.names=TRUE)

Fit_rh_Summ <- data.frame(summary(Fit_rh))
write.csv(Fit_rh_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rh.csv",row.names=TRUE)

Fit_rb_Summ <- data.frame(summary(Fit_rb))
write.csv(Fit_rb_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rb.csv",row.names=TRUE)

Fit_rc_Summ <- data.frame(summary(Fit_rc))
write.csv(Fit_rc_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rc.csv",row.names=TRUE)

Fit_ru_Summ <- data.frame(summary(Fit_ru))
write.csv(Fit_ru_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_ru.csv",row.names=TRUE)

Fit_Oh_Summ <- data.frame(summary(Fit_Oh))
write.csv(Fit_Oh_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_Oh.csv",row.names=TRUE)

Fit_Ob_Summ <- data.frame(summary(Fit_Ob))
write.csv(Fit_Ob_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_Ob.csv",row.names=TRUE)

Fit_Oc_Summ <- data.frame(summary(Fit_Oc))
write.csv(Fit_Oc_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_Oc.csv",row.names=TRUE)

Fit_Ou_Summ <- data.frame(summary(Fit_Ou))
write.csv(Fit_Ou_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_Ou.csv",row.names=TRUE)

Fit_hb_Summ <- data.frame(summary(Fit_hb))
write.csv(Fit_hb_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_hb.csv",row.names=TRUE)

Fit_hc_Summ <- data.frame(summary(Fit_hc))
write.csv(Fit_hc_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_hc.csv",row.names=TRUE)

Fit_hu_Summ <- data.frame(summary(Fit_hu))
write.csv(Fit_hu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_hu.csv",row.names=TRUE)

Fit_bc_Summ <- data.frame(summary(Fit_bc))
write.csv(Fit_bc_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_bc.csv",row.names=TRUE)

Fit_bu_Summ <- data.frame(summary(Fit_bu))
write.csv(Fit_bu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_bu.csv",row.names=TRUE)

Fit_cu_Summ <- data.frame(summary(Fit_cu))
write.csv(Fit_cu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_cu.csv",row.names=TRUE)

Fit_rOh_Summ <- data.frame(summary(Fit_rOh))
write.csv(Fit_rOh_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rOh.csv",row.names=TRUE)

Fit_rOb_Summ <- data.frame(summary(Fit_rOb))
write.csv(Fit_rOb_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rOb.csv",row.names=TRUE)

Fit_rOc_Summ <- data.frame(summary(Fit_rOc))
write.csv(Fit_rOc_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rOc.csv",row.names=TRUE)

Fit_rOu_Summ <- data.frame(summary(Fit_rOu))
write.csv(Fit_rOu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rOu.csv",row.names=TRUE)

Fit_rhb_Summ <- data.frame(summary(Fit_rhb))
write.csv(Fit_rhb_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rhb.csv",row.names=TRUE)

Fit_rhc_Summ <- data.frame(summary(Fit_rhc))
write.csv(Fit_rhc_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rhc.csv",row.names=TRUE)

Fit_rhu_Summ <- data.frame(summary(Fit_rhu))
write.csv(Fit_rhu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rhu.csv",row.names=TRUE)

Fit_rbc_Summ <- data.frame(summary(Fit_rbc))
write.csv(Fit_rbc_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rbc.csv",row.names=TRUE)

Fit_rbu_Summ <- data.frame(summary(Fit_rbu))
write.csv(Fit_rbu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rbu.csv",row.names=TRUE)

Fit_rcu_Summ <- data.frame(summary(Fit_rcu))
write.csv(Fit_rcu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rcu.csv",row.names=TRUE)

Fit_Ohb_Summ <- data.frame(summary(Fit_Ohb))
write.csv(Fit_Ohb_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_Ohb.csv",row.names=TRUE)

Fit_Ohc_Summ <- data.frame(summary(Fit_Ohc))
write.csv(Fit_Ohc_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_Ohc.csv",row.names=TRUE)

Fit_Ohu_Summ <- data.frame(summary(Fit_Ohu))
write.csv(Fit_Ohu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_Ohu.csv",row.names=TRUE)

Fit_Obc_Summ <- data.frame(summary(Fit_Obc))
write.csv(Fit_Obc_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_Obc.csv",row.names=TRUE)

Fit_Obu_Summ <- data.frame(summary(Fit_Obu))
write.csv(Fit_Obu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_Obu.csv",row.names=TRUE)

Fit_Ocu_Summ <- data.frame(summary(Fit_Ocu))
write.csv(Fit_Ocu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_Ocu.csv",row.names=TRUE)

Fit_hbc_Summ <- data.frame(summary(Fit_hbc))
write.csv(Fit_hbc_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_hbc.csv",row.names=TRUE)

Fit_hbu_Summ <- data.frame(summary(Fit_hbu))
write.csv(Fit_hbu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_hbu.csv",row.names=TRUE)

Fit_hcu_Summ <- data.frame(summary(Fit_hcu))
write.csv(Fit_hcu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_hcu.csv",row.names=TRUE)

Fit_bcu_Summ <- data.frame(summary(Fit_bcu))
write.csv(Fit_bcu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_bcu.csv",row.names=TRUE)

Fit_rOhb_Summ <- data.frame(summary(Fit_rOhb))
write.csv(Fit_rOhb_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rOhb.csv",row.names=TRUE)

Fit_rOhc_Summ <- data.frame(summary(Fit_rOhc))
write.csv(Fit_rOhc_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rOhc.csv",row.names=TRUE)

Fit_rOhu_Summ <- data.frame(summary(Fit_rOhu))
write.csv(Fit_rOhu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rOhu.csv",row.names=TRUE)

Fit_rObc_Summ <- data.frame(summary(Fit_rObc))
write.csv(Fit_rObc_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rObc.csv",row.names=TRUE)

Fit_rObu_Summ <- data.frame(summary(Fit_rObu))
write.csv(Fit_rObu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rObu.csv",row.names=TRUE)

Fit_rOcu_Summ <- data.frame(summary(Fit_rOcu))
write.csv(Fit_rOcu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rOcu.csv",row.names=TRUE)

Fit_rhbc_Summ <- data.frame(summary(Fit_rhbc))
write.csv(Fit_rhbc_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rhbc.csv",row.names=TRUE)

Fit_rhbu_Summ <- data.frame(summary(Fit_rhbu))
write.csv(Fit_rhbu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rhbu.csv",row.names=TRUE)

Fit_rhcu_Summ <- data.frame(summary(Fit_rhcu))
write.csv(Fit_rhcu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rhcu.csv",row.names=TRUE)

Fit_rbcu_Summ <- data.frame(summary(Fit_rbcu))
write.csv(Fit_rbcu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rbcu.csv",row.names=TRUE)

Fit_Ohbc_Summ <- data.frame(summary(Fit_Ohbc))
write.csv(Fit_Ohbc_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_Ohbc.csv",row.names=TRUE)

Fit_Ohbu_Summ <- data.frame(summary(Fit_Ohbu))
write.csv(Fit_Ohbu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_Ohbu.csv",row.names=TRUE)

Fit_Ohcu_Summ <- data.frame(summary(Fit_Ohcu))
write.csv(Fit_Ohcu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_Ohcu.csv",row.names=TRUE)

Fit_Obcu_Summ <- data.frame(summary(Fit_Obcu))
write.csv(Fit_Obcu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_Obcu.csv",row.names=TRUE)

Fit_hbcu_Summ <- data.frame(summary(Fit_hbcu))
write.csv(Fit_hbcu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_hbcu.csv",row.names=TRUE)

Fit_rOhbc_Summ <- data.frame(summary(Fit_rOhbc))
write.csv(Fit_rOhbc_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rOhbc.csv",row.names=TRUE)

Fit_rOhbu_Summ <- data.frame(summary(Fit_rOhbu))
write.csv(Fit_rOhbu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rOhbu.csv",row.names=TRUE)

Fit_rOhcu_Summ <- data.frame(summary(Fit_rOhcu))
write.csv(Fit_rOhcu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rOhcu.csv",row.names=TRUE)

Fit_rObcu_Summ <- data.frame(summary(Fit_rObcu))
write.csv(Fit_rObcu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rObcu.csv",row.names=TRUE)

Fit_rhbcu_Summ <- data.frame(summary(Fit_rhbcu))
write.csv(Fit_rhbcu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_rhbcu.csv",row.names=TRUE)

Fit_Ohbcu_Summ <- data.frame(summary(Fit_Ohbcu))
write.csv(Fit_Ohbcu_Summ,"/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/csv_Files/Global_Fit_Ohbcu.csv",row.names=TRUE)

#### PDF Outputs ####

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit.pdf",onefile=TRUE)
plot(Fit)
dev.off()

#### 1 ####

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_r.pdf",onefile=TRUE)
plot(Fit_r)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_O.pdf",onefile=TRUE)
plot(Fit_O)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_h.pdf",onefile=TRUE)
plot(Fit_h)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_b.pdf",onefile=TRUE)
plot(Fit_b)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_c.pdf",onefile=TRUE)
plot(Fit_c)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_u.pdf",onefile=TRUE)
plot(Fit_u)
dev.off()

#### 2 ####

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rO.pdf",onefile=TRUE)
plot(Fit_rO)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rh.pdf",onefile=TRUE)
plot(Fit_rh)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rb.pdf",onefile=TRUE)
plot(Fit_rb)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rc.pdf",onefile=TRUE)
plot(Fit_rc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_ru.pdf",onefile=TRUE)
plot(Fit_ru)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_Oh.pdf",onefile=TRUE)
plot(Fit_Oh)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_Ob.pdf",onefile=TRUE)
plot(Fit_Ob)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_Oc.pdf",onefile=TRUE)
plot(Fit_Oc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_Ou.pdf",onefile=TRUE)
plot(Fit_Ou)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_hb.pdf",onefile=TRUE)
plot(Fit_hb)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_hc.pdf",onefile=TRUE)
plot(Fit_hc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_hu.pdf",onefile=TRUE)
plot(Fit_hu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_bc.pdf",onefile=TRUE)
plot(Fit_bc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_bu.pdf",onefile=TRUE)
plot(Fit_bu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_cu.pdf",onefile=TRUE)
plot(Fit_cu)
dev.off()

#### 3 ####

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rOh.pdf",onefile=TRUE)
plot(Fit_rOh)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rOb.pdf",onefile=TRUE)
plot(Fit_rOb)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rOc.pdf",onefile=TRUE)
plot(Fit_rOc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rOu.pdf",onefile=TRUE)
plot(Fit_rOu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rhb.pdf",onefile=TRUE)
plot(Fit_rhb)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rhc.pdf",onefile=TRUE)
plot(Fit_rhc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rhu.pdf",onefile=TRUE)
plot(Fit_rhu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rbc.pdf",onefile=TRUE)
plot(Fit_rbc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rbu.pdf",onefile=TRUE)
plot(Fit_rbu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rcu.pdf",onefile=TRUE)
plot(Fit_rcu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_Ohb.pdf",onefile=TRUE)
plot(Fit_Ohb)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_Ohc.pdf",onefile=TRUE)
plot(Fit_Ohc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_Ohu.pdf",onefile=TRUE)
plot(Fit_Ohu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_Obc.pdf",onefile=TRUE)
plot(Fit_Obc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_Obu.pdf",onefile=TRUE)
plot(Fit_Obu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_Ocu.pdf",onefile=TRUE)
plot(Fit_Ocu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_hbc.pdf",onefile=TRUE)
plot(Fit_hbc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_hbu.pdf",onefile=TRUE)
plot(Fit_hbu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_hcu.pdf",onefile=TRUE)
plot(Fit_hcu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_bcu.pdf",onefile=TRUE)
plot(Fit_bcu)
dev.off()

#### 4 ####

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rOhb.pdf",onefile=TRUE)
plot(Fit_rOhb)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rOhc.pdf",onefile=TRUE)
plot(Fit_rOhc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rOhu.pdf",onefile=TRUE)
plot(Fit_rOhu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rObc.pdf",onefile=TRUE)
plot(Fit_rObc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rObu.pdf",onefile=TRUE)
plot(Fit_rObu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rOcu.pdf",onefile=TRUE)
plot(Fit_rOcu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rhbc.pdf",onefile=TRUE)
plot(Fit_rhbc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rhbu.pdf",onefile=TRUE)
plot(Fit_rhbu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rhcu.pdf",onefile=TRUE)
plot(Fit_rhcu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rbcu.pdf",onefile=TRUE)
plot(Fit_rbcu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_Ohbc.pdf",onefile=TRUE)
plot(Fit_Ohbc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_Ohbu.pdf",onefile=TRUE)
plot(Fit_Ohbu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_Ohcu.pdf",onefile=TRUE)
plot(Fit_Ohcu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_Obcu.pdf",onefile=TRUE)
plot(Fit_Obcu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_hbcu.pdf",onefile=TRUE)
plot(Fit_hbcu)
dev.off()

#### 5 ####

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rOhbcu.pdf",onefile=TRUE)
plot(Fit_rOhbc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rOhbu.pdf",onefile=TRUE)
plot(Fit_rOhbu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rOhcu.pdf",onefile=TRUE)
plot(Fit_rOhcu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rObcu.pdf",onefile=TRUE)
plot(Fit_rObcu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_rhbcu.pdf",onefile=TRUE)
plot(Fit_rhbcu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/JAGS_Global_Outputs/Plots/Global_Fit_Ohbcu.pdf",onefile=TRUE)
plot(Fit_Ohbcu)
dev.off()
