install.packages("rjags")
install.packages("runjags")
install.packages("modeest") 
install.packages("mc2d")
install.packages("readr")
install.packages("beepr")

library(rjags)
library(runjags)
library(modeest) 
library(mc2d)
library(readr)
library(beepr)

#### Things Needed ####

#IC=list(r=2.5, O=0.012, h=0.075, b=35, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
#monitor = c('r','O','h','b','c','u')

SmallTS_Data <- read_csv("/home/maddie/SmallTS_Data.csv")
data = SmallTS_Data
L_T = length(data$t)
C_P = as.integer(data$P); C_H = as.integer(data$H)
dt <- 0.001 # Same for all subsequnt models; small step size for "better" Euler approx.

#### Give NOTHING ####
IC=list(r=2.5, O=0.012, h=0.075, b=35, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model.R",
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

#### Give 1 parm ####

#### Give r ####
IC=list(O=0.012, h=0.075, b=35, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_r = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-r.R",
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

#### Give O ####
IC=list(r=2.5, h=0.075, b=35, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_O = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-O.R",
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

#### Give h ####
IC=list(r=2.5, O=0.012, b=35, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_h = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-h.R",
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

#### Give b ####
IC=list(r=2.5, O=0.012, h=0.075, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_b = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-b.R",
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

#### Give c ####
IC=list(r=2.5, O=0.012, h=0.075, b=35, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_c = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-c.R",
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

#### Give u ####
IC=list(r=2.5, O=0.012, h=0.075, b=35, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_u = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-u.R",
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

#### Give 2 parms ####

#### Give rO ####
IC=list(h=0.075, b=35, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rO = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rO.R",
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

#### Give rh ####
IC=list(O=0.012, b=35, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rh = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rh.R",
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

#### Give rb ####
IC=list(O=0.012, h=0.075, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rb = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rb.R",
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

#### Give rc ####
IC=list(O=0.012, h=0.075, b=35, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rc = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rc.R",
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

#### Give ru ####
IC=list(O=0.012, h=0.075, b=35, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_ru = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-ru.R",
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

#### Give Oh ####
IC=list(r=2.5, b=35, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_Oh = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-Oh.R",
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

#### Give Ob ####
IC=list(r=2.5, h=0.075, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_Ob = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-Ob.R",
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

#### Give Oc ####
IC=list(r=2.5, h=0.075, b=35, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_Oc = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-Oc.R",
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

#### Give Ou ####
IC=list(r=2.5, h=0.075, b=35, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_Ou = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-Ou.R",
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

#### Give hb ####
IC=list(r=2.5, O=0.012, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_hb = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-hb.R",
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

#### Give hc ####
IC=list(r=2.5, O=0.012, b=35, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_hc = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-hc.R",
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

#### Give hu ####
IC=list(r=2.5, O=0.012, b=35, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_hu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-hu.R",
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

#### Give bc ####
IC=list(r=2.5, O=0.012, h=0.075, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_bc = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-bc.R",
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

#### Give bu ####
IC=list(r=2.5, O=0.012, h=0.075, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_bu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-bu.R",
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

#### Give cu ####
IC=list(r=2.5, O=0.012, h=0.075, b=35); inits=list(IC,IC,IC,IC)
Global_Fit_cu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-cu.R",
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

#### Give 3 parms ####

#### Give rOh ####
IC=list(b=35, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rOh = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rOh.R",
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

#### Give rOb ####
IC=list(h=0.075, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rOb = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rOb.R",
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

#### Give rOc ####
IC=list(h=0.075, b=35, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rOc = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rOc.R",
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

#### Give rOu ####
IC=list(h=0.075, b=35, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_rOu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rOu.R",
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

#### Give rhb ####
IC=list(O=0.012, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rhb = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rhb.R",
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

#### Give rhc ####
IC=list(O=0.012, b=35, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rhc = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rhc.R",
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

#### Give rhu ####
IC=list(O=0.012, b=35, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_rhu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rhu.R",
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

#### Give rbc ####
IC=list(O=0.012, h=0.075, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rbc = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rbc.R",
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

#### Give rbu ####
IC=list(O=0.012, h=0.075, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_rbu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rbu.R",
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

#### Give rcu ####
IC=list(O=0.012, h=0.075, b=35); inits=list(IC,IC,IC,IC)
Global_Fit_rcu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rcu.R",
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

#### Give Ohb ####
IC=list(r=2.5, c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_Ohb = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-Ohb.R",
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

#### Give Ohc ####
IC=list(r=2.5, b=35, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_Ohc = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-Ohc.R",
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

#### Give Ohu ####
IC=list(r=2.5, b=35, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_Ohu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-Ohu.R",
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

#### Give Obc ####
IC=list(r=2.5, h=0.075, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_Obc = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-Obc.R",
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

#### Give Obu ####
IC=list(r=2.5, h=0.075, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_Obu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-Obu.R",
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

#### Give Ocu ####
IC=list(r=2.5, h=0.075, b=35); inits=list(IC,IC,IC,IC)
Global_Fit_Ocu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-Ocu.R",
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

#### Give hbc ####
IC=list(r=2.5, O=0.012, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_hbc = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-hbc.R",
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

#### Give hbu ####
IC=list(r=2.5, O=0.012, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_hbu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-hbu.R",
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

#### Give hcu ####
IC=list(r=2.5, O=0.012, b=35); inits=list(IC,IC,IC,IC)
Global_Fit_hcu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-hcu.R",
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

#### Give bcu ####
IC=list(r=2.5, O=0.012, h=0.075); inits=list(IC,IC,IC,IC)
Global_Fit_bcu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-bcu.R",
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

#### Give 4 parms ####

#### Give rOhb ####
IC=list(c=0.3, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rOhb = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rOhb.R",
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

#### Give rOhc ####
IC=list(b=35, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rOhc = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rOhc.R",
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

#### Give rOhu ####
IC=list(b=35, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_rOhu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rOhu.R",
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

#### Give rObc ####
IC=list(h=0.075, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rObc = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rObc.R",
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

#### Give rObu ####
IC=list(h=0.075, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_rObu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rObu.R",
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

#### Give rOcu ####
IC=list(h=0.075, b=35); inits=list(IC,IC,IC,IC)
Global_Fit_rOcu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rOcu.R",
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

#### Give rhbc ####
IC=list(O=0.012, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rhbc = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rhbc.R",
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

#### Give rhbu ####
IC=list(O=0.012, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_rhbu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rhbu.R",
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

#### Give rhcu ####
IC=list(O=0.012, b=35); inits=list(IC,IC,IC,IC)
Global_Fit_rhcu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rhcu.R",
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

#### Give rbcu ####
IC=list(O=0.012, h=0.075); inits=list(IC,IC,IC,IC)
Global_Fit_rbcu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rbcu.R",
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

#### Give Ohbc ####
IC=list(r=2.5, u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_Ohbc = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-Ohbc.R",
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

#### Give Ohbu ####
IC=list(r=2.5, c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_Ohbu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-Ohbu.R",
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

#### Give Ohcu ####
IC=list(r=2.5, b=35); inits=list(IC,IC,IC,IC)
Global_Fit_Ohcu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-Ohcu.R",
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

#### Give Obcu ####
IC=list(r=2.5, h=0.075); inits=list(IC,IC,IC,IC)
Global_Fit_Obcu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-Obcu.R",
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

#### Give hbcu ####
IC=list(r=2.5, O=0.012); inits=list(IC,IC,IC,IC)
Global_Fit_hbcu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-hbcu.R",
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

#### Give 5 parms ####

#### Give rOhbc ####
IC=list(u=0.41); inits=list(IC,IC,IC,IC)
Global_Fit_rOhbc = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rOhbc.R",
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

#### Give rOhbu ####
IC=list(c=0.3); inits=list(IC,IC,IC,IC)
Global_Fit_rOhbu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rOhbu.R",
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

#### Give rOhcu ####
IC=list(b=35); inits=list(IC,IC,IC,IC)
Global_Fit_rOhcu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rOhcu.R",
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

#### Give rObcu ####
IC=list(h=0.075); inits=list(IC,IC,IC,IC)
Global_Fit_rObcu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rObcu.R",
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

#### Give rhbcu ####
IC=list(O=0.012); inits=list(IC,IC,IC,IC)
Global_Fit_rhbcu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-rhbcu.R",
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

#### Give Ohbcu ####
IC=list(r=2.5); inits=list(IC,IC,IC,IC)
Global_Fit_Ohbcu = run.jags(
  model = "JAGS_Global_Models/JAGS_Global_Model_Give-Ohbcu.R",
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

#### Make outputs in .csv ####
# (For example; I did this for every run, but deleted most of the script because don't need it here)

Global_Fit_Summ <- data.frame(summary(Global_Fit))
write.csv(Global_Fit_Summ,"JAGS_Global_Outputs/Global_Fit.csv",row.names=TRUE)

Global_Fit_r_Summ <- data.frame(summary(Global_Fit_r))
write.csv(Global_Fit_r_Summ,"JAGS_Global_Outputs/Global_Fit_r.csv",row.names=TRUE)

Global_Fit_O_Summ <- data.frame(summary(Global_Fit_O))
write.csv(Global_Fit_O_Summ,"JAGS_Global_Outputs/Global_Fit_O.csv",row.names=TRUE)

Global_Fit_h_Summ <- data.frame(summary(Global_Fit_h))
write.csv(Global_Fit_h_Summ,"JAGS_Global_Outputs/Global_Fit_h.csv",row.names=TRUE)

Global_Fit_b_Summ <- data.frame(summary(Global_Fit_b))
write.csv(Global_Fit_b_Summ,"JAGS_Global_Outputs/Global_Fit_b.csv",row.names=TRUE)

Global_Fit_c_Summ <- data.frame(summary(Global_Fit_c))
write.csv(Global_Fit_c_Summ,"JAGS_Global_Outputs/Global_Fit_c.csv",row.names=TRUE)

Global_Fit_u_Summ <- data.frame(summary(Global_Fit_u))
write.csv(Global_Fit_u_Summ,"JAGS_Global_Outputs/Global_Fit_u.csv",row.names=TRUE)

Global_Fit_rO_Summ <- data.frame(summary(Global_Fit_rO))
write.csv(Global_Fit_rO_Summ,"JAGS_Global_Outputs/Global_Fit_rO.csv",row.names=TRUE)

Global_Fit_rh_Summ <- data.frame(summary(Global_Fit_rh))
write.csv(Global_Fit_rh_Summ,"JAGS_Global_Outputs/Global_Fit_rh.csv",row.names=TRUE)

Global_Fit_rb_Summ <- data.frame(summary(Global_Fit_rb))
write.csv(Global_Fit_rb_Summ,"JAGS_Global_Outputs/Global_Fit_rb.csv",row.names=TRUE)

Global_Fit_rc_Summ <- data.frame(summary(Global_Fit_rc))
write.csv(Global_Fit_rc_Summ,"JAGS_Global_Outputs/Global_Fit_rc.csv",row.names=TRUE)

Global_Fit_ru_Summ <- data.frame(summary(Global_Fit_ru))
write.csv(Global_Fit_ru_Summ,"JAGS_Global_Outputs/Global_Fit_ru.csv",row.names=TRUE)

Global_Fit_Oh_Summ <- data.frame(summary(Global_Fit_Oh))
write.csv(Global_Fit_Oh_Summ,"JAGS_Global_Outputs/Global_Fit_Oh.csv",row.names=TRUE)

Global_Fit_Ob_Summ <- data.frame(summary(Global_Fit_Ob))
write.csv(Global_Fit_Ob_Summ,"JAGS_Global_Outputs/Global_Fit_Ob.csv",row.names=TRUE)

Global_Fit_Oc_Summ <- data.frame(summary(Global_Fit_Oc))
write.csv(Global_Fit_Oc_Summ,"JAGS_Global_Outputs/Global_Fit_Oc.csv",row.names=TRUE)

Global_Fit_Ou_Summ <- data.frame(summary(Global_Fit_Ou))
write.csv(Global_Fit_Ou_Summ,"JAGS_Global_Outputs/Global_Fit_Ou.csv",row.names=TRUE)

Global_Fit_hb_Summ <- data.frame(summary(Global_Fit_hb))
write.csv(Global_Fit_hb_Summ,"JAGS_Global_Outputs/Global_Fit_hb.csv",row.names=TRUE)

Global_Fit_hc_Summ <- data.frame(summary(Global_Fit_hc))
write.csv(Global_Fit_hc_Summ,"JAGS_Global_Outputs/Global_Fit_hc.csv",row.names=TRUE)

Global_Fit_hu_Summ <- data.frame(summary(Global_Fit_hu))
write.csv(Global_Fit_hu_Summ,"JAGS_Global_Outputs/Global_Fit_hu.csv",row.names=TRUE)

Global_Fit_bc_Summ <- data.frame(summary(Global_Fit_bc))
write.csv(Global_Fit_bc_Summ,"JAGS_Global_Outputs/Global_Fit_bc.csv",row.names=TRUE)

Global_Fit_bu_Summ <- data.frame(summary(Global_Fit_bu))
write.csv(Global_Fit_bu_Summ,"JAGS_Global_Outputs/Global_Fit_bu.csv",row.names=TRUE)

Global_Fit_cu_Summ <- data.frame(summary(Global_Fit_cu))
write.csv(Global_Fit_cu_Summ,"JAGS_Global_Outputs/Global_Fit_cu.csv",row.names=TRUE)

Global_Fit_rOh_Summ <- data.frame(summary(Global_Fit_rOh))
write.csv(Global_Fit_rOh_Summ,"JAGS_Global_Outputs/Global_Fit_rOh.csv",row.names=TRUE)

Global_Fit_rOb_Summ <- data.frame(summary(Global_Fit_rOb))
write.csv(Global_Fit_rOb_Summ,"JAGS_Global_Outputs/Global_Fit_rOb.csv",row.names=TRUE)

Global_Fit_rOc_Summ <- data.frame(summary(Global_Fit_rOc))
write.csv(Global_Fit_rOc_Summ,"JAGS_Global_Outputs/Global_Fit_rOc.csv",row.names=TRUE)

Global_Fit_rOu_Summ <- data.frame(summary(Global_Fit_rOu))
write.csv(Global_Fit_rOu_Summ,"JAGS_Global_Outputs/Global_Fit_rOu.csv",row.names=TRUE)

Global_Fit_rhb_Summ <- data.frame(summary(Global_Fit_rhb))
write.csv(Global_Fit_rhb_Summ,"JAGS_Global_Outputs/Global_Fit_rhb.csv",row.names=TRUE)

Global_Fit_rhc_Summ <- data.frame(summary(Global_Fit_rhc))
write.csv(Global_Fit_rhc_Summ,"JAGS_Global_Outputs/Global_Fit_rhc.csv",row.names=TRUE)

Global_Fit_rhu_Summ <- data.frame(summary(Global_Fit_rhu))
write.csv(Global_Fit_rhu_Summ,"JAGS_Global_Outputs/Global_Fit_rhu.csv",row.names=TRUE)

Global_Fit_rbc_Summ <- data.frame(summary(Global_Fit_rbc))
write.csv(Global_Fit_rbc_Summ,"JAGS_Global_Outputs/Global_Fit_rbc.csv",row.names=TRUE)

Global_Fit_rbu_Summ <- data.frame(summary(Global_Fit_rbu))
write.csv(Global_Fit_rbu_Summ,"JAGS_Global_Outputs/Global_Fit_rbu.csv",row.names=TRUE)

Global_Fit_rcu_Summ <- data.frame(summary(Global_Fit_rcu))
write.csv(Global_Fit_rcu_Summ,"JAGS_Global_Outputs/Global_Fit_rcu.csv",row.names=TRUE)

Global_Fit_Ohb_Summ <- data.frame(summary(Global_Fit_Ohb))
write.csv(Global_Fit_Ohb_Summ,"JAGS_Global_Outputs/Global_Fit_Ohb.csv",row.names=TRUE)

Global_Fit_Ohc_Summ <- data.frame(summary(Global_Fit_Ohc))
write.csv(Global_Fit_Ohc_Summ,"JAGS_Global_Outputs/Global_Fit_Ohc.csv",row.names=TRUE)

Global_Fit_Ohu_Summ <- data.frame(summary(Global_Fit_Ohu))
write.csv(Global_Fit_Ohu_Summ,"JAGS_Global_Outputs/Global_Fit_Ohu.csv",row.names=TRUE)

Global_Fit_Obc_Summ <- data.frame(summary(Global_Fit_Obc))
write.csv(Global_Fit_Obc_Summ,"JAGS_Global_Outputs/Global_Fit_Obc.csv",row.names=TRUE)

Global_Fit_Obu_Summ <- data.frame(summary(Global_Fit_Obu))
write.csv(Global_Fit_Obu_Summ,"JAGS_Global_Outputs/Global_Fit_Obu.csv",row.names=TRUE)

Global_Fit_Ocu_Summ <- data.frame(summary(Global_Fit_Ocu))
write.csv(Global_Fit_Ocu_Summ,"JAGS_Global_Outputs/Global_Fit_Ocu.csv",row.names=TRUE)

Global_Fit_hbc_Summ <- data.frame(summary(Global_Fit_hbc))
write.csv(Global_Fit_hbc_Summ,"JAGS_Global_Outputs/Global_Fit_hbc.csv",row.names=TRUE)

Global_Fit_hbu_Summ <- data.frame(summary(Global_Fit_hbu))
write.csv(Global_Fit_hbu_Summ,"JAGS_Global_Outputs/Global_Fit_hbu.csv",row.names=TRUE)

Global_Fit_hcu_Summ <- data.frame(summary(Global_Fit_hcu))
write.csv(Global_Fit_hcu_Summ,"JAGS_Global_Outputs/Global_Fit_hcu.csv",row.names=TRUE)

Global_Fit_bcu_Summ <- data.frame(summary(Global_Fit_bcu))
write.csv(Global_Fit_bcu_Summ,"JAGS_Global_Outputs/Global_Fit_bcu.csv",row.names=TRUE)

Global_Fit_rOhb_Summ <- data.frame(summary(Global_Fit_rOhb))
write.csv(Global_Fit_rOhb_Summ,"JAGS_Global_Outputs/Global_Fit_rOhb.csv",row.names=TRUE)

Global_Fit_rOhc_Summ <- data.frame(summary(Global_Fit_rOhc))
write.csv(Global_Fit_rOhc_Summ,"JAGS_Global_Outputs/Global_Fit_rOhc.csv",row.names=TRUE)

Global_Fit_rOhu_Summ <- data.frame(summary(Global_Fit_rOhu))
write.csv(Global_Fit_rOhu_Summ,"JAGS_Global_Outputs/Global_Fit_rOhu.csv",row.names=TRUE)

Global_Fit_rObc_Summ <- data.frame(summary(Global_Fit_rObc))
write.csv(Global_Fit_rObc_Summ,"JAGS_Global_Outputs/Global_Fit_rObc.csv",row.names=TRUE)

Global_Fit_rObu_Summ <- data.frame(summary(Global_Fit_rObu))
write.csv(Global_Fit_rObu_Summ,"JAGS_Global_Outputs/Global_Fit_rObu.csv",row.names=TRUE)

Global_Fit_rOcu_Summ <- data.frame(summary(Global_Fit_rOcu))
write.csv(Global_Fit_rOcu_Summ,"JAGS_Global_Outputs/Global_Fit_rOcu.csv",row.names=TRUE)

Global_Fit_rhbc_Summ <- data.frame(summary(Global_Fit_rhbc))
write.csv(Global_Fit_rhbc_Summ,"JAGS_Global_Outputs/Global_Fit_rhbc.csv",row.names=TRUE)

Global_Fit_rhbu_Summ <- data.frame(summary(Global_Fit_rhbu))
write.csv(Global_Fit_rhbu_Summ,"JAGS_Global_Outputs/Global_Fit_rhbu.csv",row.names=TRUE)

Global_Fit_rhcu_Summ <- data.frame(summary(Global_Fit_rhcu))
write.csv(Global_Fit_rhcu_Summ,"JAGS_Global_Outputs/Global_Fit_rhcu.csv",row.names=TRUE)

Global_Fit_rbcu_Summ <- data.frame(summary(Global_Fit_rbcu))
write.csv(Global_Fit_rbcu_Summ,"JAGS_Global_Outputs/Global_Fit_rbcu.csv",row.names=TRUE)

Global_Fit_Ohbc_Summ <- data.frame(summary(Global_Fit_Ohbc))
write.csv(Global_Fit_Ohbc_Summ,"JAGS_Global_Outputs/Global_Fit_Ohbc.csv",row.names=TRUE)

Global_Fit_Ohbu_Summ <- data.frame(summary(Global_Fit_Ohbu))
write.csv(Global_Fit_Ohbu_Summ,"JAGS_Global_Outputs/Global_Fit_Ohbu.csv",row.names=TRUE)

Global_Fit_Ohcu_Summ <- data.frame(summary(Global_Fit_Ohcu))
write.csv(Global_Fit_Ohcu_Summ,"JAGS_Global_Outputs/Global_Fit_Ohcu.csv",row.names=TRUE)

Global_Fit_Obcu_Summ <- data.frame(summary(Global_Fit_Obcu))
write.csv(Global_Fit_Obcu_Summ,"JAGS_Global_Outputs/Global_Fit_Obcu.csv",row.names=TRUE)

Global_Fit_hbcu_Summ <- data.frame(summary(Global_Fit_hbcu))
write.csv(Global_Fit_hbcu_Summ,"JAGS_Global_Outputs/Global_Fit_hbcu.csv",row.names=TRUE)

Global_Fit_rOhbc_Summ <- data.frame(summary(Global_Fit_rOhbc))
write.csv(Global_Fit_rOhbc_Summ,"JAGS_Global_Outputs/Global_Fit_rOhbc.csv",row.names=TRUE)

Global_Fit_rOhbu_Summ <- data.frame(summary(Global_Fit_rOhbu))
write.csv(Global_Fit_rOhbu_Summ,"JAGS_Global_Outputs/Global_Fit_rOhbu.csv",row.names=TRUE)

Global_Fit_rOhcu_Summ <- data.frame(summary(Global_Fit_rOhcu))
write.csv(Global_Fit_rOhcu_Summ,"JAGS_Global_Outputs/Global_Fit_rOhcu.csv",row.names=TRUE)

Global_Fit_rObcu_Summ <- data.frame(summary(Global_Fit_rObcu))
write.csv(Global_Fit_rObcu_Summ,"JAGS_Global_Outputs/Global_Fit_rObcu.csv",row.names=TRUE)

Global_Fit_rhbcu_Summ <- data.frame(summary(Global_Fit_rhbcu))
write.csv(Global_Fit_rhbcu_Summ,"JAGS_Global_Outputs/Global_Fit_rhbcu.csv",row.names=TRUE)

Global_Fit_Ohbcu_Summ <- data.frame(summary(Global_Fit_Ohbcu))
write.csv(Global_Fit_Ohbcu_Summ,"JAGS_Global_Outputs/Global_Fit_Ohbcu.csv",row.names=TRUE)

#### Save fit objects as .rds files ####

saveRDS(Global_Fit,"Global_Fit.rds")

saveRDS(Global_Fit_r,"Global_Fit_r.rds")
saveRDS(Global_Fit_O,"Global_Fit_O.rds")
saveRDS(Global_Fit_h,"Global_Fit_h.rds")
saveRDS(Global_Fit_b,"Global_Fit_b.rds")
saveRDS(Global_Fit_c,"Global_Fit_c.rds")
saveRDS(Global_Fit_u,"Global_Fit_u.rds")

saveRDS(Global_Fit_rO,"Global_Fit_rO.rds")
saveRDS(Global_Fit_rh,"Global_Fit_rh.rds")
saveRDS(Global_Fit_rb,"Global_Fit_rb.rds")
saveRDS(Global_Fit_rc,"Global_Fit_rc.rds")
saveRDS(Global_Fit_ru,"Global_Fit_ru.rds")
saveRDS(Global_Fit_Oh,"Global_Fit_Oh.rds")
saveRDS(Global_Fit_Ob,"Global_Fit_Ob.rds")
saveRDS(Global_Fit_Oc,"Global_Fit_Oc.rds")
saveRDS(Global_Fit_Ou,"Global_Fit_Ou.rds")
saveRDS(Global_Fit_hb,"Global_Fit_hb.rds")
saveRDS(Global_Fit_hc,"Global_Fit_hc.rds")
saveRDS(Global_Fit_hu,"Global_Fit_hu.rds")
saveRDS(Global_Fit_bc,"Global_Fit_bc.rds")
saveRDS(Global_Fit_bu,"Global_Fit_bu.rds")
saveRDS(Global_Fit_cu,"Global_Fit_cu.rds")

saveRDS(Global_Fit_rOh,"Global_Fit_rOh.rds")
saveRDS(Global_Fit_rOb,"Global_Fit_rOb.rds")
saveRDS(Global_Fit_rOc,"Global_Fit_rOc.rds")
saveRDS(Global_Fit_rOu,"Global_Fit_rOu.rds")
saveRDS(Global_Fit_rhb,"Global_Fit_rhb.rds")
saveRDS(Global_Fit_rhc,"Global_Fit_rhc.rds")
saveRDS(Global_Fit_rhu,"Global_Fit_rhu.rds")
saveRDS(Global_Fit_rbc,"Global_Fit_rbc.rds")
saveRDS(Global_Fit_rbu,"Global_Fit_rbu.rds")
saveRDS(Global_Fit_rcu,"Global_Fit_rcu.rds")
saveRDS(Global_Fit_Ohb,"Global_Fit_Ohb.rds")
saveRDS(Global_Fit_Ohc,"Global_Fit_Ohc.rds")
saveRDS(Global_Fit_Ohu,"Global_Fit_Ohu.rds")
saveRDS(Global_Fit_Obc,"Global_Fit_Obc.rds")
saveRDS(Global_Fit_Obu,"Global_Fit_Obu.rds")
saveRDS(Global_Fit_Ocu,"Global_Fit_Ocu.rds")
saveRDS(Global_Fit_hbc,"Global_Fit_hbc.rds")
saveRDS(Global_Fit_hbu,"Global_Fit_hbu.rds")
saveRDS(Global_Fit_hcu,"Global_Fit_hcu.rds")
saveRDS(Global_Fit_bcu,"Global_Fit_bcu.rds")

saveRDS(Global_Fit_rOhb,"Global_Fit_rOhb.rds")
saveRDS(Global_Fit_rOhc,"Global_Fit_rOhc.rds")
saveRDS(Global_Fit_rOhu,"Global_Fit_rOhu.rds")
saveRDS(Global_Fit_rObc,"Global_Fit_rObc.rds")
saveRDS(Global_Fit_rObu,"Global_Fit_rObu.rds")
saveRDS(Global_Fit_rOcu,"Global_Fit_rOcu.rds")
saveRDS(Global_Fit_rhbc,"Global_Fit_rhbc.rds")
saveRDS(Global_Fit_rhbu,"Global_Fit_rhbu.rds")
saveRDS(Global_Fit_rhcu,"Global_Fit_rhcu.rds")
saveRDS(Global_Fit_rbcu,"Global_Fit_rbcu.rds")
saveRDS(Global_Fit_Ohbc,"Global_Fit_Ohbc.rds")
saveRDS(Global_Fit_Ohbu,"Global_Fit_Ohbu.rds")
saveRDS(Global_Fit_Ohcu,"Global_Fit_Ohcu.rds")
saveRDS(Global_Fit_Obcu,"Global_Fit_Obcu.rds")
saveRDS(Global_Fit_hbcu,"Global_Fit_hbcu.rds")

saveRDS(Global_Fit_rOhbc,"Global_Fit_rOhbc.rds")
saveRDS(Global_Fit_rOhbu,"Global_Fit_rOhbu.rds")
saveRDS(Global_Fit_rOhcu,"Global_Fit_rOhcu.rds")
saveRDS(Global_Fit_rObcu,"Global_Fit_rObcu.rds")
saveRDS(Global_Fit_rhbcu,"Global_Fit_rhbcu.rds")
saveRDS(Global_Fit_Ohbcu,"Global_Fit_Ohbcu.rds")

#### Re-introduce fit objects ####

Global_Fit <- readRDS("/Users/mjarviscross/Downloads/Global_Fit.rds")

Global_Fit_r <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_r.rds")
Global_Fit_O <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_O.rds")
Global_Fit_h <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_h.rds")
Global_Fit_b <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_b.rds")
Global_Fit_c <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_c.rds")
Global_Fit_u <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_u.rds")

Global_Fit_rO <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rO.rds")
Global_Fit_rh <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rh.rds")
Global_Fit_rb <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rb.rds")
Global_Fit_rc <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rc.rds")
Global_Fit_ru <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_ru.rds")
Global_Fit_Oh <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_Oh.rds")
Global_Fit_Ob <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_Ob.rds")
Global_Fit_Oc <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_Oc.rds")
Global_Fit_Ou <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_Ou.rds")
Global_Fit_hb <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_hb.rds")
Global_Fit_hc <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_hc.rds")
Global_Fit_hu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_hu.rds")
Global_Fit_bc <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_bc.rds")
Global_Fit_bu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_bu.rds")
Global_Fit_cu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_cu.rds")

Global_Fit_rOh <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rOh.rds")
Global_Fit_rOb <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rOb.rds")
Global_Fit_rOc <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rOc.rds")
Global_Fit_rOu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rOu.rds")
Global_Fit_rhb <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rhb.rds")
Global_Fit_rhc <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rhc.rds")
Global_Fit_rhu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rhu.rds")
Global_Fit_rbc <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rbc.rds")
Global_Fit_rbu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rbu.rds")
Global_Fit_rcu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rcu.rds")
Global_Fit_Ohb <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_Ohb.rds")
Global_Fit_Ohc <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_Ohc.rds")
Global_Fit_Ohu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_Ohu.rds")
Global_Fit_Obc <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_Obc.rds")
Global_Fit_Obu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_Obu.rds")
Global_Fit_Ocu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_Ocu.rds")
Global_Fit_hbc <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_hbc.rds")
Global_Fit_hbu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_hbu.rds")
Global_Fit_hcu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_hcu.rds")
Global_Fit_bcu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_bcu.rds")

Global_Fit_rOhb <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rOhb.rds")
Global_Fit_rOhc <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rOhc.rds")
Global_Fit_rOhu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rOhu.rds")
Global_Fit_rObc <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rObc.rds")
Global_Fit_rObu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rObu.rds")
Global_Fit_rOcu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rOcu.rds")
Global_Fit_rhbc <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rhbc.rds")
Global_Fit_rhbu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rhbu.rds")
Global_Fit_rhcu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rhcu.rds")
Global_Fit_rbcu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rbcu.rds")
Global_Fit_Ohbc <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_Ohbc.rds")
Global_Fit_Ohbu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_Ohbu.rds")
Global_Fit_Ohcu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_Ohcu.rds")
Global_Fit_Obcu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_Obcu.rds")
Global_Fit_hbcu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_hbcu.rds")

Global_Fit_rOhbc <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rOhbc.rds")
Global_Fit_rOhbu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rOhbu.rds")
Global_Fit_rOhcu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rOhcu.rds")
Global_Fit_rObcu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rObcu.rds")
Global_Fit_rhbcu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_rhbcu.rds")
Global_Fit_Ohbcu <- readRDS("/Users/mjarviscross/Downloads/Global_Fit_Ohbcu.rds")

#### And now save as PDF ####

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit.pdf",onefile=TRUE)
plot(Global_Fit)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_r.pdf",onefile=TRUE)
plot(Global_Fit_r)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_O.pdf",onefile=TRUE)
plot(Global_Fit_O)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_h.pdf",onefile=TRUE)
plot(Global_Fit_h)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_b.pdf",onefile=TRUE)
plot(Global_Fit_b)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_c.pdf",onefile=TRUE)
plot(Global_Fit_c)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_u.pdf",onefile=TRUE)
plot(Global_Fit_u)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rO.pdf",onefile=TRUE)
plot(Global_Fit_rO)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rh.pdf",onefile=TRUE)
plot(Global_Fit_rh)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rb.pdf",onefile=TRUE)
plot(Global_Fit_rb)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rc.pdf",onefile=TRUE)
plot(Global_Fit_rc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_ru.pdf",onefile=TRUE)
plot(Global_Fit_ru)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_Oh.pdf",onefile=TRUE)
plot(Global_Fit_Oh)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_Ob.pdf",onefile=TRUE)
plot(Global_Fit_Ob)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_Oc.pdf",onefile=TRUE)
plot(Global_Fit_Oc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_Ou.pdf",onefile=TRUE)
plot(Global_Fit_Ou)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_hb.pdf",onefile=TRUE)
plot(Global_Fit_hb)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_hc.pdf",onefile=TRUE)
plot(Global_Fit_hc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_hu.pdf",onefile=TRUE)
plot(Global_Fit_hu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_bc.pdf",onefile=TRUE)
plot(Global_Fit_bc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_bu.pdf",onefile=TRUE)
plot(Global_Fit_bu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_cu.pdf",onefile=TRUE)
plot(Global_Fit_cu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rOh.pdf",onefile=TRUE)
plot(Global_Fit_rOh)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rOb.pdf",onefile=TRUE)
plot(Global_Fit_rOb)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rOc.pdf",onefile=TRUE)
plot(Global_Fit_rOc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rOu.pdf",onefile=TRUE)
plot(Global_Fit_rOu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rhb.pdf",onefile=TRUE)
plot(Global_Fit_rhb)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rhc.pdf",onefile=TRUE)
plot(Global_Fit_rhc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rhu.pdf",onefile=TRUE)
plot(Global_Fit_rhu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rbc.pdf",onefile=TRUE)
plot(Global_Fit_rbc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rbu.pdf",onefile=TRUE)
plot(Global_Fit_rbu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rcu.pdf",onefile=TRUE)
plot(Global_Fit_rcu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_Ohb.pdf",onefile=TRUE)
plot(Global_Fit_Ohb)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_Ohc.pdf",onefile=TRUE)
plot(Global_Fit_Ohc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_Ohu.pdf",onefile=TRUE)
plot(Global_Fit_Ohu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_Obc.pdf",onefile=TRUE)
plot(Global_Fit_Obc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_Obu.pdf",onefile=TRUE)
plot(Global_Fit_Obu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_Ocu.pdf",onefile=TRUE)
plot(Global_Fit_Ocu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_hbc.pdf",onefile=TRUE)
plot(Global_Fit_hbc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_hbu.pdf",onefile=TRUE)
plot(Global_Fit_hbu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_hcu.pdf",onefile=TRUE)
plot(Global_Fit_hcu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_bcu.pdf",onefile=TRUE)
plot(Global_Fit_bcu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rOhb.pdf",onefile=TRUE)
plot(Global_Fit_rOhb)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rOhc.pdf",onefile=TRUE)
plot(Global_Fit_rOhc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rOhu.pdf",onefile=TRUE)
plot(Global_Fit_rOhu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rObc.pdf",onefile=TRUE)
plot(Global_Fit_rObc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rObu.pdf",onefile=TRUE)
plot(Global_Fit_rObu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rOcu.pdf",onefile=TRUE)
plot(Global_Fit_rOcu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rhbc.pdf",onefile=TRUE)
plot(Global_Fit_rhbc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rhbu.pdf",onefile=TRUE)
plot(Global_Fit_rhbu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rhcu.pdf",onefile=TRUE)
plot(Global_Fit_rhcu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rbcu.pdf",onefile=TRUE)
plot(Global_Fit_rbcu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_Ohbc.pdf",onefile=TRUE)
plot(Global_Fit_Ohbc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_Ohbu.pdf",onefile=TRUE)
plot(Global_Fit_Ohbu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_Ohcu.pdf",onefile=TRUE)
plot(Global_Fit_Ohcu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_Obcu.pdf",onefile=TRUE)
plot(Global_Fit_Obcu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_hbcu.pdf",onefile=TRUE)
plot(Global_Fit_hbcu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rOhbcu.pdf",onefile=TRUE)
plot(Global_Fit_rOhbc)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rOhbu.pdf",onefile=TRUE)
plot(Global_Fit_rOhbu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rOhcu.pdf",onefile=TRUE)
plot(Global_Fit_rOhcu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rObcu.pdf",onefile=TRUE)
plot(Global_Fit_rObcu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_rhbcu.pdf",onefile=TRUE)
plot(Global_Fit_rhbcu)
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj/Output/JAGS_Global_Outputs/Plots/Global_Fit_Ohbcu.pdf",onefile=TRUE)
plot(Global_Fit_Ohbcu)
dev.off()
