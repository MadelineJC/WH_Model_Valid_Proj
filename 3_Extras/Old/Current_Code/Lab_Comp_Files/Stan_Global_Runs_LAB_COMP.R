install.packages("rstan")
install.packages("gdata")
install.packages("bayesplot")
install.packages("tidyverse")
install.packages("brms")
install.packages("MCMCvis")
install.packages("readr")

library(rstan)
library(gdata)
library(bayesplot)
library(tidyverse)
library(brms)
library(MCMCvis)
library(readr)

#### Data org. ####
NewData = read_csv('/home/maddie/NewData.csv')
N <- length(NewData$t) - 1
ts <- 1:N
y_init <- c(NewData$P[1], NewData$H[1])
y <- as.matrix(NewData[2:(N + 1), 2:3])
y <- cbind(y[ , 1], y[ , 2]); 
Data <- list(N = N, ts = ts, y_init = y_init, y = y)

#### Model runs: compile and run models, save outputs, save fit objects ####

#### Est. all ####

# Compile and run model; see summary of output
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model.stan")
# Global_Fit_Stan <- sampling(Model, data = Data,
#                             chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, 
#                             seed = 123)
# Global_Fit_Stan_Summ <- print(Global_Fit_Stan, pars=c("theta", "sigma", "z_init"),
#                               probs=c(0.1, 0.5, 0.9), digits = 3)
# # Export the output as a .csv
# Global_Fit_Stan_Summ <- data.frame(summary(Global_Fit_Stan))
# write.csv(Global_Fit_Stan_Summ,"Stan_Global_Outputs/Global_Fit_Stan.csv",row.names=TRUE)
# # Save the fit object as an .rds
# saveRDS(Global_Fit_Stan,"Stan_Global_FitObjects/Global_Fit_Stan.rds")
# 
# #### Give 1 parm ####
# 
# #### Give r ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-r.stan")
# Global_Fit_r_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 123)
# Global_Fit_r_Stan_Summ <- print(Global_Fit_r_Stan, pars=c("theta", "sigma", "z_init"),
#                                 probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_r_Stan_Summ <- data.frame(summary(Global_Fit_r_Stan))
# write.csv(Global_Fit_r_Stan_Summ,"Stan_Global_Outputs/Global_Fit_r_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_r_Stan,"Stan_Global_FitObjects/Global_Fit_r_Stan.rds")
# 
# #### Give O ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-O.stan")
# Global_Fit_O_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
# Global_Fit_O_Stan_Summ <- print(Global_Fit_O_Stan, pars=c("theta", "sigma", "z_init"),
#                                 probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_O_Stan_Summ <- data.frame(summary(Global_Fit_O_Stan))
# write.csv(Global_Fit_O_Stan_Summ,"Stan_Global_Outputs/Global_Fit_O_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_O_Stan,"Stan_Global_FitObjects/Global_Fit_O_Stan.rds")
# 
# #### Give h ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-h.stan")
# Global_Fit_h_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
# Global_Fit_h_Stan_Summ <- print(Global_Fit_h_Stan, pars=c("theta", "sigma", "z_init"),
#                                 probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_h_Stan_Summ <- data.frame(summary(Global_Fit_h_Stan))
# write.csv(Global_Fit_h_Stan_Summ,"Stan_Global_Outputs/Global_Fit_h_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_h_Stan,"Stan_Global_FitObjects/Global_Fit_h_Stan.rds")
# 
# #### Give b ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-b.stan")
# Global_Fit_b_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 123)
# Global_Fit_b_Stan_Summ <- print(Global_Fit_b_Stan, pars=c("theta", "sigma", "z_init"),
#                                 probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_b_Stan_Summ <- data.frame(summary(Global_Fit_b_Stan))
# write.csv(Global_Fit_b_Stan_Summ,"Stan_Global_Outputs/Global_Fit_b_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_b_Stan,"Stan_Global_FitObjects/Global_Fit_b_Stan.rds")
# 
# #### Give c ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-c.stan")
# Global_Fit_c_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
# Global_Fit_c_Stan_Summ <- print(Global_Fit_c_Stan, pars=c("theta", "sigma", "z_init"),
#                                 probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_c_Stan_Summ <- data.frame(summary(Global_Fit_c_Stan))
# write.csv(Global_Fit_c_Stan_Summ,"Stan_Global_Outputs/Global_Fit_c_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_c_Stan,"Stan_Global_FitObjects/Global_Fit_c_Stan.rds")
# 
# #### Give u ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-u.stan")
# Global_Fit_u_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 123)
# Global_Fit_u_Stan_Summ <- print(Global_Fit_u_Stan, pars=c("theta", "sigma", "z_init"),
#                                 probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_u_Stan_Summ <- data.frame(summary(Global_Fit_u_Stan))
# write.csv(Global_Fit_u_Stan_Summ,"Stan_Global_Outputs/Global_Fit_u_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_u_Stan,"Stan_Global_FitObjects/Global_Fit_u_Stan.rds")
# 
# #### Give 2 parms ####
# 
# #### Give rO ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rO.stan")
# Global_Fit_rO_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
# Global_Fit_rO_Stan_Summ <- print(Global_Fit_rO_Stan, pars=c("theta", "sigma", "z_init"),
#                                  probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_rO_Stan_Summ <- data.frame(summary(Global_Fit_rO_Stan))
# write.csv(Global_Fit_rO_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rO_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_rO_Stan,"Stan_Global_FitObjects/Global_Fit_rO_Stan.rds")
# 
# #### Give rh ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rh.stan")
# Global_Fit_rh_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
# Global_Fit_rh_Stan_Summ <- print(Global_Fit_rh_Stan, pars=c("theta", "sigma", "z_init"),
#                                  probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_rh_Stan_Summ <- data.frame(summary(Global_Fit_rh_Stan))
# write.csv(Global_Fit_rh_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rh_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_rh_Stan,"Stan_Global_FitObjects/Global_Fit_rh_Stan.rds")
# 
# #### Give rb ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rb.stan")
# Global_Fit_rb_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
# Global_Fit_rb_Stan_Summ <- print(Global_Fit_rb_Stan, pars=c("theta", "sigma", "z_init"),
#                                  probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_rb_Stan_Summ <- data.frame(summary(Global_Fit_rb_Stan))
# write.csv(Global_Fit_rb_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rb_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_rb_Stan,"Stan_Global_FitObjects/Global_Fit_rb_Stan.rds")
# 
# #### Give rc ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rc.stan")
# Global_Fit_rc_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
# Global_Fit_rc_Stan_Summ <- print(Global_Fit_rc_Stan, pars=c("theta", "sigma", "z_init"),
#                                  probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_rc_Stan_Summ <- data.frame(summary(Global_Fit_rc_Stan))
# write.csv(Global_Fit_rc_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rc_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_rc_Stan,"Stan_Global_FitObjects/Global_Fit_rc_Stan.rds")
# 
# #### Give ru ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-ru.stan")
# Global_Fit_ru_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
# Global_Fit_ru_Stan_Summ <- print(Global_Fit_ru_Stan, pars=c("theta", "sigma", "z_init"),
#                                  probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_ru_Stan_Summ <- data.frame(summary(Global_Fit_ru_Stan))
# write.csv(Global_Fit_ru_Stan_Summ,"Stan_Global_Outputs/Global_Fit_ru_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_ru_Stan,"Stan_Global_FitObjects/Global_Fit_ru_Stan.rds")
# 
# #### Give Oh ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-Oh.stan")
# Global_Fit_Oh_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
# Global_Fit_Oh_Stan_Summ <- print(Global_Fit_Oh_Stan, pars=c("theta", "sigma", "z_init"),
#                                  probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_Oh_Stan_Summ <- data.frame(summary(Global_Fit_Oh_Stan))
# write.csv(Global_Fit_Oh_Stan_Summ,"Stan_Global_Outputs/Global_Fit_Oh_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_Oh_Stan,"Stan_Global_FitObjects/Global_Fit_Oh_Stan.rds")
# 
# #### Give Ob ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-Ob.stan")
# Global_Fit_Ob_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
# Global_Fit_Ob_Stan_Summ <- print(Global_Fit_Ob_Stan, pars=c("theta", "sigma", "z_init"),
#                                  probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_Ob_Stan_Summ <- data.frame(summary(Global_Fit_Ob_Stan))
# write.csv(Global_Fit_Ob_Stan_Summ,"Stan_Global_Outputs/Global_Fit_Ob_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_Ob_Stan,"Stan_Global_FitObjects/Global_Fit_Ob_Stan.rds")
# 
# #### Give Oc ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-Oc.stan")
# Global_Fit_Oc_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
# Global_Fit_Oc_Stan_Summ <- print(Global_Fit_Oc_Stan, pars=c("theta", "sigma", "z_init"),
#                                  probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_Oc_Stan_Summ <- data.frame(summary(Global_Fit_Oc_Stan))
# write.csv(Global_Fit_Oc_Stan_Summ,"Stan_Global_Outputs/Global_Fit_Oc_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_Oc_Stan,"Stan_Global_FitObjects/Global_Fit_Oc_Stan.rds")
# 
# #### Give Ou ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-Ou.stan")
# Global_Fit_Ou_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
# Global_Fit_Ou_Stan_Summ <- print(Global_Fit_Ou_Stan, pars=c("theta", "sigma", "z_init"),
#                                  probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_Ou_Stan_Summ <- data.frame(summary(Global_Fit_Ou_Stan))
# write.csv(Global_Fit_Ou_Stan_Summ,"Stan_Global_Outputs/Global_Fit_Ou_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_Ou_Stan,"Stan_Global_FitObjects/Global_Fit_Ou_Stan.rds")
# 
# #### Give hb ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-hb.stan")
# Global_Fit_hb_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
# Global_Fit_hb_Stan_Summ <- print(Global_Fit_hb_Stan, pars=c("theta", "sigma", "z_init"),
#                                  probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_hb_Stan_Summ <- data.frame(summary(Global_Fit_hb_Stan))
# write.csv(Global_Fit_hb_Stan_Summ,"Stan_Global_Outputs/Global_Fit_hb_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_hb_Stan,"Stan_Global_FitObjects/Global_Fit_hb_Stan.rds")
# 
# #### Give hc ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-hc.stan")
# Global_Fit_hc_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
# Global_Fit_hc_Stan_Summ <- print(Global_Fit_hc_Stan, pars=c("theta", "sigma", "z_init"),
#                                  probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_hc_Stan_Summ <- data.frame(summary(Global_Fit_hc_Stan))
# write.csv(Global_Fit_hc_Stan_Summ,"Stan_Global_Outputs/Global_Fit_hc_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_hc_Stan,"Stan_Global_FitObjects/Global_Fit_hc_Stan.rds")
# 
# #### Give hu ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-hu.stan")
# Global_Fit_hu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
# Global_Fit_hu_Stan_Summ <- print(Global_Fit_hu_Stan, pars=c("theta", "sigma", "z_init"),
#                                  probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_hu_Stan_Summ <- data.frame(summary(Global_Fit_hu_Stan))
# write.csv(Global_Fit_hu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_hu_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_hu_Stan,"Stan_Global_FitObjects/Global_Fit_hu_Stan.rds")
# 
# #### Give bc ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-bc.stan")
# Global_Fit_bc_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
# Global_Fit_bc_Stan_Summ <- print(Global_Fit_bc_Stan, pars=c("theta", "sigma", "z_init"),
#                                  probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_bc_Stan_Summ <- data.frame(summary(Global_Fit_bc_Stan))
# write.csv(Global_Fit_bc_Stan_Summ,"Stan_Global_Outputs/Global_Fit_bc_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_bc_Stan,"Stan_Global_FitObjects/Global_Fit_bc_Stan.rds")
# 
# #### Give bu ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-bu.stan")
# Global_Fit_bu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
# Global_Fit_bu_Stan_Summ <- print(Global_Fit_bu_Stan, pars=c("theta", "sigma", "z_init"),
#                                  probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_bu_Stan_Summ <- data.frame(summary(Global_Fit_bu_Stan))
# write.csv(Global_Fit_bu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_bu_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_bu_Stan,"Stan_Global_FitObjects/Global_Fit_bu_Stan.rds")
# 
# #### Give cu ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-cu.stan")
# Global_Fit_cu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
# Global_Fit_cu_Stan_Summ <- print(Global_Fit_cu_Stan, pars=c("theta", "sigma", "z_init"),
#                                  probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_cu_Stan_Summ <- data.frame(summary(Global_Fit_cu_Stan))
# write.csv(Global_Fit_cu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_cu_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_cu_Stan,"Stan_Global_FitObjects/Global_Fit_cu_Stan.rds")
# 
# #### Give 3 parms ####
# 
# #### Give rOh ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rOh.stan")
# Global_Fit_rOh_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
# Global_Fit_rOh_Stan_Summ <- print(Global_Fit_rOh_Stan, pars=c("theta", "sigma", "z_init"),
#                                   probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_rOh_Stan_Summ <- data.frame(summary(Global_Fit_rOh_Stan))
# write.csv(Global_Fit_rOh_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rOh_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_rOh_Stan,"Stan_Global_FitObjects/Global_Fit_rOh_Stan.rds")
# 
# #### Give rOb ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rOb.stan")
# Global_Fit_rOb_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
# Global_Fit_rOb_Stan_Summ <- print(Global_Fit_rOb_Stan, pars=c("theta", "sigma", "z_init"),
#                                   probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_rOb_Stan_Summ <- data.frame(summary(Global_Fit_rOb_Stan))
# write.csv(Global_Fit_rOb_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rOb_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_rOb_Stan,"Stan_Global_FitObjects/Global_Fit_rOb_Stan.rds")
# 
# #### Give rOc ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rOc.stan")
# Global_Fit_rOc_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
# Global_Fit_rOc_Stan_Summ <- print(Global_Fit_rOc_Stan, pars=c("theta", "sigma", "z_init"),
#                                   probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_rOc_Stan_Summ <- data.frame(summary(Global_Fit_rOc_Stan))
# write.csv(Global_Fit_rOc_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rOc_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_rOc_Stan,"Stan_Global_FitObjects/Global_Fit_rOc_Stan.rds")
# 
# #### Give rOu ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rOu.stan")
# Global_Fit_rOu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
# Global_Fit_rOu_Stan_Summ <- print(Global_Fit_rOu_Stan, pars=c("theta", "sigma", "z_init"),
#                                   probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_rOu_Stan_Summ <- data.frame(summary(Global_Fit_rOu_Stan))
# write.csv(Global_Fit_rOu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rOu_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_rOu_Stan,"Stan_Global_FitObjects/Global_Fit_rOu_Stan.rds")
# 
# #### Give rhb ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rhb.stan")
# Global_Fit_rhb_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
# Global_Fit_rhb_Stan_Summ <- print(Global_Fit_rhb_Stan, pars=c("theta", "sigma", "z_init"),
#                                   probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_rhb_Stan_Summ <- data.frame(summary(Global_Fit_rhb_Stan))
# write.csv(Global_Fit_rhb_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rhb_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_rhb_Stan,"Stan_Global_FitObjects/Global_Fit_rhb_Stan.rds")
# 
# #### Give rhc ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rhc.stan")
# Global_Fit_rhc_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
# Global_Fit_rhc_Stan_Summ <- print(Global_Fit_rhc_Stan, pars=c("theta", "sigma", "z_init"),
#                                   probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_rhc_Stan_Summ <- data.frame(summary(Global_Fit_rhc_Stan))
# write.csv(Global_Fit_rhc_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rhc_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_rhc_Stan,"Stan_Global_FitObjects/Global_Fit_rhc_Stan.rds")
# 
# #### Give rhu ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rhu.stan")
# Global_Fit_rhu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
# Global_Fit_rhu_Stan_Summ <- print(Global_Fit_rhu_Stan, pars=c("theta", "sigma", "z_init"),
#                                   probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_rhu_Stan_Summ <- data.frame(summary(Global_Fit_rhu_Stan))
# write.csv(Global_Fit_rhu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rhu_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_rhu_Stan,"Stan_Global_FitObjects/Global_Fit_rhu_Stan.rds")
# 
# #### Give rbc ####
# 
# Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rbc.stan")
# Global_Fit_rbc_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
# Global_Fit_rbc_Stan_Summ <- print(Global_Fit_rbc_Stan, pars=c("theta", "sigma", "z_init"),
#                                   probs=c(0.1, 0.5, 0.9), digits = 3)
# Global_Fit_rbc_Stan_Summ <- data.frame(summary(Global_Fit_rbc_Stan))
# write.csv(Global_Fit_rbc_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rbc_Stan.csv",row.names=TRUE)
# saveRDS(Global_Fit_rbc_Stan,"Stan_Global_FitObjects/Global_Fit_rbc_Stan.rds")

#### Give rbu ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rbu.stan")
Global_Fit_rbu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_rbu_Stan_Summ <- print(Global_Fit_rbu_Stan, pars=c("theta", "sigma", "z_init"),
                                  probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_rbu_Stan_Summ <- data.frame(summary(Global_Fit_rbu_Stan))
write.csv(Global_Fit_rbu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rbu_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_rbu_Stan,"Stan_Global_FitObjects/Global_Fit_rbu_Stan.rds")

#### Give rcu ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rcu.stan")
Global_Fit_rcu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_rcu_Stan_Summ <- print(Global_Fit_rcu_Stan, pars=c("theta", "sigma", "z_init"),
                                  probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_rcu_Stan_Summ <- data.frame(summary(Global_Fit_rcu_Stan))
write.csv(Global_Fit_rcu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rcu_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_rcu_Stan,"Stan_Global_FitObjects/Global_Fit_rcu_Stan.rds")

#### Give Ohb ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-Ohb.stan")
Global_Fit_Ohb_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_Ohb_Stan_Summ <- print(Global_Fit_Ohb_Stan, pars=c("theta", "sigma", "z_init"),
                                  probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_Ohb_Stan_Summ <- data.frame(summary(Global_Fit_Ohb_Stan))
write.csv(Global_Fit_Ohb_Stan_Summ,"Stan_Global_Outputs/Global_Fit_Ohb_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_Ohb_Stan,"Stan_Global_FitObjects/Global_Fit_Ohb_Stan.rds")


#### Give Ohc ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-Ohc.stan")
Global_Fit_Ohc_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_Ohc_Stan_Summ <- print(Global_Fit_Ohc_Stan, pars=c("theta", "sigma", "z_init"),
                                  probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_Ohc_Stan_Summ <- data.frame(summary(Global_Fit_Ohc_Stan))
write.csv(Global_Fit_Ohc_Stan_Summ,"Stan_Global_Outputs/Global_Fit_Ohc_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_Ohc_Stan,"Stan_Global_FitObjects/Global_Fit_Ohc_Stan.rds")

#### Give Ohu ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-Ohu.stan")
Global_Fit_Ohu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_Ohu_Stan_Summ <- print(Global_Fit_Ohu_Stan, pars=c("theta", "sigma", "z_init"),
                                  probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_Ohu_Stan_Summ <- data.frame(summary(Global_Fit_Ohu_Stan))
write.csv(Global_Fit_Ohu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_Ohu_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_Ohu_Stan,"Stan_Global_FitObjects/Global_Fit_Ohu_Stan.rds")

#### Give Obc ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-Obc.stan")
Global_Fit_Obc_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_Obc_Stan_Summ <- print(Global_Fit_Obc_Stan, pars=c("theta", "sigma", "z_init"),
                                  probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_Obc_Stan_Summ <- data.frame(summary(Global_Fit_Obc_Stan))
write.csv(Global_Fit_Obc_Stan_Summ,"Stan_Global_Outputs/Global_Fit_Obc_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_Obc_Stan,"Stan_Global_FitObjects/Global_Fit_Obc_Stan.rds")

#### Give Obu ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-Obu.stan")
Global_Fit_Obu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_Obu_Stan_Summ <- print(Global_Fit_Obu_Stan, pars=c("theta", "sigma", "z_init"),
                                  probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_Obu_Stan_Summ <- data.frame(summary(Global_Fit_Obu_Stan))
write.csv(Global_Fit_Obu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_Obu_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_Obu_Stan,"Stan_Global_FitObjects/Global_Fit_Obu_Stan.rds")

#### Give Ocu ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-Ocu.stan")
Global_Fit_Ocu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_Ocu_Stan_Summ <- print(Global_Fit_Ocu_Stan, pars=c("theta", "sigma", "z_init"),
                                  probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_Ocu_Stan_Summ <- data.frame(summary(Global_Fit_Ocu_Stan))
write.csv(Global_Fit_Ocu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_Ocu_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_Ocu_Stan,"Stan_Global_FitObjects/Global_Fit_Ocu_Stan.rds")

#### Give hbc ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-hbc.stan")
Global_Fit_hbc_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_hbc_Stan_Summ <- print(Global_Fit_hbc_Stan, pars=c("theta", "sigma", "z_init"),
                                  probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_hbc_Stan_Summ <- data.frame(summary(Global_Fit_hbc_Stan))
write.csv(Global_Fit_hbc_Stan_Summ,"Stan_Global_Outputs/Global_Fit_hbc_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_hbc_Stan,"Stan_Global_FitObjects/Global_Fit_hbc_Stan.rds")

#### Give hbu ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-hbu.stan")
Global_Fit_hbu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_hbu_Stan_Summ <- print(Global_Fit_hbu_Stan, pars=c("theta", "sigma", "z_init"),
                                  probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_hbu_Stan_Summ <- data.frame(summary(Global_Fit_hbu_Stan))
write.csv(Global_Fit_hbu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_hbu_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_hbu_Stan,"Stan_Global_FitObjects/Global_Fit_hbu_Stan.rds")

#### Give hcu ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-hcu.stan")
Global_Fit_hcu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_hcu_Stan_Summ <- print(Global_Fit_hcu_Stan, pars=c("theta", "sigma", "z_init"),
                                  probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_hcu_Stan_Summ <- data.frame(summary(Global_Fit_hcu_Stan))
write.csv(Global_Fit_hcu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_hcu_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_hcu_Stan,"Stan_Global_FitObjects/Global_Fit_hcu_Stan.rds")

#### Give bcu ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-bcu.stan")
Global_Fit_bcu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_bcu_Stan_Summ <- print(Global_Fit_bcu_Stan, pars=c("theta", "sigma", "z_init"),
                                  probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_bcu_Stan_Summ <- data.frame(summary(Global_Fit_bcu_Stan))
write.csv(Global_Fit_bcu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_bcu_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_bcu_Stan,"Stan_Global_FitObjects/Global_Fit_bcu_Stan.rds")

#### Give 4 parms ####

#### Give rOhb ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rOhb.stan")
Global_Fit_rOhb_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_rOhb_Stan_Summ <- print(Global_Fit_rOhb_Stan, pars=c("theta", "sigma", "z_init"),
                                   probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_rOhb_Stan_Summ <- data.frame(summary(Global_Fit_rOhb_Stan))
write.csv(Global_Fit_rOhb_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rOhb_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_rOhb_Stan,"Stan_Global_FitObjects/Global_Fit_rOhb_Stan.rds")

#### Give rOhc ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rOhc.stan")
Global_Fit_rOhc_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_rOhc_Stan_Summ <- print(Global_Fit_rOhc_Stan, pars=c("theta", "sigma", "z_init"),
                                   probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_rOhc_Stan_Summ <- data.frame(summary(Global_Fit_rOhc_Stan))
write.csv(Global_Fit_rOhc_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rOhc_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_rOhc_Stan,"Stan_Global_FitObjects/Global_Fit_rOhc_Stan.rds")

#### Give rOhu ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rOhu.stan")
Global_Fit_rOhu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_rOhu_Stan_Summ <- print(Global_Fit_rOhu_Stan, pars=c("theta", "sigma", "z_init"),
                                   probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_rOhu_Stan_Summ <- data.frame(summary(Global_Fit_rOhu_Stan))
write.csv(Global_Fit_rOhu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rOhu_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_rOhu_Stan,"Stan_Global_FitObjects/Global_Fit_rOhu_Stan.rds")

#### Give rObc ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rObc.stan")
Global_Fit_rObc_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_rObc_Stan_Summ <- print(Global_Fit_rObc_Stan, pars=c("theta", "sigma", "z_init"),
                                   probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_rObc_Stan_Summ <- data.frame(summary(Global_Fit_rObc_Stan))
write.csv(Global_Fit_rObc_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rObc_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_rObc_Stan,"Stan_Global_FitObjects/Global_Fit_rObc_Stan.rds")

#### Give rObu ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rObu.stan")
Global_Fit_rObu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_rObu_Stan_Summ <- print(Global_Fit_rObu_Stan, pars=c("theta", "sigma", "z_init"),
                                   probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_rObu_Stan_Summ <- data.frame(summary(Global_Fit_rObu_Stan))
write.csv(Global_Fit_rObu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rObu_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_rObu_Stan,"Stan_Global_FitObjects/Global_Fit_rObu_Stan.rds")

#### Give rOcu ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rOcu.stan")
Global_Fit_rOcu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_rOcu_Stan_Summ <- print(Global_Fit_rOcu_Stan, pars=c("theta", "sigma", "z_init"),
                                   probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_rOcu_Stan_Summ <- data.frame(summary(Global_Fit_rOcu_Stan))
write.csv(Global_Fit_rOcu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rOcu_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_rOcu_Stan,"Stan_Global_FitObjects/Global_Fit_rOcu_Stan.rds")

#### Give rhbc ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rhbc.stan")
Global_Fit_rhbc_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_rhbc_Stan_Summ <- print(Global_Fit_rhbc_Stan, pars=c("theta", "sigma", "z_init"),
                                   probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_rhbc_Stan_Summ <- data.frame(summary(Global_Fit_rhbc_Stan))
write.csv(Global_Fit_rhbc_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rhbc_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_rhbc_Stan,"Stan_Global_FitObjects/Global_Fit_rhbc_Stan.rds")

#### Give rhbu ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rhbu.stan")
Global_Fit_rhbu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_rhbu_Stan_Summ <- print(Global_Fit_rhbu_Stan, pars=c("theta", "sigma", "z_init"),
                                   probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_rhbu_Stan_Summ <- data.frame(summary(Global_Fit_rhbu_Stan))
write.csv(Global_Fit_rhbu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rhbu_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_rhbu_Stan,"Stan_Global_FitObjects/Global_Fit_rhbu_Stan.rds")

#### Give rhcu ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rhcu.stan")
Global_Fit_rhcu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_rhcu_Stan_Summ <- print(Global_Fit_rhcu_Stan, pars=c("theta", "sigma", "z_init"),
                                   probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_rhcu_Stan_Summ <- data.frame(summary(Global_Fit_rhcu_Stan))
write.csv(Global_Fit_rhcu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rhcu_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_rhcu_Stan,"Stan_Global_FitObjects/Global_Fit_rhcu_Stan.rds")

#### Give rbcu ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rbcu.stan")
Global_Fit_rbcu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_rbcu_Stan_Summ <- print(Global_Fit_rbcu_Stan, pars=c("theta", "sigma", "z_init"),
                                   probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_rbcu_Stan_Summ <- data.frame(summary(Global_Fit_rbcu_Stan))
write.csv(Global_Fit_rbcu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rbcu_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_rbcu_Stan,"Stan_Global_FitObjects/Global_Fit_rbcu_Stan.rds")

#### Give Ohbc ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-Ohbc.stan")
Global_Fit_Ohbc_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_Ohbc_Stan_Summ <- print(Global_Fit_Ohbc_Stan, pars=c("theta", "sigma", "z_init"),
                                   probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_Ohbc_Stan_Summ <- data.frame(summary(Global_Fit_Ohbc_Stan))
write.csv(Global_Fit_Ohbc_Stan_Summ,"Stan_Global_Outputs/Global_Fit_Ohbc_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_Ohbc_Stan,"Stan_Global_FitObjects/Global_Fit_Ohbc_Stan.rds")

#### Give Ohbu ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-Ohbu.stan")
Global_Fit_Ohbu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_Ohbu_Stan_Summ <- print(Global_Fit_Ohbu_Stan, pars=c("theta", "sigma", "z_init"),
                                   probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_Ohbu_Stan_Summ <- data.frame(summary(Global_Fit_Ohbu_Stan))
write.csv(Global_Fit_Ohbu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_Ohbu_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_Ohbu_Stan,"Stan_Global_FitObjects/Global_Fit_Ohbu_Stan.rds")

#### Give Ohcu ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-Ohcu.stan")
Global_Fit_Ohcu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_Ohcu_Stan_Summ <- print(Global_Fit_Ohcu_Stan, pars=c("theta", "sigma", "z_init"),
                                   probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_Ohcu_Stan_Summ <- data.frame(summary(Global_Fit_Ohcu_Stan))
write.csv(Global_Fit_Ohcu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_Ohcu_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_Ohcu_Stan,"Stan_Global_FitObjects/Global_Fit_Ohcu_Stan.rds")

#### Give Obcu ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-Obcu.stan")
Global_Fit_Obcu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_Obcu_Stan_Summ <- print(Global_Fit_Obcu_Stan, pars=c("theta", "sigma", "z_init"),
                                   probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_Obcu_Stan_Summ <- data.frame(summary(Global_Fit_Obcu_Stan))
write.csv(Global_Fit_Obcu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_Obcu_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_Obcu_Stan,"Stan_Global_FitObjects/Global_Fit_Obcu_Stan.rds")

#### Give hbcu ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-hbcu.stan")
Global_Fit_hbcu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_hbcu_Stan_Summ <- print(Global_Fit_hbcu_Stan, pars=c("theta", "sigma", "z_init"),
                                   probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_hbcu_Stan_Summ <- data.frame(summary(Global_Fit_hbcu_Stan))
write.csv(Global_Fit_hbcu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_hbcu_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_hbcu_Stan,"Stan_Global_FitObjects/Global_Fit_hbcu_Stan.rds")

#### Give 5 parms ####

#### Give rOhbc ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rOhbc.stan")
Global_Fit_rOhbc_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_rOhbc_Stan_Summ <- print(Global_Fit_rOhbc_Stan, pars=c("theta", "sigma", "z_init"),
                                    probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_rOhbc_Stan_Summ <- data.frame(summary(Global_Fit_rOhbc_Stan))
write.csv(Global_Fit_rOhbc_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rOhbc_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_rOhbc_Stan,"Stan_Global_FitObjects/Global_Fit_rOhbc_Stan.rds")

#### Give rOhbu ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rOhbu.stan")
Global_Fit_rOhbu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_rOhbu_Stan_Summ <- print(Global_Fit_rOhbu_Stan, pars=c("theta", "sigma", "z_init"),
                                    probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_rOhbu_Stan_Summ <- data.frame(summary(Global_Fit_rOhbu_Stan))
write.csv(Global_Fit_rOhbu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rOhbu_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_rOhbu_Stan,"Stan_Global_FitObjects/Global_Fit_rOhbu_Stan.rds")

#### Give rOhcu ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rOhcu.stan")
Global_Fit_rOhcu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_rOhcu_Stan_Summ <- print(Global_Fit_rOhcu_Stan, pars=c("theta", "sigma", "z_init"),
                                    probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_rOhcu_Stan_Summ <- data.frame(summary(Global_Fit_rOhcu_Stan))
write.csv(Global_Fit_rOhcu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rOhcu_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_rOhcu_Stan,"Stan_Global_FitObjects/Global_Fit_rOhcu_Stan.rds")

#### Give rObcu ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rObcu.stan")
Global_Fit_rObcu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_rObcu_Stan_Summ <- print(Global_Fit_rObcu_Stan, pars=c("theta", "sigma", "z_init"),
                                    probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_rObcu_Stan_Summ <- data.frame(summary(Global_Fit_rObcu_Stan))
write.csv(Global_Fit_rObcu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rObcu_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_rObcu_Stan,"Stan_Global_FitObjects/Global_Fit_rObcu_Stan.rds")

#### Give rhbcu ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-rhbcu.stan")
Global_Fit_rhbcu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_rhbcu_Stan_Summ <- print(Global_Fit_rhbcu_Stan, pars=c("theta", "sigma", "z_init"),
                                    probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_rhbcu_Stan_Summ <- data.frame(summary(Global_Fit_rhbcu_Stan))
write.csv(Global_Fit_rhbcu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_rhbcu_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_rhbcu_Stan,"Stan_Global_FitObjects/Global_Fit_rhbcu_Stan.rds")

#### Give Ohbcu ####

Model <- stan_model("Stan_Global_Models/Global_Stan_Model_Give-Ohbcu.stan")
Global_Fit_Ohbcu_Stan <- sampling(Model, data = Data, chains = 4, iter = 20000, warmup = 10000, thin = 1, cores = 4, seed = 12)
Global_Fit_Ohbcu_Stan_Summ <- print(Global_Fit_Ohbcu_Stan, pars=c("theta", "sigma", "z_init"),
                                    probs=c(0.1, 0.5, 0.9), digits = 3)
Global_Fit_Ohbcu_Stan_Summ <- data.frame(summary(Global_Fit_Ohbcu_Stan))
write.csv(Global_Fit_Ohbcu_Stan_Summ,"Stan_Global_Outputs/Global_Fit_Ohbcu_Stan.csv",row.names=TRUE)
saveRDS(Global_Fit_Ohbcu_Stan,"Stan_Global_FitObjects/Global_Fit_Ohbcu_Stan.rds")

#### Re-introduce fit objects ####

Global_Fit_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_Stan.rds")

Global_Fit_r_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_r_Stan.rds")
Global_Fit_O_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_O_Stan.rds")
Global_Fit_h_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_h_Stan.rds")
Global_Fit_b_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_b_Stan.rds")
Global_Fit_c_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_c_Stan.rds")
Global_Fit_u_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_u_Stan.rds")

Global_Fit_rO_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rO_Stan.rds")
Global_Fit_rh_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rh_Stan.rds")
Global_Fit_rb_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rb_Stan.rds")
Global_Fit_rc_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rc_Stan.rds")
Global_Fit_ru_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_ru_Stan.rds")
Global_Fit_Oh_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_Oh_Stan.rds")
Global_Fit_Ob_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_Ob_Stan.rds")
Global_Fit_Oc_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_Oc_Stan.rds")
Global_Fit_Ou_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_Ou_Stan.rds")
Global_Fit_hb_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_hb_Stan.rds")
Global_Fit_hc_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_hc_Stan.rds")
Global_Fit_hu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_hu_Stan.rds")
Global_Fit_bc_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_bc_Stan.rds")
Global_Fit_bu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_bu_Stan.rds")
Global_Fit_cu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_cu_Stan.rds")

Global_Fit_rOh_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rOh_Stan.rds")
Global_Fit_rOb_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rOb_Stan.rds")
Global_Fit_rOc_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rOc_Stan.rds")
Global_Fit_rOu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rOu_Stan.rds")
Global_Fit_rhb_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rhb_Stan.rds")
Global_Fit_rhc_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rhc_Stan.rds")
Global_Fit_rhu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rhu_Stan.rds")
Global_Fit_rbc_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rbc_Stan.rds")
Global_Fit_rbu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rbu_Stan.rds")
Global_Fit_rcu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rcu_Stan.rds")
Global_Fit_Ohb_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_Ohb_Stan.rds")
Global_Fit_Ohc_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_Ohc_Stan.rds")
Global_Fit_Ohu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_Ohu_Stan.rds")
Global_Fit_Obc_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_Obc_Stan.rds")
Global_Fit_Obu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_Obu_Stan.rds")
Global_Fit_Ocu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_Ocu_Stan.rds")
Global_Fit_hbc_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_hbc_Stan.rds")
Global_Fit_hbu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_hbu_Stan.rds")
Global_Fit_hcu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_hcu_Stan.rds")
Global_Fit_bcu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_bcu_Stan.rds")

Global_Fit_rOhb_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rOhb_Stan.rds")
Global_Fit_rOhc_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rOhc_Stan.rds")
Global_Fit_rOhu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rOhu_Stan.rds")
Global_Fit_rObc_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rObc_Stan.rds")
Global_Fit_rObu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rObu_Stan.rds")
Global_Fit_rOcu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rOcu_Stan.rds")
Global_Fit_rhbc_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rhbc_Stan.rds")
Global_Fit_rhbu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rhbu_Stan.rds")
Global_Fit_rhcu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rhcu_Stan.rds")
Global_Fit_rbcu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rbcu_Stan.rds")
Global_Fit_Ohbc_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_Ohbc_Stan.rds")
Global_Fit_Ohbu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_Ohbu_Stan.rds")
Global_Fit_Ohcu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_Ohcu_Stan.rds")
Global_Fit_Obcu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_Obcu_Stan.rds")
Global_Fit_hbcu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_hbcu_Stan.rds")

Global_Fit_rOhbc_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rOhbc_Stan.rds")
Global_Fit_rOhbu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rOhbu_Stan.rds")
Global_Fit_rOhcu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rOhcu_Stan.rds")
Global_Fit_rObcu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rObcu_Stan.rds")
Global_Fit_rhbcu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_rhbcu_Stan.rds")
Global_Fit_Ohbcu_Stan <- readRDS("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Current_Code/Stan_Global_FitObjects/Global_Fit_Ohbcu_Stan.rds")

#### And now save as PDF ####

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]","theta[5]","theta[6]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]","theta[5]","theta[6]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_r_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_r_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]","theta[5]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_r_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_r_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]","theta[5]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_O_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_O_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]","theta[5]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_O_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_O_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]","theta[5]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_h_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_h_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]","theta[5]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_h_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_h_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]","theta[5]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_b_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_b_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]","theta[5]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_b_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_b_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]","theta[5]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_c_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_c_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]","theta[5]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_c_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_c_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]","theta[5]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_u_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_u_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]","theta[5]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_u_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_u_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]","theta[5]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rO_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rO_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rO_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_rO_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rh_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rh_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rh_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_rh_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rb_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rb_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rb_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_rb_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rc_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rc_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rc_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_rc_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_ru_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_ru_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_ru_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_ru_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Oh_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_Oh_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Oh_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_Oh_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Ob_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_Ob_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Ob_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_Ob_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Oc_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_Oc_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Oc_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_Oc_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Ou_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_Ou_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Ou_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_Ou_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_hb_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_hb_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_hb_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_hb_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_hc_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_hc_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_hc_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_hc_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_hu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_hu_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_hu_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_hu_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_bc_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_bc_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_bc_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_bc_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_bu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_bu_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_bu_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_bu_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_cu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_cu_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_cu_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_cu_Stan,pars=c("theta[1]","theta[2]","theta[3]","theta[4]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rOh_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rOh_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rOh_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_rOh_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rOb_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rOb_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rOb_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_rOb_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

##

#### rOc: THIS IS EMPTY; SOMETHING WRONG WITH MODEL!! ####
pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rOc_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rOc_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rOc_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_rOc_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rOu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rOu_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rOu_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_rOu_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rhb_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rhb_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rhb_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_rhb_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rhc_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rhc_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rhc_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_rhc_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rhu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rhu_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rhu_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_rhu_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rbc_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rbc_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rbc_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_rbc_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rbu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rbu_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rbu_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_rbu_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rcu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rcu_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rcu_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_rcu_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Ohb_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_Ohb_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Ohb_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_Ohb_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Ohc_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_Ohc_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Ohc_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_Ohc_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Ohu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_Ohu_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Ohu_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_Ohu_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

##

#### Obc: THIS IS EMPTY; SOMETHING WRONG WITH MODEL!! ####
pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Obc_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_Obc_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Obc_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_Obc_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Obu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_Obu_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Obu_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_Obu_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Ocu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_Ocu_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Ocu_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_Ocu_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_hbc_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_hbc_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_hbc_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_hbc_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

##

#### hbu: THIS IS EMPTY; SOMETHING WRONG WITH MODEL!! ####
pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_hbu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_hbu_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_hbu_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_hbu_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_hcu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_hcu_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_hcu_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_hcu_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_bcu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_bcu_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_bcu_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_bcu_Stan,pars=c("theta[1]","theta[2]","theta[3]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rOhb_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rOhb_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rOhb_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_rOhb_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rOhc_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rOhc_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rOhc_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_rOhc_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rOhu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rOhu_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rOhu_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_rOhu_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rObc_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rObc_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rObc_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_rObc_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rObu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rObu_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rObu_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_rObu_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rOcu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rOcu_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rOcu_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_rOcu_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rhbc_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rhbc_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rhbc_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_rhbc_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

##

#### rhbu: THIS IS EMPTY; SOMETHING WRONG WITH MODEL!! ####
pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rhbu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rhbu_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rhbu_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_rhbu_Stan,pars=c("theta[1]","theta[2]"))

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rhcu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rhcu_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rhcu_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_rhcu_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rbcu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rbcu_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rbcu_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_rbcu_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Ohbc_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_Ohbc_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Ohbc_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_Ohbc_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Ohbu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_Ohbu_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Ohbu_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_Ohbu_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Ohcu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_Ohcu_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Ohcu_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_Ohcu_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Obcu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_Obcu_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Obcu_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_Obcu_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

##

#### hbcu: THIS IS EMPTY; SOMETHING WRONG WITH MODEL!! ####
pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_hbcu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_hbcu_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_hbcu_Stan_Pairs.pdf",onefile=TRUE)
mcmc_pairs(Global_Fit_hbcu_Stan,pars=c("theta[1]","theta[2]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rOhbc_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rOhbc_Stan,pars=c("theta[1]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rOhbu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rOhbu_Stan,pars=c("theta[1]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rOhcu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rOhcu_Stan,pars=c("theta[1]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rObcu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rObcu_Stan,pars=c("theta[1]"))
dev.off()

##

#### rhbcu: THIS IS EMPTY; SOMETHING WRONG WITH MODEL!! ####
pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_rhbcu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_rhbcu_Stan,pars=c("theta[1]"))
dev.off()

##

pdf("/Users/mjarviscross/Desktop/GitHub/FcnalRespProj_2/Output/Stan_Global_Outputs/Plots/Global_Fit_Ohbcu_Stan_Trace.pdf",onefile=TRUE)
stan_trace(Global_Fit_Ohbcu_Stan,pars=c("theta[1]"))
dev.off()
