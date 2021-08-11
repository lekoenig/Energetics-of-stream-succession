# Project: Semi-mechanistic modeling of river metabolism recovery following storm disturbances
# Script: Model simulations
# LE Koenig
# last updated August 2021

# Load packages:
library(dplyr)         # general data cleaning and manipulation
library(ggplot2)       # create plots
library(cowplot)       # plot formatting
library(patchwork)     # plot formatting
library(grid)          # plot formatting
library(lubridate)     # format timestamps
library(dataRetrieval) # interface with NWIS 
library(tidybayes)     # easily extract draws from Bayesian models 
library(sbtools)       # interface with ScienceBase
library(rstan)         # fit stan models

source("./R/Analysis_Functions.R")
source("./R/Create_plots.R")
options(mc.cores = 2)


##===================================================================##
##     Focus on one river as a reference for data simulations        ##
##===================================================================##

storm <- read.csv("./data/out/SouthBrPotomac_storm_sub.csv",header=TRUE)

# For now, pretend light is at maximum during entire recovery trajectory:
storm$light_rel <- 1


##===================================================================##
##      Generate simulated time series by 3 different models         ##
##===================================================================##

set.seed(1013)

# 1) Simulate GPP given underlying Ricker model:
source("./R/simGPP_mod1_Ricker.R")
mod1_data <- simGPP_mod1_Ricker(df = storm, storms_n = 10,  
                                sig_obs = 0.25, sig_proc = 0.1, 
                                r = 0.4, K = 15)
# Plot simulated time series:
sim.plot(filter(mod1_data$data,storm_id==1))

# Create a list housing simulated data: 
df_mod1 <- mod1_data$data
fake_data_mod1 <- list(N = length(df_mod1$storm_id),
                       storm_N = length(unique(df_mod1$storm_id)),
                       storm_ID = df_mod1$storm_id,
                       light = df_mod1$light,
                       GPP = df_mod1$GPPsim,
                       GPP_sd = df_mod1$GPPsim_sd)


# 2) Simulate GPP given underlying GPP-ER model:
source("./R/simGPP_mod2_Odum.R")
mod2_data <- simGPP_mod2_ldensdep(df = storm, storms_n = 10,
                                  sig_obs_GPP = 0.25, sig_obs_ER = 0.5,
                                  sig_proc = 0.1,
                                  mup = 0.3,mur = 0.08, K = 100)
# Plot simulated time series:
sim.plot(filter(mod2_data$data,storm_id==1))

# Create a list housing the simulated data:
df_mod2 <- filter(mod2_data$data,storm_id == 1)
fake_data_mod2 <- list(N = length(df_mod2$storm_id),
                       light=df_mod2$light,
                       GPP=df_mod2$GPPsim,
                       ER = df_mod2$ERsim)

# 3) Biomass AR(1) model:
source("./R/simGPP_mod3_AR1.R")
mod3_data <- simGPP_mod3_biomass(df = storm, storms_n = 10,
                                 sig_obs_GPP = 0.25, sig_proc = 0.1,
                                 a = 1, b = 0.1)
# Plot simulated time series:
sim.plot(filter(mod3_data$data,storm_id==1))

# Create a list housing the simulated data:
df_mod3 <- filter(mod3_data$data,storm_id == 1)
fake_data_mod3 <- list(N = length(df_mod3$storm_id),
                       light=df_mod3$light,
                       GPP=df_mod3$GPPsim)

# 4) GPP AR(1) model:
mod4_data <- simGPP_mod3(df = storm, storms_n = 10,
                         sig_obs_GPP = 0.45, sig_proc = 0.1,
                         alpha = 0.9, phi = 0.5)
# Plot simulated time series:
sim.plot(filter(mod4_data$data,storm_id==1))

# Create a list housing the simulated data:
df_mod4 <- filter(mod4_data$data,storm_id == 1)
fake_data_mod4 <- list(N = length(df_mod4$storm_id),
                       light=df_mod4$light,
                       GPP=df_mod4$GPPsim,GPP_sd = df_mod4$GPPsim_sd)



##===================================================================##
##                Run 3 models on simulated data                     ##
##===================================================================##

# 1) Ricker model:
mod1 <- rstan::stan_model("./stan/Ricker_mod_obserr_hierarchical_reparam.stan")
fit_mod1 <- rstan::sampling(mod1, data=fake_data_mod1, iter=4000, chains=4,
                            control=list(stepsize=0.5,adapt_delta=0.99,max_treedepth=12))

# 2) GPP-ER model:
mod2 <- rstan::stan_model("./stan/GPP_ER_log_densdep.stan")
fit_mod2  <- rstan::sampling(mod2, data=fake_data_mod2, iter = 6000, chains=4, 
                             control=list(stepsize=0.05,adapt_delta=0.99,max_treedepth=12))

# 3) Biomass AR(1) model:
mod3 <- rstan::stan_model("./stan/AR1_biomass.stan")
fit_mod3  <- rstan::sampling(mod3, data=fake_data_mod3, iter = 6000, chains=4, 
                             control=list(stepsize=0.05,adapt_delta=0.99,max_treedepth=12))

# 4) GPP AR(1) model:
mod4 <- rstan::stan_model("./stan/AR1_GPP.stan")
fit_mod4  <- rstan::sampling(mod4, data=fake_data_mod4, iter = 6000, chains=4, 
                             control=list(stepsize=0.05,adapt_delta=0.99,max_treedepth=12))


##===================================================================##
##        Examine modeled posteriors & extract parameters            ##
##===================================================================##

# 1) Ricker model:

# Extract parameter posterior probability distributions:
rstan::traceplot(fit_mod1, pars= c("mu_r", "mu_b", "sigma_proc","sigma_obs"))
fit_extract1 <- rstan::extract(fit_mod1) 

# Plot parameters:
params_mod1 <- mod1_data$parameters

rplot <- post.plot(fit_extract1$mu_r) + labs(x="mu_r") + ggtitle(label ="Param: mu_r") +
         geom_vline(xintercept=params_mod1$r,color="black",lty=2) 
bplot <- post.plot(fit_extract1$mu_b) + labs(x="mu_b") + ggtitle(label ="Param: mu_b") + 
         geom_vline(xintercept=(-params_mod1$r/params_mod1$K),color="black",lty=2) 
sigplot <- post.plot(fit_extract1$sigma_proc) + labs(x="sigma_proc") + ggtitle(label ="Process error") + 
           geom_vline(xintercept=params_mod1$sig_proc,color="black",lty=2) 
obsplot <- post.plot(fit_extract1$sigma_obs) + labs(x="sigma_obs") + ggtitle(label ="Obs error") + 
           geom_vline(xintercept=params_mod1$sig_obs,color="black",lty=2) 

print(rplot + bplot + sigplot + obsplot + plot_layout(ncol=2))

# Plot posterior predictive distribution:
y <- mod1_data$data$GPPsim
yrep <- fit_extract1[["GPP_tilde"]]
samp <- sample(nrow(yrep), size = 100)
yrep <- yrep[samp, ]
bayesplot::ppc_dens_overlay(y, yrep)


# 2) GPP-ER model:

# Extract parameter posterior probability distributions:
rstan::traceplot(fit_mod2, pars= c("mup", "mur", "sigma_proc","sigma_GPP","sigma_ER"))
fit_extract2 <- rstan::extract(fit_mod2) 

# Plot parameters:
params_mod2 <- mod2_data$parameters

mupplot <- post.plot(fit_extract2$mup) + labs(x="mup") + ggtitle(label ="Param: mup") +
  geom_vline(xintercept=params_mod2$mup,color="black",lty=2) 
murplot <- post.plot(fit_extract2$mur) + labs(x="mur") + ggtitle(label ="Param: mur") + 
  geom_vline(xintercept=params_mod2$mur,color="black",lty=2) 
Kplot <- post.plot(fit_extract2$K) + labs(x="K") + ggtitle(label ="Param: K") + 
  geom_vline(xintercept=params_mod2$K,color="black",lty=2) 
sigplot <- post.plot(fit_extract2$sigma_proc) + labs(x="sigma_proc") + ggtitle(label ="Process error") + 
  geom_vline(xintercept=params_mod2$sigma_proc,color="black",lty=2) 
obsplotGPP <- post.plot(fit_extract2$sigma_GPP) + labs(x="sigma_obs_GPP") + ggtitle(label ="Obs error: GPP") + 
  geom_vline(xintercept=params_mod2$sigma_obs_GPP,color="black",lty=2) 
obsplotER <- post.plot(fit_extract2$sigma_ER) + labs(x="sigma_obs_ER") + ggtitle(label ="Obs error: ER") + 
  geom_vline(xintercept=params_mod2$sigma_obs_ER,color="black",lty=2) 

print(mupplot + murplot + Kplot + sigplot + obsplotGPP + obsplotER + plot_layout(ncol=2))

# Plot posterior predictive distributions:
y <- df_mod2$GPPsim
yrep <- fit_extract2[["GPP_tilde"]]
samp <- sample(nrow(yrep), size = 100)
yrep <- yrep[samp, ]
bayesplot::ppc_dens_overlay(y, yrep) + ggplot2::ggtitle("GPP posterior predicted")

y <- df_mod2$ERsim
yrep <- fit_extract2[["ER_tilde"]]
samp <- sample(nrow(yrep), size = 100)
yrep <- yrep[samp, ]
bayesplot::ppc_dens_overlay(y, yrep) + ggplot2::ggtitle("ER posterior predicted")

# Plot joint distributions of parameters:
pairs(fit_mod2, pars = c("sigma_GPP", "sigma_proc","sigma_ER","mup","mur")) 


# 3) Biomass AR(1) model:

# Extract parameter posterior probability distributions:
rstan::traceplot(fit_mod3, pars= c("a", "b", "sigma_proc","sigma_GPP"))
fit_extract3 <- rstan::extract(fit_mod3) 

# Plot parameters:
params_mod3 <- mod3_data$parameters

aplot <- post.plot(fit_extract3$a) + labs(x="a") + ggtitle(label ="Param: a") +
  geom_vline(xintercept=params_mod3$a,color="black",lty=2) 
bplot <- post.plot(fit_extract3$b) + labs(x="b") + ggtitle(label ="Param: b") + 
  geom_vline(xintercept=params_mod3$b,color="black",lty=2) 
sigplot <- post.plot(fit_extract3$sigma_proc) + labs(x="sigma_proc") + ggtitle(label ="Process error") + 
  geom_vline(xintercept=params_mod3$sigma_proc,color="black",lty=2) 
obsplotGPP <- post.plot(fit_extract3$sigma_GPP) + labs(x="sigma_obs_GPP") + ggtitle(label ="Obs error: GPP") + 
  geom_vline(xintercept=params_mod3$sigma_obs_GPP,color="black",lty=2) 

print(aplot + bplot + sigplot + obsplotGPP + plot_layout(ncol=2))

# Plot posterior predictive distributions:
y <- df_mod3$GPPsim
yrep <- fit_extract3[["GPP_tilde"]]
samp <- sample(nrow(yrep), size = 100)
yrep <- yrep[samp, ]
bayesplot::ppc_dens_overlay(y, yrep) + ggplot2::ggtitle("GPP posterior predicted")

# Plot joint distributions of parameters:
pairs(fit_mod3, pars = c("sigma_GPP", "sigma_proc","a","b")) 


# 4) GPP AR(1) model:

# Extract parameter posterior probability distributions:
rstan::traceplot(fit_mod4, pars= c("alpha", "phi", "sigma_proc","sigma_GPP"))
fit_extract4 <- rstan::extract(fit_mod4) 

# Plot parameters:
params_mod4 <- mod4_data$parameters

alphaplot <- post.plot(fit_extract4$alpha) + labs(x="alpha") + ggtitle(label ="Param: alpha") +
  geom_vline(xintercept=params_mod4$alpha,color="black",lty=2) 
phiplot <- post.plot(fit_extract4$phi) + labs(x="phi") + ggtitle(label ="Param: phi") + 
  geom_vline(xintercept=params_mod4$phi,color="black",lty=2) 
sigplot <- post.plot(fit_extract4$sigma_proc) + labs(x="sigma_proc") + ggtitle(label ="Process error") + 
  geom_vline(xintercept=params_mod4$sigma_proc,color="black",lty=2) 
obsplotGPP <- post.plot(fit_extract4$sigma_GPP) + labs(x="sigma_obs_GPP") + ggtitle(label ="Obs error: GPP") + 
  geom_vline(xintercept=params_mod4$sigma_obs_GPP,color="black",lty=2) 

print(alphaplot + phiplot + sigplot + obsplotGPP + plot_layout(ncol=2))

# Plot posterior predictive distributions:
y <- df_mod4$GPPsim
yrep <- fit_extract4[["GPP_tilde"]]
samp <- sample(nrow(yrep), size = 100)
yrep <- yrep[samp, ]
bayesplot::ppc_dens_overlay(y, yrep) + ggplot2::ggtitle("GPP posterior predicted")

# Plot joint distributions of parameters:
pairs(fit_mod4, pars = c("sigma_GPP", "sigma_proc","alpha","phi")) 




