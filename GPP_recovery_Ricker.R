
## Project: Semi-mechanistic modeling of river metabolism recovery following storm disturbances
## Script: Run Ricker model - GPP recovery
## LE Koenig
## last updated June 2021

## The objective of this script is to retrieve GPP recovery parameters from simulated and empirical data using a modified Ricker model to represent biomass dynamics following storms.

# Load packages:
library(dplyr)         # general data cleaning and manipulation
library(ggplot2)       # create plots
library(cowplot)       # plot formatting
library(patchwork)     # plot formatting
library(grid)          # plot formatting
library(lubridate)     # format timestamps
library(dataRetrieval) # interface with NWIS 
library(rstan)         # fit bayesian models using stan
library(tidybayes)     # easily extract draws from Bayesian models 
library(sbtools)       # interface with ScienceBase

source("./R/Analysis_Functions.R")

##===================================================================##
##                  Read in Powell Center data                       ##
##===================================================================##

# load metabolism data (manually added):
#metab_dat_manual <- load_filtered_PC_data() 

# load filtered sites info:
sites <- read.csv("./data/out/site_data_filtered.csv",header=TRUE)

# load filtered sites metab. estimates:
metab_dat <- read.csv("./data/out/daily_predictions_filtered.csv",header=TRUE)

# Select South Branch Potomac, WV as a test case:
pot <- metab_dat[which(metab_dat$site_name=="nwis_01608500"),] %>%
       # convert light units from W m-2 to umol m-2 s-1, and relativize:
       mutate(light_umolm2s = shortwave*0.21739,
              light_rel = min_max_norm(light_umolm2s),
              date = as.Date(as.POSIXct(as.character(date), format="%Y-%m-%d")))

# are light data in the range we would expect for umol m-2 s-1? (0-100)
quantile(pot$light_umolm2s,na.rm=T)

# Pull observation error from model fits: # let's say avg. sd of posterior distribution is .5
obs_err <- readRDS("./data/out/model_fits_filtered.rds") %>% .[[which(names(.)==pot$site_name[1])]] %>%
           select(date,GPP_daily_mean,GPP_daily_sd)


##===================================================================##
##                         Identify storms                           ##
##===================================================================##

## 1. For now just find one storm from South Branch Potomac and label it:

# Add a variable that is the change in GPP from yesterday:
pot <- pot %>% mutate(deltaGPP = abs((GPP-lag(GPP,1))/lag(GPP-1)),
                      storm_id = NA)

# Find the storm in Jan/Feb 2012:
pot$storm_id[c(which(pot$date=="2012-01-13"):which(pot$date=="2012-02-25"))] <- 1
storm <- filter(pot,storm_id==1) %>% left_join(.,obs_err[,c("date","GPP_daily_mean","GPP_daily_sd")],by=c("date"="date"))

# For now, filter out days where GPP < 0:
storm <- storm %>% filter(storm$GPP>0) %>%
         mutate(GPP = zoo::na.approx(GPP),
                GPP_sd = zoo::na.approx(GPP_daily_sd),
                light_rel_smooth = zoo::rollmean(light_rel,k=3,fill=NA,align="left"))
storm$light_rel_smooth[which(is.na(storm$light_rel_smooth))] <- 0.33
which(storm$GPP<0)
head(storm)


##===================================================================##
##             Generate *simulated* GPP time series                  ##
##===================================================================##

# Simulate data according to Ricker model to see if we can retrieve parameters:

# Simulate data across how many storms?
no_storms <- 10

# Define parameter values:
r <- 0.53                                    # specific growth rate
r.true.err <- r + rnorm(10,0,0.05)
K <- 15                                      # biomass carrying capacity (g C m-2)
K.true.err <- K + rnorm(10,0,1)
b <- -r.true.err/K.true.err                  # substitute for r/K ("b")
sigma_proc <- 0.1                            # process error (based on log biomass data so this is really a scale parameter)
#sigma_obs <- 0.05
sigma_obs <- rnorm(1,mean(obs_err$GPP_daily_sd,na.rm = TRUE),sd(obs_err$GPP_daily_sd,na.rm=TRUE))  # observation/measurement error
#sigma_obs_rel <- sigma_obs/mean(obs_err$GPP_daily_mean,na.rm=TRUE) # relative observation error (for log-GPP)
  
# initialize time series:
dT <- 1                                             # time step
t <- seq(from=1,to=length(storm$date),by=dT)        # for now, all simulated storms are the same length
L <- storm$light_rel_smooth                         # relativized light - 3 day moving average (unitless)
GPP.pred <- rep(NA,length(t))
GPP.pred[1] <- 0.46                                 # Starting GPP; g O2 m^-2 d^-1
GPP.pred.l <- rep(NA,length(t))                     # logged GPP
GPP.pred.l <- log(0.46)
B <- rep(NA,length(t))
B[1] <- log(GPP.pred[1]/L[1])                               

# Create data frame to hold simulated data across j number of storms (defined by no_storms above):
sim_dat <- data.frame(storm_id = numeric(),time = integer(),light=numeric(),B = numeric(),GPP.pred = numeric())

# Simulate GPP given above parameters:
for(j in 1:no_storms){
  
  for(i in 1:length(t)){
    B[i+1] <- B[i] + r.true.err[j] + b[j]*exp(B[i]) + rnorm(1,mean=0,sd = sigma_proc)
    GPP.pred[i] <- L[i] * (exp(B[i])) + rnorm(1,mean=0,sd=sigma_obs)
    #GPP.pred.l[i] <- log(L[i]) + B[i] + rnorm(1,mean=0,sd = sigma_obs_rel) # add *relative* observation error (~sigma_obs/mean(GPP))
  }
  
  df <- data.frame(storm_id = rep(j,length(t)), time = t,light=L,B = B[c(1:length(t))], GPPsim = GPP.pred)
  df2 <- cbind(df,storm$GPP_daily_sd) %>% rename("GPP_sd" = "storm$GPP_daily_sd")
  sim_dat <- rbind(sim_dat,df2)  
  
}

# Plot simulated GPP given above parameters:
plot(sim_dat$t[which(sim_dat$storm_id==1)],sim_dat$GPPsim[which(sim_dat$storm_id==1)],type="l",main = "Simulated GPP with assigned parameters",xlab="days",ylab="GPP (g O2 m-2 d-1)")

sim_dat %>% ggplot() + geom_line(aes(x=time,y=GPPsim),color="darkgreen",alpha=.5) + 
            geom_point(aes(x=time,y=GPPsim),color="darkgreen") + facet_wrap(~ storm_id, ncol = 2,scales="free_y") +
            labs(x=expression(Days~since~storm),y=expression(Simulated~GPP~(g~O[2]~m^-2~d^-1))) + theme_bw()

# For now, isolate one storm to test non-hierarchical version of Ricker model:
storm1 <- sim_dat[which(sim_dat$storm_id==1),]
ypred <- storm1$GPPsim
print(ypred) # Is GPPsim negative?

# Does simulated GPP *look like* the data?
storm$GPP_sim <- GPP.pred
storm %>% ggplot() + geom_point(aes(x=date,y=GPP),color="darkgreen") + geom_line(aes(x=date,y=GPP),alpha=.5,color="darkgreen") + theme_classic() + 
          geom_point(aes(x=date,y=GPP.pred),color="darkblue")


##===================================================================##
##        STAN: Non-hierarchical version of Ricker model             ##
##===================================================================##

## Write out stan instructions:  

sink("./stan/GPPrecovery_Ricker.stan")

cat("
  
  data{
  int <lower=1> N;  // number of time steps (days)
  vector[N] light;  // relativized to max value
  vector[N] GPP;    // for now, simulated GPP values
  vector[N] GPP_sd; // sd estimates from Appling GPP posterior prob. distribution
  }
  
  //transformed data{

  parameters{
  
  //biomass growth parameters
  real <lower=0> r;       // specific growth rate
  real B[N];              // log-transformed latent biomass
  real b;                 // r/K term in Ricker model

  //error parameters
  real<lower=0> sigma_proc;  //process error
  real<lower=0> sigma_obs;   // obs error for gpp
  }

  transformed parameters{
  vector[N] GPPmod;
  
  for(i in 1:N){
    GPPmod[i] = light[i] * exp(B[i]);
  }
  }
    
  model{
  
  // Initialize biomass:
  if(GPP[1] < 0)
    B[1] ~ normal(log(0.01/light[1]),0.005);
  else
    B[1]~normal(log(GPP[1]/light[1]),0.005);   //small sd around initialized biomass 
  
 // Process model:
   for(i in 2:N){
       B[i] ~ normal(B[i-1] + r + b*exp(B[i-1]), sigma_proc);
   }

  // GPP observation model:
  for(i in 2:N){
      GPP[i] ~ normal(GPPmod[i],sigma_obs); 
  }
  
  // Priors on model parameters:
  r ~ normal(0,1);            // prior on growth rate, r
  b ~ normal(0,0.25);            // prior on ricker model term r/k
  sigma_proc ~ normal(0,1);     // prior on process error
  sigma_obs ~ normal(mean(GPP_sd),sd(GPP_sd));     // strong prior on observation error that corresponds w/ abs. sd on GPP posterior
  }
    
  generated quantities{
  //real GPP_tilde[N];               // posterior predictive check on GPP

  //for(i in 1:N){
  //  GPP_tilde[i] = normal_rng(light[i] * exp(B[i]),sigma_obs);
  //}    

  }
    
  ",fill=TRUE)
sink()

## Create a list housing the simulated data:
fake_data <- list(N=length(storm1$time),GPP=storm1$GPPsim,light=storm1$light,GPP_sd = storm1$GPP_sd)

## Run stan model:
fit <- rstan::stan("./stan/GPPrecovery_Ricker.stan",data=fake_data,iter=6000,chains=4,
                   control = list(stepsize = 0.5,adapt_delta=0.99,max_treedepth=12))
#shinystan::launch_shinystan(fit)

## Inspect traceplots:
traceplot(fit, pars= c("r", "b", "sigma_proc","sigma_obs"))

## Inspect pairs:
pairs(fit, pars = c("r", "lp__","b"))

## Extract posterior probability distributions for parameters:
fit_extract<-extract(fit) # pulls out our mcmc chains

## Plot parameters:
rplot <- post.plot(fit_extract$r) + geom_vline(xintercept=r.true.err[1],color="black",lty=2) + labs(x="r") + ggtitle(label ="Recovery rate")
bplot <- post.plot(fit_extract$b) + geom_vline(xintercept=(b[1]),color="black",lty=2) + labs(x="b") + ggtitle(label ="Parameter equal to r/K")
sigmaplot <- post.plot(fit_extract$sigma_proc) + geom_vline(xintercept=sigma_proc,color="black",lty=2) + labs(x="sigma_proc") + ggtitle(label ="Process error")
obsplot <- post.plot(fit_extract$sigma_obs) + geom_vline(xintercept=sigma_obs,color="black",lty=2) + labs(x="sigma_obs") + ggtitle(label ="Obs error")

sim_plot <- rplot + bplot + sigmaplot + obsplot + plot_layout(ncol=2)
print(sim_plot)

## Plot joint distribution r + b:
post_df <- data.frame(r = fit_extract$r,
                      b = fit_extract$b,
                      sigma_obs = fit_extract$sigma_obs,
                      sigma_proc = fit_extract$sigma_proc)

joint_param_plot <- ggplot(post_df, aes(r, b)) + geom_point() + theme_classic() + labs(x="r",y="b") 
joint_param_plot_dens <- ggExtra::ggMarginal(joint_param_plot, type = "density",fill="darkblue",alpha=.25)

joint_error_plot <- ggplot(post_df, aes(sigma_obs, sigma_proc)) + geom_point() + theme_classic() + labs(x="sigma_obs",y="sigma_proc") 
joint_error_plot_dens <- ggExtra::ggMarginal(joint_error_plot, type = "density",fill="darkblue",alpha=.25)

joint_plot <- plot_grid(joint_param_plot_dens,joint_error_plot_dens,ncol=2)
print(joint_plot)

