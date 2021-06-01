
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

source("./R/Analysis_Functions.R")

##===================================================================##
##                  Read in Powell Center data                       ##
##===================================================================##

# Appling data repo: https://www.sciencebase.gov/catalog/item/59bff507e4b091459a5e0982

# Read in list of filtered sites:
filtered.sites <- read.csv("data/out/PC_site_data_filtered.csv",header=TRUE,stringsAsFactors = FALSE)

# Read in Powell Center dataset from Maite (contains all sites):
dailytot <- readRDS("./data/in/Appling_data/dailytot.rds")
dailytot$date <- as.POSIXct(as.character(dailytot$date), format="%Y-%m-%d")
dailytot$doy <- as.numeric(strftime(dailytot$date, format = "%j"))
dailytot$year <- lubridate::year(dailytot$date)

# Filter Powell Center dataset based on list of 35 filtered sites above:
dailytot.filtered <- filter(dailytot,dailytot$site %in% unique(filtered.sites$site))

# Create a new column for storm number:
dailytot.filtered$storm_no <- NA

# Split filtered dataset into a list based on site name:
metab.dat <- split(dailytot.filtered, dailytot.filtered$site)

# Select South Branch Potomac, WV as a test case:
pot <- metab.dat$nwis_01608500

# bring in light data:
unzip(zipfile = "./data/in/Appling_data/timeseries/nwis_01608500_timeseries.zip", exdir = "./data/in/Appling_data/timeseries/")
pot.light <- read.table(file = './data/in/Appling_data/timeseries/nwis_01608500_timeseries/nwis_01608500-ts_par_calcLatSw.tsv', sep = '\t', header = TRUE)
unlink(x = "./data/in/Appling_data/timeseries/nwis_01608500_timeseries/",recursive=T)
unlink(x = "./data/in/Appling_data/timeseries/__MACOSX/",recursive=T)

# sum light (daily):
pot.light2 <- pot.light %>% mutate(date = as.Date(DateTime)) %>% 
  mutate(par_molm2_30min = par * 60 * 30 / 1000000) %>%  # par is in umol/m2s
  group_by(date) %>% 
  summarize(par_molm2d = sum(par_molm2_30min,na.rm=T)) %>%
  mutate(rel.par = min_max_norm(par_molm2d))   

# join light data with metabolism parameter estimates:
pot <- left_join(pot,pot.light2,by="date")

# are light data in the range we would expect for mol m-2 d-1? (0-100)
quantile(pot$par_molm2d,na.rm=T)

# Pull observation error from model fits: # let's say avg. sd of posterior distribution is .5
pot_fit <- read.csv("./data/in/Appling_data/fits/nwis_01608500_30min_fit/daily.tsv.resaved.csv",header=TRUE) %>%
           mutate(date2 = as.Date(date,format="%m/%d/%y"))
obs_err <- pot_fit %>% summarize(mean_gpp_sd = mean(GPP_daily_sd,na.rm=T))


##===================================================================##
##                         Identify storms                           ##
##===================================================================##

## 1. For now just find one storm from South Branch Potomac and label it:

# Add a variable that is the change in GPP from yesterday:
pot <- pot %>% mutate(deltaGPP = abs((GPP-lag(GPP,1))/lag(GPP-1)))

# Find the storm in Jan/Feb 2012:
xstart1 <- as.POSIXct("2012-01-09")
which(pot$deltaGPP>1)[which.min(abs(as.POSIXct(pot$date[which(pot$deltaGPP>1)]) - xstart1))]
pot$storm_no[c(1505:1552)] <- 1

storm <- filter(pot,storm_no==1)
storm <- left_join(storm,pot_fit[,c("date2","GPP_daily_mean","GPP_daily_sd")],by=c("date"="date2"))

# For now, filter out days where GPP < 0:
storm <- storm[-which(storm$GPP<0),]
which(storm$GPP<0)

storm$GPP <- zoo::na.approx(storm$GPP)
storm$GPP_sd <- zoo::na.approx(storm$GPP_daily_sd)
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
K <- 12                                      # biomass carrying capacity (g C m-2)
K.true.err <- K + rnorm(10,0,1)
b <- -r.true.err/K.true.err                  # substitute for r/K ("b")
sigma_proc <- 0.1                            # process error (based on log biomass data so this is really a scale parameter)
sigma_obs <- as.numeric(obs_err)             # observation error 

# initialize time series:
dT <- 1                                             # time step
t <- seq(from=1,to=length(storm$date),by=dT)        # for now, all simulated storms are the same length
L <- storm$par_molm2d/max(storm$par_molm2d,na.rm=T) # relativized light (unitless)
GPP.pred <- rep(NA,length(t))
GPP.pred[1] <- 0.46                                 # Starting GPP; g O2 m^-2 d^-1
B <- rep(NA,length(t))
B[1] <- log(GPP.pred[1]/L[1])                               

# Create data frame to hold simulated data across j number of storms (defined by no_storms above):
sim_dat <- data.frame(storm_id = numeric(),time = integer(),light=numeric(),B = numeric(),GPP.pred = numeric())

# Simulate GPP given above parameters:
for(j in 1:no_storms){
  
  for(i in 1:length(t)){
    B[i+1] <- B[i] + r.true.err[j] + b[j]*exp(B[i]) + rnorm(1,mean=0,sd = sigma_proc)
    GPP.pred[i] <- L[i] * (exp(B[i])) + rnorm(1,mean=0,sd=sigma_obs)
  }
  
  df <- data.frame(storm_id = rep(j,length(t)), time = t,light=L,B = B[c(1:length(t))], GPPsim = GPP.pred)
  sim_dat <- rbind(sim_dat,df)  
  
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

##===================================================================##
##        STAN: Non-hierarchical version of Ricker model             ##
##===================================================================##

## Write out stan instructions:  

sink("./stan/GPPrecovery_Ricker.stan")

cat("
  
  data{
  int <lower=1> N; // number of time steps (days)
  vector[N] light; // relativized to max value
  vector[N] GPP;   // for now, simulated GPP values
  }

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
  r ~ normal(0,0.5);            // prior on growth rate, r
  b ~ normal(0,0.1);            // prior on ricker model term r/k
  sigma_proc ~ normal(0,1);     // prior on process error
  sigma_obs ~ normal(0,1);     // strong prior on observation error that corresponds w/ abs. sd on GPP posterior
  }
    
  generated quantities{
  real GPP_tilde[N];               // posterior predictive check on GPP

  for(i in 1:N){
    GPP_tilde[i] = normal_rng(light[i] * exp(B[i]),sigma_obs);
  }
  }
    
  ",fill=TRUE)
sink()

## Create a list housing the simulated data:
fake_data <- list(N=length(storm1$time),GPP=storm1$GPPsim,light=storm1$light)

## Run stan model:
fit <- rstan::stan("./stan/GPPrecovery_Ricker.stan",data=fake_data,iter=4000,chains=4,
                   control = list(stepsize = 0.5,adapt_delta=0.99,max_treedepth=12))
#shinystan::launch_shinystan(fit)

## Inspect traceplots:
traceplot(fit, pars= c("r", "b", "sigma_proc","sigma_obs"))

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






















