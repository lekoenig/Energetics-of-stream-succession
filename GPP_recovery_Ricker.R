
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
library(tidybayes)     # easily extract draws from Bayesian models 
library(sbtools)       # interface with ScienceBase
library(rstan)
options(mc.cores = 2)
rstan_options(auto_write = TRUE)

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
#pot$storm_id[c(which(pot$date=="2012-01-13"):which(pot$date=="2012-02-25"))] <- 1
pot$storm_id[c(which(pot$date=="2012-01-13"):which(pot$date=="2012-04-04"))] <- 1
storm <- filter(pot,storm_id==1) %>% left_join(.,obs_err[,c("date","GPP_daily_mean","GPP_daily_sd")],by=c("date"="date"))

# For now, filter out days where GPP < 0:
storm <- storm %>% 
         #filter(storm$GPP>0) %>%
         mutate(GPP = zoo::na.approx(GPP),
                GPP_sd = zoo::na.approx(GPP_daily_sd),
                light_rel_smooth = zoo::rollmean(light_rel,k=3,fill=NA,align="left"))
storm$light_rel_smooth[which(is.na(storm$light_rel_smooth))] <- rnorm(length(which(is.na(storm$light_rel_smooth))),mean(tail(storm$light_rel_smooth,n = 7),na.rm=TRUE),0.02)
#which(storm$GPP<0)
#head(storm)

##===================================================================##
##             Generate *simulated* GPP time series                  ##
##===================================================================##

# Simulate data according to Ricker model to see if we can retrieve parameters:

set.seed(1013)

# Create a function that simulates data according to Ricker model (and input parameter values):
simulate_Ricker <- function(no_storms,length_recovery,sigma_obs,sigma_proc,r,K){
  
  # Define parameter values:
  b <- -r/K                                   # substitute for r/K ("b")
  
  # Create data frame to hold simulated data across j number of storms (defined by no_storms above):
  sim_dat <- data.frame(storm_id = numeric(),time = integer(),light=numeric(),B = numeric(),GPP.pred = numeric())
  
  # Simulate GPP given above parameters:
  for(j in 1:no_storms){
    
    # Initialize time series:
    dT <- 1                                             # time step
    t <- seq(from=1,to=length_recovery,by=dT)           # for now, all simulated storms are the same length
    L <- storm$light_rel_smooth[c(1:length_recovery)]   # relativized light - 3 day moving average (unitless)
    GPP.pred <- rep(NA,length(t))
    GPP.pred[1] <- 0.45                                 # starting GPP; g O2 m^-2 d^-1
    GPPsd <- rnorm(length(t),mean(obs_err$GPP_daily_sd,na.rm=TRUE),sd(obs_err$GPP_daily_sd,na.rm=TRUE))          # simulate GPP observation error as we would have in Appling data release
    #GPPsd <- rnorm(length(t),0.1,0.2)
    #GPP.pred <- log(GPP.pred)
    B <- rep(NA,length(t))
    B[1] <- log(GPP.pred[1]/L[1])
    #B[1] <- GPP.pred[1] - log(L[1])
    
    for(i in 1:length(t)){
      B[i+1] <- B[i] + r + b*exp(B[i]) + rnorm(1,mean=0,sd = sigma_proc)
      GPP.pred[i] <- L[i] * (exp(B[i])) + rnorm(1,mean=0,sd = sigma_obs)
      #GPP.pred[i] <- log(L[i]) + B[i] + rnorm(1,mean=0,sd=sigma_obs)
    }
    
    df <- data.frame(storm_id = rep(j,length(t)), time = t,light=L,B = B[c(1:length(t))], GPPsim = GPP.pred,GPP_sd = GPPsd)
    sim_dat <- rbind(sim_dat,df)  
  }
  
  # Save parameter values and simulated data:
  pars <- data.frame(no_storms = no_storms,
                     length_recovery = length_recovery,
                     sigma_obs = sigma_obs,
                     sigma_proc = sigma_proc,
                     r = r,
                     K = K)
  out <- list(parameters = pars,data = sim_dat)
  
  return(out)
}

sim_dat <- simulate_Ricker(no_storms=10,length_recovery=21,sigma_obs=0.25,sigma_proc=0.1,r=0.35,K=15)

# For now, isolate one storm to test non-hierarchical version of Ricker model:
storm1 <- sim_dat$data[which(sim_dat$data$storm_id==1),]
ypred <- storm1$GPPsim

storm1 %>% ggplot() + 
           geom_point(aes(x=time,y=(GPPsim)),color="darkgreen") + 
           geom_line(aes(x=time,y=(GPPsim)),alpha=.5,color="darkgreen") + 
           labs(x="Days",y=expression(Simulated~GPP~(g~O[2]~m^-2~d^-1))) + 
           theme(axis.title=element_text(size=11),axis.text=element_text(size=10),
                 panel.border = element_rect(fill=NA),panel.background = element_blank())

# Plot population of simulated storms:                                        
#sim_dat %>% ggplot() + geom_line(aes(x=time,y=GPPsim),color="darkgreen",alpha=.5) + 
#  geom_point(aes(x=time,y=GPPsim),color="darkgreen") + facet_wrap(~ storm_id, ncol = 2,scales="free_y") +
#  labs(x=expression(Days~since~storm),y=expression(Simulated~GPP~(g~O[2]~m^-2~d^-1))) + theme_bw()

   
##===================================================================##
##         Run non-hierarchical version of Ricker model              ##
##===================================================================##

# Create a list housing the simulated data:
fake_data <- list(N=length(storm1$time),GPP=storm1$GPPsim,light=storm1$light,GPP_sd = storm1$GPP_sd)

# Run stan model:
mod <- rstan::stan_model("./stan/Ricker_mod_obserr.stan")
fit <- rstan::sampling(mod,data=fake_data,iter=2000,chains=4,
                   control = list(stepsize = 0.5,adapt_delta=0.99,max_treedepth=10))
#shinystan::launch_shinystan(fit)

# Inspect traceplots:
rstan::traceplot(fit, pars= c("r", "b", "sigma_proc","sigma_obs"))

# Pull out MCMC chains + extract posterior probability distributions for parameters:
fit_extract<-rstan::extract(fit) 

# Plot parameters:
rplot <- post.plot(fit_extract$r) + geom_vline(xintercept=sim_dat$parameters$r,color="black",lty=2) + labs(x="r") + ggtitle(label ="Recovery rate")
bplot <- post.plot(fit_extract$b) + geom_vline(xintercept=(-sim_dat$parameters$r/sim_dat$parameters$K),color="black",lty=2) + labs(x="b") + ggtitle(label ="Parameter equal to r/K")
sigmaplot <- post.plot(fit_extract$sigma_proc) + geom_vline(xintercept=sim_dat$parameters$sigma_proc,color="black",lty=2) + labs(x="sigma_proc") + ggtitle(label ="Process error")
obsplot <- post.plot(fit_extract$sigma_obs) + geom_vline(xintercept=sim_dat$parameters$sigma_obs,color="black",lty=2) + labs(x="sigma_obs") + ggtitle(label ="Obs error")
sim_plot <- rplot + bplot + sigmaplot + obsplot + plot_layout(ncol=2)
print(sim_plot)

# Plot joint distribution r + b:
bayesplot::color_scheme_set("blue")
joint_param_plot <- bayesplot::mcmc_scatter(as.array(fit),pars=c("r","b"),size=1.5,alpha=0.5) + theme_classic()
joint_param_plot_dens <- ggExtra::ggMarginal(joint_param_plot, type = "density",fill="darkblue",alpha=.25)

joint_error_plot <- bayesplot::mcmc_scatter(as.array(fit),pars=c("sigma_obs","sigma_proc"),size=1.5,alpha=0.5) + theme_classic()
joint_error_plot_dens <- ggExtra::ggMarginal(joint_error_plot, type = "density",fill="darkblue",alpha=.25)

joint_plot <- plot_grid(joint_param_plot_dens,joint_error_plot_dens,ncol=2)
print(joint_plot)

# Plot actual fit of the model to the data:
#post_GPP <- data.frame(GPPmod = fit_extract$GPPmod)  %>%
#            tidyr::pivot_longer(data = .,cols=starts_with("GPP"),names_to = "day", values_to = "value") %>%
#            group_by(day) %>% summarize(GPPmod = median(value,na.rm=T)) %>%
#            mutate(timestep = as.integer(substr(day,8,11))) %>% arrange(.,timestep) %>% 
#            left_join(.,storm1[,c("time","GPPsim")],by=c("timestep"="time"))
#post_GPP2 <- data.frame(GPPmod = fit_extract$GPPmod)  %>%
#            tidyr::pivot_longer(data = .,cols=starts_with("GPP"),names_to = "day", values_to = "value") %>%
#            mutate(timestep = as.integer(substr(day,8,11))) %>% arrange(.,timestep) %>% 
#            left_join(.,storm1[,c("time","GPPsim")],by=c("timestep"="time"))

#plot_grid(
#post_GPP %>% ggplot() + geom_point(aes(x=GPPsim,y=GPPmod),color="black") + theme_classic() + 
#             geom_abline(intercept=0,slope=1,lty=2) + labs(y=expression(GPP-mod[median~posterior]),x=expression(Simulated~(known)~GPP)),
#post_GPP2 %>% ggplot() + geom_point(aes(x=timestep,y=value,color="Posterior GPP-mod")) + 
#             geom_point(aes(x=timestep,y=GPPsim,color="Simulated GPP")) +
#             scale_color_manual(values=c("gray","blue"),name="") + 
#             theme_classic() + labs(x="timestep",y="GPP"),
#ncol=2,rel_widths = c(0.5,1))

# posterior predictive:
fit_pars <- fit_extract[ c('r','b','sigma_proc','sigma_obs')] %>% 
  purrr::map_df(as_tibble, .id = 'variable')

pp_gpp <- fit_extract[ 'GPP_tilde'] %>% 
  purrr::map_df(as_tibble, .id = 'variable') %>% 
  tidyr::gather(observation,value, -variable)

ggplot() + 
  geom_density(data = pp_gpp, aes(value,fill = 'Posterior Predictive'), alpha = 0.5) + 
  geom_density(data = storm1, aes(GPPsim, fill = 'Observed'), alpha = 0.5) + 
  scale_fill_manual(name="",values=c("darkorange","darkblue")) +
  labs(x=expression(GPP~(g~O[2]~m^-2~d^-1))) +
  theme_classic(base_size = 13) + theme(legend.position = c(0.85,0.85))



##===================================================================##
##           Run hierarchical version of Ricker model                ##
##===================================================================##

# Create a list housing the simulated data:
fake_data_pool <- list(N=length(sim_dat$data$storm_id),storm_N=length(unique(sim_dat$data$storm_id)),
                       storm_ID = sim_dat$data$storm_id,light=sim_dat$data$light,
                       GPP=sim_dat$data$GPPsim,GPP_sd = sim_dat$data$GPP_sd)

# Run stan model:
mod_pool <- rstan::stan_model("./stan/Ricker_mod_obserr_hierarchical_reparam.stan")
fit_pool <- rstan::sampling(mod_pool,data=fake_data_pool,iter=4000,chains=4,
                       control = list(stepsize = 0.5,adapt_delta=0.99,max_treedepth=10))
#shinystan::launch_shinystan(fit_pool)

# Inspect traceplots:
rstan::traceplot(fit_pool, pars= c("mu_r", "mu_b", "sigma_obs", "sigma_proc"))

# Extract posterior probability distributions for parameters:
fit_extract_pool<-rstan::extract(fit_pool) # pulls out our mcmc chains
pars <- tibble(.variable = c("mu_r","mu_b","sigma_obs","sigma_proc"),
               values=c(sim_dat$parameters$r,
                        (-sim_dat$parameters$r/sim_dat$parameters$K),
                        sim_dat$parameters$sigma_obs,
                        sim_dat$parameters$sigma_proc))

# Plot population-level parameters:
pop_parms <- fit_pool %>% gather_draws(mu_r,mu_b,sigma_obs,sigma_proc) %>% ggplot(aes(x=.value))+stat_halfeye(.width=.95)+
  geom_vline(data=pars,aes(xintercept=values),pars,color="blue") + facet_wrap(~.variable,nrow=2,scales="free_x") + theme_bw() + 
  labs(x="Parameter value")+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(),
        axis.title=element_text(size=11),axis.text=element_text(size=10),
        plot.title = element_text(size=12))+ggtitle(label ="Population level parameters",subtitle = "")
print(pop_parms)

# Plot group-level parameters:
rplot <- post.plot(fit_extract_pool$mu_r) + geom_vline(xintercept=sim_dat$parameters$r,color="black",lty=2) + labs(x="mu_r") + ggtitle(label ="Param: mu_r")
bplot <- post.plot(fit_extract_pool$mu_b) + geom_vline(xintercept=(-sim_dat$parameters$r/sim_dat$parameters$K),color="black",lty=2) + labs(x="mu_b") + ggtitle(label ="Param: mu_b")
sigmaplot <- post.plot(fit_extract_pool$sigma_proc) + geom_vline(xintercept=sim_dat$parameters$sigma_proc,color="black",lty=2) + labs(x="sigma_proc") + ggtitle(label ="Process error")
obsplot <- post.plot(fit_extract_pool$sigma_obs) + geom_vline(xintercept=sim_dat$parameters$sigma_obs,color="black",lty=2) + labs(x="sigma_obs") + ggtitle(label ="Obs error")
sim_plot_pool <- rplot + bplot + sigmaplot + obsplot + plot_layout(ncol=2)
print(sim_plot_pool)
multistorm_r <- fit_pool %>% spread_draws(r[N]) %>% ggplot(aes(x=r))+stat_halfeye(.width=.95,fill="blue",alpha=.4) + 
  geom_vline(aes(xintercept=r),data=sim_dat$parameters,color="blue")+facet_wrap(~N,nrow=2,scales="free") + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(),
                     axis.title=element_text(size=11),axis.text=element_text(size=10),
                     plot.title = element_text(size=12))+ggtitle(label ="Recovery rate",subtitle = "")
print(multistorm_r)
