
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

# use fake data where GPP saturates:
storm <- read.csv("./data/fakedata/GPP_sat_testing_fakedata.csv",header=TRUE) %>% mutate(date = as.Date(as.POSIXct(as.character(date), format="%m/%d/%y")))
storm %>% ggplot() + geom_point(aes(x=date,y=GPP_daily_mean),color="darkgreen") + theme_classic()
  
##===================================================================##
##             Generate *simulated* GPP time series                  ##
##===================================================================##

# Simulate data according to Ricker model to see if we can retrieve parameters:
set.seed(2000)

# Simulate data across how many storms?
no_storms <- 10

# Define parameter values:
r <- 0.53                                    # specific growth rate
r.true.err <- (r + rnorm(10,0,0.05))
K <- 15                                      # biomass carrying capacity (g C m-2)
K.true.err <- K + rnorm(10,0,1)
b <- -r.true.err/K.true.err                  # substitute for r/K ("b")
sigma_proc <- 0.1                            # process error (based on log biomass data so this is really a scale parameter)
sigma_obs <- rnorm(1,mean(obs_err$GPP_daily_sd,na.rm = TRUE),sd(obs_err$GPP_daily_sd,na.rm=TRUE))  # observation/measurement error
#sigma_obs <- 0.1                             # effectively "turn off" observation error
#sigma_proc <- 0.01                          # effectively "turn off" process error

# initialize time series:
dT <- 1                                             # time step
t <- seq(from=1,to=length(storm$date),by=dT)        # for now, all simulated storms are the same length
L <- storm$light_rel_smooth                         # relativized light - 3 day moving average (unitless)
GPP.pred <- rep(NA,length(t))
GPP.pred[1] <- 0.45                                 # Starting GPP; g O2 m^-2 d^-1
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
storm$GPP_sim <- ypred
storm %>% ggplot() + geom_point(aes(x=date,y=GPP_daily_mean),color="darkgreen") + geom_line(aes(x=date,y=GPP_daily_mean),alpha=.5,color="darkgreen") + theme_classic() + 
          geom_point(aes(x=date,y=GPP_sim),color="darkblue")


##===================================================================##
##         Run non-hierarchical version of Ricker model              ##
##===================================================================##

# Create a list housing the simulated data:
fake_data <- list(N=length(storm1$time),GPP=storm1$GPPsim,light=storm1$light,GPP_sd = storm1$GPP_sd)

# Run stan model:
fit <- rstan::stan("./stan/Ricker_mod_obserr_daily.stan",data=fake_data,iter=4000,chains=4,
                   control = list(stepsize = 0.5,adapt_delta=0.99,max_treedepth=10))
#shinystan::launch_shinystan(fit)

# Inspect traceplots:
rstan::traceplot(fit, pars= c("r", "b", "sigma_proc","sigma_obs"))

# Extract posterior probability distributions for parameters:
fit_extract<-rstan::extract(fit) # pulls out our mcmc chains

# Plot parameters:
rplot <- post.plot(fit_extract$r) + geom_vline(xintercept=r.true.err[1],color="black",lty=2) + labs(x="r") + ggtitle(label ="Recovery rate")
bplot <- post.plot(fit_extract$b) + geom_vline(xintercept=(b[1]),color="black",lty=2) + labs(x="b") + ggtitle(label ="Parameter equal to r/K")
sigmaplot <- post.plot(fit_extract$sigma_proc) + geom_vline(xintercept=sigma_proc,color="black",lty=2) + labs(x="sigma_proc") + ggtitle(label ="Process error")
obsplot <- post.plot(fit_extract$sigma_obs) + geom_vline(xintercept=sigma_obs,color="black",lty=2) + labs(x="sigma_obs") + ggtitle(label ="Obs error")
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
post_GPP <- data.frame(GPPmod = fit_extract$GPPmod)  %>%
            tidyr::pivot_longer(data = .,cols=starts_with("GPP"),names_to = "day", values_to = "value") %>%
            group_by(day) %>% summarize(GPPmod = median(value,na.rm=T)) %>%
            mutate(timestep = as.integer(substr(day,8,11))) %>% arrange(.,timestep) %>% 
            left_join(.,storm1[,c("time","GPPsim")],by=c("timestep"="time"))
post_GPP2 <- data.frame(GPPmod = fit_extract$GPPmod)  %>%
            tidyr::pivot_longer(data = .,cols=starts_with("GPP"),names_to = "day", values_to = "value") %>%
            mutate(timestep = as.integer(substr(day,8,11))) %>% arrange(.,timestep) %>% 
            left_join(.,storm1[,c("time","GPPsim")],by=c("timestep"="time"))

plot_grid(
post_GPP %>% ggplot() + geom_point(aes(x=GPPsim,y=GPPmod),color="black") + theme_classic() + 
             geom_abline(intercept=0,slope=1,lty=2) + labs(y=expression(GPP-mod[median~posterior]),x=expression(Simulated~(known)~GPP)),

post_GPP2 %>% ggplot() + geom_point(aes(x=timestep,y=value,color="Posterior GPP-mod")) + 
             geom_point(aes(x=timestep,y=GPPsim,color="Simulated GPP")) +
             scale_color_manual(values=c("gray","blue"),name="") + 
             theme_classic() + labs(x="timestep",y="GPP"),
ncol=2,rel_widths = c(0.5,1))


##======================================================================##
##   Run non-hierarchical version of Ricker model: model corr. params   ##
##======================================================================##

# Run stan model:
fit_corrparams <- rstan::stan("./stan/Ricker_mod_corrparams.stan",data=fake_data,iter=6000,chains=4,
                   control = list(stepsize = 0.5,adapt_delta=0.99,max_treedepth=13))

summary_corrparams <- rstan::summary(fit_corrparams, 
                      pars = c("beta[1]", "beta[2]","beta_tilde[1]","beta_tilde[2]","sigma_b[1]","sigma_b[2]",
                               "sigma_proc","sigma_obs","L"), 
                      probs = c(0.1, 0.9))$summary

# Inspect traceplots:
rstan::traceplot(fit_corrparams, pars= c("beta[1]", "beta[2]", "sigma_proc","sigma_obs"))

# Inspect pairs:
#pairs(fit_corrparams, pars = c("beta_tilde[1]", "lp__","beta_tilde[2]"))

# Transformed parameters beta are highly correlated, but *sampled parameters beta_tilde* are not:
draws_beta <- as.matrix(fit_corrparams, pars = "beta")
draws_betatilde <- as.matrix(fit_corrparams, pars = "beta_tilde")
plot_grid(
bayesplot::mcmc_scatter(draws_beta),
bayesplot::mcmc_scatter(draws_betatilde),
ncol=2)

# Extract posterior probability distributions for parameters:
fit_extract_corrparams <- rstan::extract(fit_corrparams) # pulls out our mcmc chains

# Plot parameters:
rplot <- post.plot(fit_extract_corrparams$beta[,1]) + geom_vline(xintercept=r.true.err[1],color="black",lty=2) + labs(x="r") + ggtitle(label ="Recovery rate")
bplot <- post.plot(fit_extract_corrparams$beta[,2]) + geom_vline(xintercept=(b[1]),color="black",lty=2) + labs(x="b") + ggtitle(label ="Parameter equal to r/K")
sigmaplot <- post.plot(fit_extract_corrparams$sigma_proc) + geom_vline(xintercept=sigma_proc,color="black",lty=2) + labs(x="sigma_proc") + ggtitle(label ="Process error")
obsplot <- post.plot(fit_extract_corrparams$sigma_obs) + geom_vline(xintercept=sigma_obs,color="black",lty=2) + labs(x="sigma_obs") + ggtitle(label ="Obs error")

sim_plot <- rplot + bplot + sigmaplot + obsplot + plot_layout(ncol=2)
print(sim_plot)

