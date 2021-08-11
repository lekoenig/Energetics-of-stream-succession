## Simulate data - recovery Model 2 (Odum GPP, ER)

##==========================================================##
##     Logged biomass with density-dependence on biomass    ##
##==========================================================##

simGPP_mod2_ldensdep <- function(storms_n,df,mup,mur,K,sig_obs_GPP,sig_obs_ER,sig_proc){
  
  # Create data frame to hold simulated data across j number of storms (defined by storms_n above):
  sim_dat <- data.frame(storm_id = numeric(),
                        time = integer(),
                        light=numeric(),
                        B = numeric(),
                        GPPsim = numeric(),
                        ERsim = numeric())
  
  # Simulate GPP given above parameters:
  for(j in 1:storms_n){
    
    # Initialize data:
    dT <- 1                                         # time step
    t <- seq(from=1,to=dim(df)[1],by=dT)            # for now, all simulated storms are the same length
    light <- df$light_rel                           # relativized light 
    GPP <- df$GPP_daily_mean
    GPP_sd <- df$GPP_daily_sd
    ER <- df$ER_daily_mean*-1
    ER_sd <- df$ER_daily_sd
    
    # Initialize vectors for model output of B, GPP, and ER:
    B <- rep(NA,length(t))
    B[1] <- log(GPP[1]/(mup*light[1]))  # approx (assume B/K on day 1 is very small)
    GPP_pred <- rep(NA,length(t))
    GPP_pred[1] <- GPP[1] 
    GPP_pred_sd <- rnorm(length(t),mean(df$GPP_daily_sd,na.rm=TRUE),sd(df$GPP_daily_sd,na.rm=TRUE))          # simulate GPP observation error as we would have in Appling data release
    ER_pred <- rep(NA,length(t))
    ER_pred[1] <- ER[1]  
    ER_pred_sd <- rnorm(length(t),mean(df$ER_daily_sd,na.rm=TRUE),sd(df$ER_daily_sd,na.rm=TRUE))          # simulate ER observation error as we would have in Appling data release
    
    # Model log-biomass:
    for(i in 2:length(t)){
      
      # process model:
      B[i] <- B[i-1] + log(1+(mup*light[i])*(1-exp(B[i-1])/K)-mur) + rnorm(1,mean = 0,sd = sig_proc)
      
      # observation models:
      #GPP_pred[i] <- mup * light[i] * exp(B[i]) * (1-exp(B[i])/K) + rnorm(1,mean=0,sd = sig_obs_GPP) 
      #ER_pred[i] <- mur * exp(B[i]) + rnorm(1,mean=0,sd = sig_obs_ER) 
      
      # truncated normal distribution:
      GPP_pred[i] <- MCMCglmm::rtnorm(1, mean = mup * light[i] * exp(B[i]) * (1-exp(B[i])/K), sd = sig_obs_GPP, lower = 0)
      ER_pred[i] <- MCMCglmm::rtnorm(1, mean = mur * exp(B[i]), sd = sig_obs_ER, lower = 0)
      
    }
    
  # save simulations for j storms:
  out_data <- data.frame(storm_id = rep(j,length(t)), time = t,light=light,B = B[c(1:length(t))], 
                   GPPsim = GPP_pred,ERsim = ER_pred)
  sim_dat <- rbind(sim_dat,out_data)  
  }
      
  # Save parameter values and simulated data:
  pars <- data.frame(no_storms = storms_n,
                     sigma_obs_GPP = sig_obs_GPP,
                     sigma_obs_ER = sig_obs_ER,
                     sigma_proc = sig_proc,
                     mup = mup,
                     mur = mur,
                     K = K)
  out <- list(parameters = pars,data = sim_dat)
  
  return(out)
}



##==========================================================##
##    Unlogged biomass with density-dependence on biomass   ##
##==========================================================##

simGPP_mod2_densdep <- function(storms_n,df,mup,mur,K,sig_obs_GPP,sig_obs_ER,sig_proc){
  
  # Create data frame to hold simulated data across j number of storms (defined by no_storms above):
  sim_dat <- data.frame(storm_id = numeric(),
                        time = integer(),
                        light=numeric(),
                        B = numeric(),
                        GPP.pred = numeric(),
                        ER.pred = numeric())
  
  # Simulate GPP given above parameters:
  for(j in 1:storms_n){
    
    # Initialize data:
    dT <- 1                                         # time step
    t <- seq(from=1,to=dim(df)[1],by=dT)            # for now, all simulated storms are the same length
    light <- df$light_rel                           # relativized light 
    GPP <- df$GPP_daily_mean
    GPP_sd <- df$GPP_daily_sd
    
    # Initialize vectors for model output of B, GPP, and ER:
    B <- rep(NA,length(t))
    B[1] <- GPP[1]/light[1]
    GPP_pred <- rep(NA,length(t))
    GPP_pred[1] <- mup * light[1] * B[1] * (1-B[1]/K) # approx using t = 1
    ER_pred <- rep(NA,length(t))
    ER_pred[1] <- mur * B[1]  # approx using t = 1
    
    # Model log-biomass:
    for(i in 2:length(t)){
      
      # process model:
      B[i] <- B[i-1] + B[i-1] * mup * light[i] * (1-B[i-1]/K) - mur * B[i-1] + rnorm(1,mean = 0,sd = sig_proc)
      
      # observation models:
      GPP_pred[i] <- mup * light[i] * B[i] * (1-B[i]/K) + rnorm(1,mean=0,sd = sig_obs_GPP)
      ER_pred[i] <- mur * B[i] + rnorm(1,mean=0,sd = sig_obs_ER) 
    }
    
    # save simulations for j storms:
    out_data <- data.frame(storm_id = rep(j,length(t)), time = t,light=light,B = B, 
                           GPPsim = GPP_pred,ERsim = ER_pred)
    sim_dat <- rbind(sim_dat,out_data)  
  }
  
  # Save parameter values and simulated data:
  pars <- data.frame(no_storms = storms_n,
                     sigma_obs_GPP = sig_obs_GPP,
                     sigma_obs_ER = sig_obs_ER,
                     sigma_proc = sig_proc,
                     mup = mup,
                     mur = mur,
                     K = K)
  out <- list(parameters = pars,data = sim_dat)
  
  return(out)
}
 


##==========================================================##
##     Logged biomass - no density-dependence on biomass    ##
##==========================================================##

simGPP_mod2_l <- function(storms_n,df,mup,mur,sig_obs_GPP,sig_obs_ER,sig_proc){
  
  # Create data frame to hold simulated data across j number of storms (defined by storms_n above):
  sim_dat <- data.frame(storm_id = numeric(),
                        time = integer(),
                        light=numeric(),
                        B = numeric(),
                        GPPsim = numeric(),
                        ERsim = numeric())
  
  # Simulate GPP given above parameters:
  for(j in 1:storms_n){
    
    # Initialize data:
    dT <- 1                                         # time step
    t <- seq(from=1,to=dim(df)[1],by=dT)            # for now, all simulated storms are the same length
    light <- df$light_rel                           # relativized light 
    GPP <- df$GPP_daily_mean
    GPP_sd <- df$GPP_daily_sd
    ER <- df$ER_daily_mean*-1
    ER_sd <- df$ER_daily_sd
    
    # Initialize vectors for model output of B, GPP, and ER:
    B <- rep(NA,length(t))
    B[1] <- log(GPP[1]/(mup*light[1]))  # approx (assume B/K on day 1 is very small)
    GPP_pred <- rep(NA,length(t))
    GPP_pred[1] <- GPP[1] 
    GPP_pred_sd <- rnorm(length(t),mean(df$GPP_daily_sd,na.rm=TRUE),sd(df$GPP_daily_sd,na.rm=TRUE))          # simulate GPP observation error as we would have in Appling data release
    ER_pred <- rep(NA,length(t))
    ER_pred[1] <- ER[1]  
    ER_pred_sd <- rnorm(length(t),mean(df$ER_daily_sd,na.rm=TRUE),sd(df$ER_daily_sd,na.rm=TRUE))          # simulate ER observation error as we would have in Appling data release
    
    # Model log-biomass:
    for(i in 2:length(t)){
      
      # process model:
      B[i] <- B[i-1] + log(1+(mup*light[i])-mur) + rnorm(1,mean = 0,sd = sig_proc)
      
      # observation models:

      # truncated normal distribution:
      GPP_pred[i] <- MCMCglmm::rtnorm(1, mean = mup * light[i] * exp(B[i]), sd = sig_obs_GPP, lower = 0)
      ER_pred[i] <- MCMCglmm::rtnorm(1, mean = mur * exp(B[i]), sd = sig_obs_ER, lower = 0)
      
    }
    
    # save simulations for j storms:
    out_data <- data.frame(storm_id = rep(j,length(t)), time = t,light=light,B = B[c(1:length(t))], 
                           GPPsim = GPP_pred,ERsim = ER_pred)
    sim_dat <- rbind(sim_dat,out_data)  
  }
  
  # Save parameter values and simulated data:
  pars <- data.frame(no_storms = storms_n,
                     sigma_obs_GPP = sig_obs_GPP,
                     sigma_obs_ER = sig_obs_ER,
                     sigma_proc = sig_proc,
                     mup = mup,
                     mur = mur)
  out <- list(parameters = pars,data = sim_dat)
  
  return(out)
}



##==========================================================##
##   Unlogged biomass - no density-dependence on biomass    ##
##==========================================================##

simGPP_mod2 <- function(storms_n,df,mup,mur,K,sig_obs_GPP,sig_obs_ER,sig_proc){
  
# Create data frame to hold simulated data across j number of storms (defined by no_storms above):
sim_dat <- data.frame(storm_id = numeric(),
                      time = integer(),
                      light=numeric(),
                      B = numeric(),
                      GPP.pred = numeric(),
                      ER.pred = numeric())
  
# Simulate GPP given above parameters:
for(j in 1:storms_n){
    
  # Initialize data:
  dT <- 1                                         # time step
  t <- seq(from=1,to=dim(df)[1],by=dT)            # for now, all simulated storms are the same length
  light <- df$light_rel                           # relativized light 
  GPP <- df$GPP_daily_mean
  GPP_sd <- df$GPP_daily_sd
    
  # Initialize vectors for model output of B, GPP, and ER:
  B <- rep(NA,length(t))
  B[1] <- GPP[1]/light[1]
  GPP_pred <- rep(NA,length(t))
  GPP_pred[1] <- mup * light[1] * B[1]  # approx using t = 1
  ER_pred <- rep(NA,length(t))
  ER_pred[1] <- mur * B[1]  # approx using t = 1
    
  # Model log-biomass:
  for(i in 2:length(t)){
    
    # process model:
    B[i] <- B[i-1] + B[i-1] * mup * light[i] - mur * B[i-1] + rnorm(1,mean = 0,sd = sig_proc)
    
    # observation models:
    GPP_pred[i] <- mup * light[i] * B[i] + rnorm(1,mean=0,sd = sig_obs_GPP)
    ER_pred[i] <- mur * B[i] + rnorm(1,mean=0,sd = sig_obs_ER) 
  }
    
  # save simulations for j storms:
  out_data <- data.frame(storm_id = rep(j,length(t)), time = t,light=light,B = B, 
                         GPPsim = GPP_pred,ERsim = ER_pred)
  sim_dat <- rbind(sim_dat,out_data)  
  }
  
# Save parameter values and simulated data:
pars <- data.frame(no_storms = storms_n,
                   sigma_obs_GPP = sig_obs_GPP,
                   sigma_obs_ER = sig_obs_ER,
                   sigma_proc = sig_proc,
                   mup = mup,
                   mur = mur,
                   K = K)
out <- list(parameters = pars,data = sim_dat)

return(out)
}
