## Simulate data - recovery Model 3 (AR1)


##==========================================================##
##                    AR1 model of biomass                  ##
##==========================================================##

simGPP_mod3_biomass <- function(df,storms_n,a,b,sig_obs_GPP,sig_proc){
  
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
    B[1] <- log(GPP[1]/light[1]) 
    GPP_pred <- rep(NA,length(t))
    GPP_pred[1] <- GPP[1] 
    GPP_pred_sd <- rnorm(length(t),mean(df$GPP_daily_sd,na.rm=TRUE),sd(df$GPP_daily_sd,na.rm=TRUE))          # simulate GPP observation error as we would have in Appling data release
    ER_pred <- rep(NA,length(t))
    ER_pred[1] <- ER[1]  
    ER_pred_sd <- rnorm(length(t),mean(df$ER_daily_sd,na.rm=TRUE),sd(df$ER_daily_sd,na.rm=TRUE))          # simulate ER observation error as we would have in Appling data release
    
    # Model log-biomass:
    for(i in 2:length(t)){
      
      # process model:
      B[i] <- a + b*B[i-1] + rnorm(1,mean = 0, sd = sig_proc)
      
      # observation model:
      GPP_pred[i] <- light[i] * exp(B[i]) + rnorm(1,mean = 0, sd = sig_obs_GPP)
    }
    
    # save simulations for j storms:
    out_data <- data.frame(storm_id = rep(j,length(t)), time = t,light=light,B = B[c(1:length(t))], 
                           GPPsim = GPP_pred,ERsim = ER_pred)
    sim_dat <- rbind(sim_dat,out_data)  
  }                                      
  
  # Save parameter values and simulated data:
  pars <- data.frame(no_storms = storms_n,
                     sigma_obs_GPP = sig_obs_GPP,
                     sigma_proc = sig_proc,
                     a = a,
                     b = b)
  out <- list(parameters = pars,data = sim_dat)
  
  return(out)
  
  
}



##==========================================================##
##                      AR1 model of GPP                    ##
##==========================================================##

simGPP_mod3 <- function(df,storms_n,phi,alpha,sig_obs_GPP,sig_proc){
  
  # Create data frame to hold simulated data across j number of storms (defined by storms_n above):
  sim_dat <- data.frame(storm_id = numeric(),
                        time = integer(),
                        light=numeric(),
                        GPPsim = numeric(),
                        GPPsim_sd = numeric())
  
  # Simulate GPP given above parameters:
  for(j in 1:storms_n){
    
    # Initialize data:
    dT <- 1                                         # time step
    t <- seq(from=1,to=dim(df)[1],by=dT)            # for now, all simulated storms are the same length
    light <- df$light_rel                           # relativized light 
    GPP <- df$GPP_daily_mean
    GPP_sd <- df$GPP_daily_sd
  
    # Initialize vectors for model output of B, GPP, and ER:
    GPP_pred <- rep(NA,length(t))
    GPP_pred[1] <- GPP[1] 
    GPP_pred_sd <- rnorm(length(t),mean(df$GPP_daily_sd,na.rm=TRUE),sd(df$GPP_daily_sd,na.rm=TRUE))          # simulate GPP observation error as we would have in Appling data release
    GPP_l <- log(GPP_pred)
    
    # Model log-GPP:
    for(i in 2:length(t)){
      
      # process model:
      GPP_l[i] <- phi*GPP_l[i-1] + alpha*light[i] + rnorm(1,mean = 0, sd = sig_proc)
      
      # observation model:
      GPP_pred[i] <- exp(GPP_l[i]) + rnorm(1,mean = 0, sd = sig_obs_GPP)
    }
    
    # save simulations for j storms:
    out_data <- data.frame(storm_id = rep(j,length(t)), time = t,light=light, 
                           GPPsim = GPP_pred, GPPsim_sd = GPP_pred_sd)
    sim_dat <- rbind(sim_dat,out_data)  
  }                                      
  
  # Save parameter values and simulated data:
  pars <- data.frame(no_storms = storms_n,
                     sigma_obs_GPP = sig_obs_GPP,
                     sigma_proc = sig_proc,
                     phi = phi,
                     alpha = alpha)
  out <- list(parameters = pars,data = sim_dat)
  
  return(out)
  
  
}


