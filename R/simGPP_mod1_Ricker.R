## Simulate data - recovery Model 1 (Ricker GPP, ER)

##==========================================================##
##                       Logged biomass                     ##
##==========================================================##

simGPP_mod1_Ricker <- function(storms_n,df,r,K,sig_obs,sig_proc){
  
  # Define parameter values:
  b <- -r/K                                   # substitute for r/K ("b")
  
  # Create data frame to hold simulated data across j number of storms (defined by storms_n above):
  sim_dat <- data.frame(storm_id = numeric(),time = integer(),light=numeric(),B = numeric(),GPPsim = numeric(),GPPsim_sd = numeric())
  
  # Simulate GPP given above parameters:
  for(j in 1:storms_n){
    
    # Initialize data:
    dT <- 1                                        # time step
    t <- seq(from=1,to=dim(df)[1],by=dT)           # for now, all simulated storms are the same length
    light <- df$light_rel                          # relativized light (unitless)
    GPP <- df$GPP_daily_mean
    GPP_sd <- df$GPP_daily_sd
    
    # Initialize vectors for model output of B, GPP, and ER:
    B <- rep(NA,length(t))
    B[1] <- log(GPP[1]/light[1])
    GPP_pred <- rep(NA,length(t))
    GPP_pred[1] <- GPP[1]                          # starting GPP; g O2 m^-2 d^-1
    GPP_pred_sd <- rnorm(length(t),mean(df$GPP_daily_sd,na.rm=TRUE),sd(df$GPP_daily_sd,na.rm=TRUE))          # simulate GPP observation error as we would have in Appling data release

    # Model log-biomass:
    for(i in 1:length(t)){
      
      # process model:
      B[i+1] <- B[i] + r + b*exp(B[i]) + rnorm(1,mean=0,sd = sig_proc)
      
      # observation model:
      GPP_pred[i] <- light[i] * (exp(B[i])) + rnorm(1,mean=0,sd = sig_obs)
    }
    
    # save simulations for j storm:
    out_data <- data.frame(storm_id = rep(j,length(t)), time = t,light=light,B = B[c(1:length(t))], GPPsim = GPP_pred,GPPsim_sd = GPP_pred_sd)
    sim_dat <- rbind(sim_dat,out_data)  
  }
  
  # Save parameter values and simulated data:
  pars <- data.frame(storms_n = storms_n,
                     sig_obs = sig_obs,
                     sig_proc = sig_proc,
                     r = r,
                     K = K)
  out <- list(parameters = pars,data = sim_dat)
  
  return(out)
}