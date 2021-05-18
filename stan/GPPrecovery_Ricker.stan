
  
  data{
  int <lower=1> N; // number of time steps (days)
  vector[N] light; // relativized to max value
  vector[N] GPP;   // for now, simulated GPP values
  }
  
  transformed data{
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
  B[1]~normal(log(GPP[1]/light[1]),0.05);   //small sd around initialized biomass? 
  
 // Process model:
   for(i in 2:N){
       B[i] ~ normal(B[i-1] + r + b*exp(B[i-1]), sigma_proc);
   }

  // GPP observation model:
  for(i in 2:N){
      GPP[i] ~ normal(GPPmod[i],sigma_obs);   
  }
  
  // Priors on model parameters:
  r ~ normal(0,1);             // prior on growth rate, r
  b ~ normal(0,1);             // prior on ricker model term r/k
  sigma_proc ~ normal(0,2);    // prior on process error
  sigma_obs ~ normal(0,0.5);   // prior on observation error
  }
    
  generated quantities{
  real GPP_tilde[N];              // posterior predictive check on GPP
  
  for(i in 1:N){
    GPP_tilde[i] = normal_rng(light[i] * exp(B[i]),sigma_obs);
  }
  }
    
  
