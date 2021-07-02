data{
  int <lower=1> N;  // number of time steps (days)
  vector[N] light;  // relativized to max value
  vector[N] GPP;    // for now, simulated GPP values
  vector[N] GPP_sd; // sd on daily GPP estimates, from Appling et al. 2018
}

transformed data{
}

parameters{
  
  //biomass growth parameters
  real B[N];                   // log-transformed latent biomass
  vector[2] beta_tilde;
  cholesky_factor_corr[2] L;
  
  //error parameters
  vector<lower=0>[2] sigma_b; // sd on non-centered parameters (I think!)
  real<lower=0> sigma_proc;  //process error
  real<lower=0> sigma_obs;   // obs error for gpp
}

transformed parameters{
  vector[N] GPPmod;
  vector[2] beta;
  
  beta = [0.5, -0.05]' + sigma_b .* (L * beta_tilde);

  for(i in 1:N){
    GPPmod[i] = light[i] * exp(B[i]);
  }
}
    
model{
  
  // Initialize biomass:
  if(GPP[1] < 0)
    B[1] ~ normal(log(0.01/light[1]),0.005);
  else
    B[1] ~ normal(log(GPP[1]/light[1]),0.005);   //small sd around initialized biomass 
  
 // Process model:
   for(i in 2:N){
       B[i] ~ normal(B[i-1] + beta[1] + beta[2]*exp(B[i-1]), sigma_proc);
   }

  // GPP observation model:
  for(i in 2:N){
      GPP[i] ~ normal(GPPmod[i],sigma_obs); 
  }
  
  // Priors on model parameters:
  beta_tilde ~ normal(0,1);
  sigma_b ~ student_t(7, 0, 1);
  L ~ lkj_corr_cholesky(2);
  sigma_proc ~ normal(0,1);                        // prior on process error
  sigma_obs ~ normal(mean(GPP_sd),sd(GPP_sd));     // strong prior on observation error that corresponds w/ abs. sd on GPP posterior
  }
    
generated quantities{
  real GPP_tilde[N];              // posterior predictive check on GPP
  
  for(i in 1:N){
    GPP_tilde[i] = normal_rng(light[i] * exp(B[i]),sigma_obs);
  }
}

