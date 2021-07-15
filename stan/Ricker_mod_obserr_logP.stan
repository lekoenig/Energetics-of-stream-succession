data{
  int <lower=1> N;  // number of time steps (days)
  vector[N] light;  // relativized to max value
  vector[N] GPP;    // for now, simulated GPP values
  vector[N] GPP_sd; // sd on daily GPP estimates, from Appling et al. 2018
}

transformed data{
//vector[N] P;        // model log-GPP
//P = log(GPP);
}

parameters{
  
  //biomass growth parameters
  real B[N];              // log-transformed latent biomass
  real r;       // specific growth rate
  real b;                 // r/K term in Ricker model
  
  //error parameters
  real<lower=0> sigma_proc;  //process error
  real<lower=0> sigma_obs;   // obs error for gpp
}

transformed parameters{
  vector[N] Pmod;

  for(i in 1:N){
    Pmod[i] = log(light[i]) + B[i];
  }
}
    
model{
  
  // Initialize biomass:
    B[1] ~ normal((GPP[1] - log(light[1])),1);   
  
 // Process model:
   for(i in 2:N){
       B[i] ~ normal(B[i-1] + r + b*exp(B[i-1]), sigma_proc);
   }

  // GPP observation model:
  for(i in 2:N){
      GPP[i] ~ normal(Pmod[i],sigma_obs); 
  }
  
  // Priors on model parameters:
  r ~ normal(0,1);              // prior on growth rate, r
  b ~ normal(0,1)T[,0];              // prior on ricker model term r/k
  sigma_proc ~ normal(0,1)T[0,];     // prior on process error
  sigma_obs ~ normal(0,1)T[0,];     // strong prior on observation error that corresponds w/ abs. sd on GPP posterior
  }
    
generated quantities{
  real GPP_tilde[N];              // posterior predictive check on GPP
  
  for(i in 1:N){
    GPP_tilde[i] = normal_rng(log(light[i]) + B[i],sigma_obs);
  }
}

