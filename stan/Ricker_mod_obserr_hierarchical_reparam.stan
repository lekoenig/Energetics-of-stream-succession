data{
  int <lower=1> N;           // number of time steps (days)
  int <lower=1> storm_N;     // number of storms
  int <lower=1> storm_ID[N]; //storm ID
  vector[N] light;           // relativized to max value
  vector[N] GPP;             // for now, simulated GPP values
  //vector[N] GPP_sd;        // sd on daily GPP estimates, from Appling et al. 2018
}

transformed data{
}

parameters{
  //biomass growth parameters and hyperparameters:
  real B[N];                  // log-transformed latent biomass
  vector[storm_N] r_raw;         // vector of storm r's (reparameterized below)
  vector[storm_N] b_raw;         // vector of storm b's (reparameterized below)
  
  real<lower=0> mu_r;        // hyperparameters
  real<upper=0> mu_b;        // hyperparameters
  real<lower=0> sd_r;        // hyperparameters
  real<lower=0> sd_b;        // hyperparameters
  
  //error parameters
  real<lower=0> sigma_proc;  //process error
  real<lower=0> sigma_obs;   // obs error for gpp
}

transformed parameters{
  vector[N] GPPmod;
  vector[storm_N] r;
  vector[storm_N] b;
  //implies: r ~ normal(mu_r,sd_r):
  r = mu_r + sd_r * r_raw;
  //implies: b ~ normal(mu_b,sd_b):
  b = mu_b + sd_b * b_raw;

  for(i in 1:N){
    GPPmod[i] = light[i] * exp(B[i]);
  }
}

model{
  
  // Initialize biomass:
  if(GPP[1] < 0)
    B[1] ~ normal(log(0.1/light[1]),1);
  else
    B[1] ~ normal(log(GPP[1]/light[1]),1);   
  
 // Process model:
 
 // population model:

   for(i in 2:N){
     
     if(storm_ID[i] > storm_ID[i-1]){
        
        if(GPP[i] < 0){
          B[i] ~ normal(log(0.1/light[i]),1);
        }
        else {
          B[i]~normal(log(GPP[1]/light[i]),1);   
      }}
     else {
       B[i] ~ normal(B[i-1] + r[storm_ID[i]] + b[storm_ID[i]]*exp(B[i-1]), sigma_proc);
   }}

  // GPP observation model:
  for(i in 2:N){
      GPP[i] ~ normal(GPPmod[i],sigma_obs); 
  }
  
  // Priors on model parameters:
  r_raw ~ normal(0,1);              // prior on growth rate, r
  b_raw ~ normal(0,1);              // prior on ricker model term r/k
  sigma_proc ~ normal(0,1)T[0,];     // prior on process error
  sigma_obs ~ normal(0,1)T[0,];
  
  //Hyper-priors:
  mu_r ~ normal(0,1)T[0,];
  sd_r ~ normal(0,1)T[0,];
  mu_b ~ normal(0,1)T[,0];
  sd_b ~ normal(0,1)T[0,];
  
  }
    
generated quantities{
  real GPP_tilde[N];              // posterior predictive check on GPP
  
  for(i in 1:N){
    GPP_tilde[i] = normal_rng(light[i] * exp(B[i]),sigma_obs);
  }
}

