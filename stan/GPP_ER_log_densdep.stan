data{
  int <lower=1> N;  // number of time steps (days)
  vector[N] light;  // relativized to max value
  vector[N] GPP;    
  vector[N] ER;     
}

transformed data{
}

parameters{
  
  //biomass growth + loss parameters
  real<lower = 0> mup;    // photosynthesis rate
  real<lower = 0> mur;    // respiration rate
  real K_raw;            // carrying capacity
  real<lower = 0> B[N];  // latent biomass

  //error parameters
  real<lower=0> sigma_proc;  //process error
  real<lower=0> sigma_GPP;   // obs error for GPP
  real<lower=0> sigma_ER;   // obs error for ER

}

transformed parameters{
  real K;
  vector[N] GPPmod;
  vector[N] ERmod;

  //implies: K ~ normal(150,25):
  K = 150 + 25 * K_raw;
  
  for(i in 1:N){
    GPPmod[i] = mup * light[i] * exp(B[i]) * (1-exp(B[i])/K);
    ERmod[i] = mur * exp(B[i]);
  }
  

  
}
    
model{

  // Initialize biomass, GPP, and ER:

    //B[1] ~ normal(2.3,1);   
    B[1] ~ normal(log(GPP[1]/(light[1])),1);

 // Process model:
   for(i in 2:N){
     
       B[i] ~ normal(B[i-1] + log(1+(mup*light[i-1])*(1-exp(B[i-1])/K)-mur), sigma_proc);

   }

  // GPP observation model:
  for(i in 2:N){
      GPP[i] ~ normal(GPPmod[i],sigma_GPP); 
      ER[i] ~ normal(ERmod[i],sigma_ER);
  }
  
  // Priors on model parameters:
  mup ~ normal(0,1);           // prior on photosynthesis rate mup
  mur ~ normal(0,1);         // prior on respiration rate mur
  K_raw ~ normal(0,1);
  sigma_proc ~ normal(0,1);    // prior on process error
  sigma_GPP ~ normal(0,1);     // observation error on GPP
  sigma_ER ~ normal(0,1);      // observation error on ER
  }
    
generated quantities{
 
 real GPP_tilde[N];              // posterior predictive check on GPP
 real ER_tilde[N];
  
  for(i in 1:N){
    GPP_tilde[i] = normal_rng(mup * light[i] * exp(B[i]) * (1-exp(B[i])/K),sigma_GPP);
    ER_tilde[i] = normal_rng(mur * exp(B[i]),sigma_ER);

  }
 
 
}

