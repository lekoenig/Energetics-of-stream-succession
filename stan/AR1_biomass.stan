data{
  int <lower=1> N;  // number of time steps (days)
  vector[N] light;  // relativized to max value
  vector[N] GPP;    
}

transformed data{
}

parameters{
  
  //biomass growth + loss parameters
  real a;    // intrinsic rate of increase
  real b;    // degree of density dependence
  real<lower = 0> B[N];  // latent biomass

  //error parameters
  real<lower=0> sigma_proc;  //process error
  real<lower=0> sigma_GPP;   // obs error for GPP
}

transformed parameters{
  vector[N] GPPmod;

  for(i in 1:N){
    GPPmod[i] = light[i] * exp(B[i]);
  }
}
    
model{

  // Initialize biomass + GPP:
    B[1] ~ normal(log(GPP[1]/(light[1])),1);

 // Process model:
   for(i in 2:N){
     
       B[i] ~ normal(a + b*B[i-1], sigma_proc);

   }

  // GPP observation model:
  for(i in 2:N){
      GPP[i] ~ normal(GPPmod[i],sigma_GPP); 
  }
  
  // Priors on model parameters:
  a ~ normal(0,1);           // prior on photosynthesis rate mup
  b ~ normal(0,1);         // prior on respiration rate mur
  sigma_proc ~ normal(0,1);    // prior on process error
  sigma_GPP ~ normal(0,1);     // observation error on GPP
  }
    
generated quantities{
 
 real GPP_tilde[N];              // posterior predictive check on GPP

  for(i in 1:N){
    GPP_tilde[i] = normal_rng(light[i] * exp(B[i]), sigma_GPP);

  }
 
 
}

