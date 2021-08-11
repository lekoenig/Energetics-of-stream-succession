data{
  int <lower=1> N;  // number of time steps (days)
  vector[N] light;  // relativized to max value
  vector[N] GPP;   
  vector[N] GPP_sd; // sd on daily GPP estimates, from Appling et al. 2018

}

transformed data{
}

parameters{
  
  //biomass growth + loss parameters
  real alpha;    // intrinsic rate of increase
  real phi;    // degree of density dependence
  real GPP_mod[N];  // latent biomass

  //error parameters
  real<lower=0> sigma_proc;  //process error
  real<lower=0> sigma_GPP;   // obs error for GPP
}

transformed parameters{
}
    
model{

  // Initialize biomass + GPP:
    GPP_mod[1] ~ normal(log(GPP[1]),1);

 // Process model:
   for(i in 2:N){
      GPP_mod[i] ~ normal(phi*GPP_mod[i-1] + alpha*light[i], sigma_proc);
   }

  // GPP observation model:
  for(i in 2:N){
      GPP[i] ~ normal(exp(GPP_mod[i]),sigma_GPP)T[0,]; 
  }
  
  // Priors on model parameters:
  alpha ~ normal(0,1);           // prior on photosynthesis rate mup
  phi ~ normal(0,1);         // prior on respiration rate mur
  sigma_proc ~ normal(0,1);    // prior on process error
  sigma_GPP ~ normal(mean(GPP_sd),sd(GPP_sd))T[0,];     // observation error on GPP
  }
    
generated quantities{
 
 real GPP_tilde[N];              // posterior predictive check on GPP

  for(i in 1:N){
    GPP_tilde[i] = normal_rng(exp(GPP_mod[i]), sigma_GPP);

  }
}

