
  
  data{
  int <lower=1> N;             // ~total~ number of time steps (days)
  int <lower=1> storm_N;       // number of storms
  int <lower=1> storm_ID[N];   // storm ID
  
  vector[N] light;             // relativized to max value
  vector[N] GPP;               // for now, simulated GPP values
  }
  
  parameters{
  
  //biomass growth parameters and hyperparameters:
  vector <lower=0> [storm_N] r;     // vector of storm r's: specific growth rate
  vector [storm_N] b;               // vector of storm b's: r/K term in Ricker model
  real B[N];                        // log-transformed latent biomass
  real<lower=0> mu_r;               // mean of the population model (mean r across storms)
  real<lower=0> sd_r;               // variance of the population model (variance in r across storms)
  real mu_b;                        // mean of the population model (mean b across storms)
  real<lower=0> sd_b;               // variance of the pouplation model (variance in b across storms)


  //error parameters:
  real<lower=0> sigma_proc; //process error
  real<lower=0> sigma_obs; // obs error for gpp
  }


  model{
  
  // Initialize biomass:
  vector[N] GPPmod;
  B[1]~normal(log(GPP[1]/light[1]),0.05);
  GPPmod[1] = light[1] * exp(B[1]);
  
  // Process model:
   
   // Population model:
   r ~ normal(mu_r,sd_r);
   b ~ normal(mu_b,sd_b);
   
   for(i in 2:N){
   
   if(storm_ID[i] > storm_ID[i-1])
    {
    B[i] ~ normal(log(GPP[i]/light[i]),0.05);
    }
   else {
    B[i] ~ normal(B[i-1] + r[storm_ID[i]] + b[storm_ID[i]]*exp(B[i-1]), sigma_proc);
   }
    GPPmod[i] = light[i] * exp(B[i]);

   }
   
   //for(j in 1:storm_N){
   
   //r[j] ~ normal(mu_r,sd_r);
   //b[j] ~ normal(mu_b,sd_b);
   
   //}
   
  // GPP observation model:
  for(i in 2:N){
      GPP[i] ~ normal(GPPmod[i],sigma_obs);   //observation model
  }
  
  // Priors:
  r ~ normal(0,1);             // prior on growth rate, r
  b ~ normal(0,1);             // prior on ricker model term r/k
  sigma_proc ~ normal(0,1);    // prior on process error
  sigma_obs ~ normal(0,0.5);     // prior on observation error
  
  
  // Hyper-priors: 
  mu_r ~ normal(0,1);
  sd_r ~ normal(0,0.25);
  mu_b ~ normal(0,1);
  sd_b ~ normal(0,0.25);
  }

  generated quantities{
  real GPP_tilde[N];              // posterior predictive check on GPP
  
  for(i in 1:N){
    GPP_tilde[i] = normal_rng(light[i] * exp(B[i]),sigma_obs);
  }
  }
    
  
