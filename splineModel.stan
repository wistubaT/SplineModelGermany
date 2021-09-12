data {
  int <lower=1> N0; // number of days for which to impute infections
  int <lower=1> B2;
  int<lower=1> N;
  int<lower=1> N2;
  int cases[N2]; // reported cases
  int deaths[N2]; // reported deaths
  vector[N2] discPi;
  vector[N2] x;
  matrix[B2,N2] B;
  vector[N2] IFR;
  int EpidemicStart;
  real pop;
  real SI[N2]; // discretized serial interval based on Nishiura 
}

transformed data {
  vector[N2] SI_rev; // SI in reverse order
  vector[N2] discPi_rev; // f in reversed order
  
  for(i in 1:N2)
    SI_rev[i] = SI[N2-i+1];
    
  for(i in 1:N2) {
    discPi_rev[i] = discPi[N2-i+1];
  }
}


parameters {
  vector[B2] alpha_hier;
  real<lower=0> y;
  real<lower=0> phi;
  real<lower=0> ifr_noise;
  real<lower=0> theta;
  real<lower=0> tau1;
}

transformed parameters {
  row_vector[B2] alpha;
  vector[N2] prediction = rep_vector(0,N2);
  vector[N2] E_deaths  = rep_vector(0,N2);
  vector[N2] Rt = rep_vector(0,N2);
  vector[N2] prediction2;
  real r0 = 3.28;
  
  {
    vector[N2] cumm_sum = rep_vector(0,N2);
    alpha[1] = alpha_hier[1];
    for(b in 2:B2){
      alpha[b] = alpha[b-1] + alpha_hier[b]*theta;
    }
      
    prediction[1:N0] = rep_vector(y,N0); // learn the number of cases in the first N0 days
    cumm_sum[2:N0] = cumulative_sum(prediction[2:N0]);
    
    Rt = to_vector(alpha * B);
    
    for(j in 1:N2){
      Rt[j] = max([Rt[j], 0]);
    }
    
    for (i in (N0+1):N2) {
      real convolution = dot_product(prediction[1:(i-1)], tail(SI_rev, i-1));
      prediction[i] = Rt[i] * convolution;
    }
    E_deaths[1]= 1e-15 * prediction[1];
    prediction2 = rep_vector(0, N2);
    prediction2 = prediction .* IFR * ifr_noise;
    for (i in 2:N2){
      E_deaths[i] = max([dot_product(prediction2[1:(i-1)], tail(discPi_rev, i-1)),1e-15]);
    }
  }
}
model {
  tau1 ~ exponential(10);
  theta ~ normal(0, 1);
  alpha_hier ~ normal(0, 1);
  y ~ exponential(1/tau1);
  phi ~ normal(0,5);
  ifr_noise ~ normal(1,0.1);
  deaths[EpidemicStart:N] ~ neg_binomial_2(E_deaths[EpidemicStart:N], phi);
}

