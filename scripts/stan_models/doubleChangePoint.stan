data {
  int<lower=0> N;                 // Number of observations
  int<lower=0> count_group;       // Number of groups in N_group
  int<lower=0> count_group2;      // Number of groups in N_group2
  vector[N] x;                    // primary predictor
  vector[N] x2;                   // age predictor
  int N_group[N];                 // group membership vector 1
  int N_group2[N];                // group membership vector 2
  int<lower=0> y[N];              // Count outcome
  real beta;                      // fixed gain of logistic changepoint
  
}

parameters {
  real a1;
  real a2;
  real a3;
  real b1;
  real b2;
  real b3;
  real b4;
  vector[count_group] aa1;
  vector[count_group2] aa2;
  real<lower=0> sigma_p; // sd for intercept global
  real<lower=0> sigma_p2; // sd for intercept global
  ordered[2] alpha; //  cp locations for ordered logit
  real<lower=0> phi; //neg binom added var
  
}

transformed parameters  {
  vector[N] w1;
  vector[N] w2;
  vector[N] w3;
  w1 = inv_logit(beta + alpha[1]*x);
  w2 = inv_logit(beta + alpha[2]*x) - inv_logit(beta + alpha[1]*x);
  w3 = 1 - inv_logit(beta + alpha[2]*x);
  vector[N] lp;
  for (j in 1:N) {
    lp[j] = w1[j]*(a1+aa1[N_group[j]]+aa2[N_group2[j]]+b1*x[j]+b3*x2[j])+
      w2[j]*(a2+aa1[N_group[j]]+aa2[N_group2[j]]+b2*x[j]+b3*x2[j])+
      w3[j]*(a3+aa1[N_group[j]]+aa2[N_group2[j]]+b4*x[j]+b3*x2[j]);
  }
}

model {
  // Prior part of Bayesian inference
  a1 ~ normal(0,1);
  a2 ~ normal(0,1);
  a3 ~ normal(0,1);
  b1 ~ normal(0,1);
  b2 ~ normal(0,1);
  b3 ~ normal(0,1);
  b4 ~ normal(0,1);
  sigma_p ~ exponential(12);
  sigma_p2 ~ exponential(4);
  alpha[1] ~ normal(30,3);
  alpha[2] ~ normal(60,3);
  
  phi ~ uniform(0, 500);
  
  // Likelihood part of Bayesian inference
  aa1 ~ normal(0,sigma_p);
  aa2 ~ normal(0,sigma_p2);
  //y ~ poisson_log(lp);
  y ~ neg_binomial_2_log(lp, phi);
}
