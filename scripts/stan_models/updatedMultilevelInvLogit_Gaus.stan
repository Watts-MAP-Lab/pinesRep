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
  real b1;
  real b2;
  real b3;
  vector[count_group] aa1;
  vector[count_group2] aa2;
  real<lower=0> sigma_p; // sd for intercept global
  real<lower=0> sigma_p2; // sd for intercept global
  real<lower=1*beta, upper =-51*beta>  alpha;
  real<lower=0> phi; //neg binom added var

}

transformed parameters  {
  vector[N] w = inv_logit(alpha + beta*x);
  vector[N] lp;
  for (j in 1:N) {
    lp[j] = w[j]*(a1+aa1[N_group[j]]+aa2[N_group2[j]]+b1*x[j]+b3*x2[j])+(1-w[j])*(a2+aa1[N_group[j]]+aa2[N_group2[j]]+b2*x[j]+b3*x2[j]);
    //+aa1[N_group[j]]
  }
}

model {
  // Prior part of Bayesian inference
  a1 ~ normal(0,1);
  a2 ~ normal(0,1);
  b1 ~ normal(0,1);
  b2 ~ normal(0,1);
  b3 ~ normal(0,1);
  sigma_p ~ exponential(12);
  sigma_p2 ~ exponential(4);
  phi ~ uniform(0, 500);
  
  // Likelihood part of Bayesian inference
  aa1 ~ normal(0,sigma_p);
  aa2 ~ normal(0,sigma_p2);
  //y ~ poisson_log(lp);
  //y ~ neg_binomial_2_log(lp, phi);
  y ~ normal(lp,phi);
}

generated quantities {
  vector[N] log_lik;
  vector[N] mu;
  for (j in 1:N) { 
    mu[j] = w[j]*(a1+aa1[N_group[j]]+aa2[N_group2[j]]+b1*x[j]+b3*x2[j])+(1-w[j])*(a2+aa1[N_group[j]]+aa2[N_group2[j]]+b2*x[j]+b3*x2[j]);
  }
  for (n in 1:N) { 
    log_lik[n] = neg_binomial_2_log_lpmf(y[n] | mu[n], phi);
  }
}