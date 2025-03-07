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
  real a1; // intercept model 1
  real a2; // intercept model 2
  real a3; // intercept model 3
  real b1; // slope model 1
  real b2; // slope model 2
  real b3; // slope **all models** age 
  real b4; // slope model 3
  vector[count_group] aa1; // random intercept for site
  vector[count_group2] aa2; // random intercept for family
  real<lower=0> sigma_p; // sd for intercept global
  real<lower=0> sigma_p2; // sd for intercept global
  positive_ordered[2] alpha; //  cp locations for ordered logit
  real<lower=0> phi; //neg binom added var
  simplex[2] pi;
	simplex[3] theta;
}

transformed data{
    real<lower=0> a;
    real<lower=0> b;
    a <- min(x);
    b <- max(x);
}


transformed parameters  {
  vector[N] w = inv_logit(alpha[1] + beta*x);
  vector[N] w2 = inv_logit(alpha[2] + beta*x);
  vector[N] lp;
  for (j in 1:N) {
    lp[j] = w[j]*(a1+aa1[N_group[j]]+aa2[N_group2[j]]+b1*x[j]+b3*x2[j])+(1-w[j])*(a2+aa1[N_group[j]]+aa2[N_group2[j]]+b2*x[j]+b3*x2[j])+(1-w2[j])*(a3+aa1[N_group[j]]+aa2[N_group2[j]]+b4*x[j]+b3*x2[j]);
    //lp[j] = w[j]*(a1+b1*x[j])+(1-w[j])*(a2+b2*x[j])+(1-w2[j])*(a3+b4*x[j]);
  }
  positive_ordered[2] alpha; 
  alpha <- a + head(cumulative_sum(theta), 2) * (b - a); 
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
  phi ~ uniform(0, 500);
  
  // Likelihood part of Bayesian inference
  aa1 ~ normal(0,sigma_p);
  aa2 ~ normal(0,sigma_p2);
  //y ~ poisson_log(lp);
  //y ~ neg_binomial_2_log(lp, phi);
  y ~ normal(lp, phi);
}

generated quantities {
  vector[N] log_lik;
  vector[N] mu;
  for (j in 1:N) { 
     mu[j] = w[j]*(a1+aa1[N_group[j]]+aa2[N_group2[j]]+b1*x[j]+b3*x2[j])+(1-w[j])*(a2+aa1[N_group[j]]+aa2[N_group2[j]]+b2*x[j]+b3*x2[j])+(1-w2[j])*(a3+aa1[N_group[j]]+aa2[N_group2[j]]+b4*x[j]+b3*x2[j]);
  }
  for (n in 1:N) { 
    log_lik[n] = neg_binomial_2_log_lpmf(y[n] | mu[n], phi);
  }
}