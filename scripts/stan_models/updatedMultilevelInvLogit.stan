data {
  // Define variables in data
  // Number of observations (an integer)
  int<lower=0> N;
  int<lower=0> count_group;
  //int<lower=0> count_group2;
  
  // primary predictor
  vector[N] x;
  // age predictor
  vector[N] x2;
  // group membership vector 1
  int N_group[N];
  // group membership vector 2
  //int N_group2[N];
  // Count outcome
  array[N] int<lower=0> y;
  real beta; //fixed gain of logistic changepoint
}

parameters {
  real a1;
  real a2;
  real b1;
  real b2;
  real b3;
  vector[count_group] aa1;
  //vector[count_group2] aa2;
  real<lower=0> sigma_p; // sd for intercept global
  //real<lower=0> sigma_p2; // sd for intercept global
  real<lower=-1*beta, upper =-51*beta>  alpha;

}

transformed parameters  {
  vector[N] w = inv_logit(alpha + beta*x);
}

model {
  // Prior part of Bayesian inference
  a1 ~ normal(0,1);
  a2 ~ normal(0,1);
  b1 ~ normal(0,1);
  b2 ~ normal(0,1);
  b3 ~ normal(0,1);
  aa1 ~ normal(0,sigma_p);
  sigma_p ~ exponential(12);
  //aa2 ~ normal(0,sigma_p2);
  //sigma_p2 ~ exponential(1);
  alpha ~ normal(50,3);
  
  vector[N] lp;
  for (j in 1:N) {
    lp[j] = w[j]*(a1+aa1[N_group[j]]+b1*x[j]+b3*x2[j])+(1-w[j])*(a2+aa1[N_group[j]]+b2*x[j]+b3*x2[j]);
    //+aa1[N_group[j]]
  }

  // Likelihood part of Bayesian inference
  y ~ poisson_log(lp);
}
