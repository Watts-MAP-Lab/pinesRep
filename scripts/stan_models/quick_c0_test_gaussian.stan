data {
  int<lower=0> N;
  vector[N] x;
  array[N] int<lower=0> y;
  int<lower=0> count_group;       // Number of groups in N_group
  int<lower=0> count_group2;      // Number of groups in N_group2
  vector[N] x2;                   // age predictor
  int N_group[N];                 // group membership vector 1
  int N_group2[N];                // group membership vector 2

}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
    vector[1] alpha; // intercepts
    vector[2] beta; // slopes
    vector[count_group] aa1; // random intercept site
    vector[count_group2] aa2; // random intercept family
    real<lower=0> sigma_p; // sd for intercept global
    real<lower=0> sigma_p2; // sd for intercept global
    real<lower=0> phi; //neg binomvar
    real<lower=0, upper=max(x)> r; // change point with implicit uniform prior
}

// The model to be estimated.
model {
  alpha[1] ~ normal(0,1);
  beta[1] ~ normal(0,1);
  beta[2] ~ normal(0,1);
  phi ~ uniform(0, 500);
  vector[N] mu;
  for (n in 1:N) {
    mu[n] = alpha[1] + aa1[N_group[n]]+aa2[N_group2[n]] + beta[1] * x[n] + beta[3] * x2[n];
  }
  aa1 ~ normal(0,sigma_p);
  aa2 ~ normal(0,sigma_p2);
  //y ~ poisson_log(mu);
  //y ~ neg_binomial_2_log(mu, phi);
  y ~ normal(mu, phi);

}

generated quantities {
  vector[N] log_lik;
  vector[N] mu;
  for (n in 1:N) { 
    mu[n] = alpha[1] + aa1[N_group[n]]+aa2[N_group2[n]] + beta[1] * x[n] + beta[2] * x2[n];
  }
  for (n in 1:N) { 
    log_lik[n] = normal_lpdf(y[n] | mu[n], phi);
  }
}