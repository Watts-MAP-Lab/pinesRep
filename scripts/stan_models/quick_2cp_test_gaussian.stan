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

transformed data{
    real<lower=0> a;
    real<lower=0> b;
    a <- min(x);
    b <- max(x);
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
    vector[3] alpha; // intercepts
    vector[4] beta; // slopes
    vector[count_group] aa1; // random intercept site
    vector[count_group2] aa2; // random intercept family
    real<lower=0> sigma_p; // sd for intercept global
    real<lower=0> sigma_p2; // sd for intercept global
    real<lower=0> phi; //neg binomvar
    simplex[2] pi;
	  simplex[3] theta;
}
transformed parameters{
    positive_ordered[2] r; 
    r <- a + head(cumulative_sum(theta), 2) * (b - a); 
}

// The model to be estimated.
model {
  alpha[1] ~ normal(0,1);
  alpha[2] ~ normal(0,1);
  alpha[2] ~ normal(0,1);
  beta[1] ~ normal(0,1);
  beta[2] ~ normal(0,1);
  beta[3] ~ normal(0,1);
  beta[3] ~ normal(0,1);
  phi ~ uniform(0, 500);
  vector[N] mu;
  vector[N] mu1;
  vector[N] mu2;
  vector[N] mu3;
  for (n in 1:N) {
    mu1[n] = alpha[1] + aa1[N_group[n]]+aa2[N_group2[n]] + beta[1] * x[n] + beta[4] * x2[n];
    mu2[n] = alpha[2] + aa1[N_group[n]]+aa2[N_group2[n]] + beta[2] * x[n] + beta[4] * x2[n];
    mu3[n] = alpha[3] + aa1[N_group[n]]+aa2[N_group2[n]] + beta[3] * x[n] + beta[4] * x2[n];
    if (x[n] < r[2])
      mu[n] = mu2[n];
    else if (x[n] < r[1])
      mu[n] = mu1[n];
    else
      mu[n] = mu3[n];
  }
  aa1 ~ normal(0,sigma_p);
  aa2 ~ normal(0,sigma_p2);
  //y ~ poisson_log(mu);
  y ~ normal(mu, phi);
}

generated quantities {
  vector[N] log_lik;
  vector[N] mu;
  vector[N] mu1;
  vector[N] mu2;
  vector[N] mu3;
  for (n in 1:N) {
    mu1[n] = alpha[1] + aa1[N_group[n]]+aa2[N_group2[n]] + beta[1] * x[n] + beta[4] * x2[n];
    mu2[n] = alpha[2] + aa1[N_group[n]]+aa2[N_group2[n]] + beta[2] * x[n] + beta[4] * x2[n];
    mu3[n] = alpha[3] + aa1[N_group[n]]+aa2[N_group2[n]] + beta[3] * x[n] + beta[4] * x2[n];
    if (x[n] < r[2])
      mu[n] = mu2[n];
    else if (x[n] < r[1])
      mu[n] = mu1[n];
    else
      mu[n] = mu3[n];
  }
  for (n in 1:N) { 
    log_lik[n] = normal_lpdf(y[n] | mu[n], phi);
  }
}