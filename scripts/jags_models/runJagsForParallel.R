## This script will be used to run all of the 
## bayesian cp models for the pfactor-nonlinear mauscript

## load library(s)
library(doParallel)

## Declare the model
model_string <- "model {

  # Priors for population-level effects
  cp_0 = 0  # mcp helper value.
  cp_2 = ABC  # mcp helper value.

  cp_1 ~ dunif(0,30)
  int_1 ~ dnorm(0, 1/(10)^2) 
  int_2 ~ dnorm(0, 1/(10)^2) 
  x_1 ~ dnorm(0, 1/(10)^2) 
  x_2 ~ dnorm(0, 1/(10)^2) 
  age_beta ~ dnorm(0, 1/(10)^2)
  # Hyperprior for group-level random intercepts
  tau_u <- pow(sigma_u, -2)
  sigma_u ~ dunif(0, 10)
  # Hyperprior for group-level random intercepts
  tau_u2 <- pow(sigma_u2, -2)
  sigma_u2 ~ dunif(0, 10)

  # Model and likelihood
  for (i_ in 1:length(x)) {
    X_1_[i_] = min(x[i_], cp_1)
    X_2_[i_] = min(x[i_], cp_2) - cp_1
    
    # Fitted value
    y_[i_] = 
    
      # Segment 1: y ~ 1 + x
      (x[i_] >= cp_0) * (x[i_] < cp_1) * int_1 + 
      (x[i_] >= cp_0) * (x[i_] < cp_1) * u[group[i_]] +
      (x[i_] >= cp_0) * (x[i_] < cp_1) * u2[group2[i_]] +
      (x[i_] >= cp_0) * (x[i_] < cp_1) * age_beta * age[i_] +
      (x[i_] >= cp_0) * (x[i_] < cp_1) * x_1 * X_1_[i_] + 
    
      # Segment 2: y ~ 1 ~ 1 + x
      (x[i_] >= cp_1) * int_2 +
      (x[i_] >= cp_1) * u[group[i_]] +
      (x[i_] >= cp_1) * u2[group2[i_]] +
      (x[i_] >= cp_1) * age_beta * age[i_]+
      (x[i_] >= cp_1) * x_2 * X_2_[i_]

    # Likelihood and log-density for family = poisson()
    y[i_] ~ dpois(exp(y_[i_]))
    loglik_[i_] = logdensity.pois(y[i_], exp(y_[i_]))
  }
  for (j in 1:G) {
    u[j] ~ dnorm(0, tau_u)
  }
  for (j in 1:G2) {
    u2[j] ~ dnorm(0, tau_u2)
  }
}
"

# Now run through all of these model options in a for parallel loop
mod.dv <- c("cbcl_scr_syn_internal_r", "cbcl_scr_syn_external_r")
mod.iv <- c("cbcl_scr_syn_thought_r", "cbcl_scr_syn_attention_r", "cbcl_scr_syn_internal_r")
all.mods <- expand.grid(mod.dv, mod.iv)
all.mods$Var1 <- as.character(all.mods$Var1)
all.mods$Var2 <- as.character(all.mods$Var2)
all.mods <- all.mods[-which(as.character(all.mods$Var1)==as.character(all.mods$Var2)),]
#all.mods <- all.mods[6,]
index.val <- 1
all.dat <- list()
for(z in 1:nrow(all.mods)){
  ## Now go through the wave value
  for(i in 1:4){
    ## create the data
    if(i == 1){
      in.dat <- read.csv("./data/ABCD_CBCL_BL.csv")
    }
    if(i == 2){
      in.dat <- read.csv("./data/ABCD_CBCL_Y1.csv")
    }
    if(i == 3){
      in.dat <- read.csv("./data/ABCD_CBCL_Y2.csv")
    }
    if(i == 4){
      in.dat <- read.csv("./data/ABCD_CBCL_Y3.csv")
    }
    if(i > 1){
      in.dat.bp <- read.csv("./data/ABCD_CBCL_BL.csv")
      in.dat.bp.family <- in.dat.bp[,c("src_subject_id", "site_id_l", "rel_family_id")]
      in.dat <- merge(in.dat, in.dat.bp.family, by=c("src_subject_id", "site_id_l"), suffixes = c(".x",""))
    }
    tmp.dat <- in.dat
    tmp.dat <- in.dat[complete.cases(tmp.dat[,c(all.mods$Var1[z], all.mods$Var2[z], "interview_age")]),]
    data_jags <- list(
      x = tmp.dat[,all.mods$Var2[z]],
      y = tmp.dat[,all.mods$Var1[z]],
      age = tmp.dat$interview_age,
      group = as.numeric(factor(tmp.dat$site_id_l)),
      G = length(unique(tmp.dat$site_id_l)),
      group2 = as.numeric(factor(tmp.dat$rel_family_id)),
      G2 = length(unique(tmp.dat$rel_family_id))
    )
    all.dat[[index.val]] <- data_jags
    index.val <- index.val + 1
  }
}

## Now run through all of the models
# set up a cluster called 'cl'
ncores =6
cl = makeCluster(ncores)
# register the cluster
registerDoParallel(cl)
## Estimate models here
#model_jagsout = foreach(i=1:length(all.dat), .packages = c("rjags", "coda")) %dopar% {
model_jagsout = foreach(i=17, .packages = c("rjags", "coda")) %dopar% {
  ## Quick mod to the model statement here
  max.val <- max(all.dat[[i]]$x)
  model_stringT <- gsub(pattern = "ABC", replacement = max.val, x = model_string)
  ## init mod
  result_case = jags.model(textConnection(model_stringT), data = all.dat[[i]], n.chains = 1, n.adapt = 1000)
  ## estimate mod
  result_case = coda.samples(result_case, variable.names = c("cp_1", "int_1", "int_2", "x_1", "x_2", "age_beta","sigma_u", "sigma_u2"), n.chains = 1, n.iter = 1000, thin = 5)
  result_case
  
}
stopCluster(cl)
