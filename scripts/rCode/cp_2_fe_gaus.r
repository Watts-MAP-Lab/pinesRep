in.dat <- read.csv("./data/ABCD_CBCL_BL.csv")

## Library(s)
library(brms)
library(bayesplot)
library(rstan)
#library(mcp)
library(ggplot2)
library(ggpubr)

## Now do the 2 cp model
tmp.dat <- in.dat[complete.cases(in.dat),]
data_jags <- list(
  N = nrow(tmp.dat),
  x = tmp.dat[,"cbcl_scr_syn_external_r"],
  y = tmp.dat[,"cbcl_scr_syn_internal_r"],
  x2 = tmp.dat$interview_age,
  # N_group = as.numeric(factor(tmp.dat$site_id_l)),
  # count_group = length(unique(as.numeric(factor(tmp.dat$site_id_l)))),
  # N_group2 = as.numeric(factor(tmp.dat$rel_family_id)),
  # count_group2 = length(unique(as.numeric(factor(tmp.dat$rel_family_id)))),
  beta = -6
)
stanmonitor <- c("a1","a2","a3","b1","b2","b3","b4","alpha","phi")

# result_case_DE = stan(file="./scripts/stan_models/doubleChangePointFE.stan",
#                    data = data_jags, cores=2,chains=2,
#                    pars = stanmonitor,
#                    iter=50000, warmup = 25000, thin = 5,control = list(max_treedepth=12))
# 
# result_case_SI = stan(file="./scripts/stan_models/updatedMultilevelInvLogit.stan",
#                       data = data_jags, cores=2,chains=2,
#                       pars = stanmonitor,
#                       iter=50000, warmup = 25000, thin = 5,control = list(max_treedepth=12))

result_case_DEGAU = stan(file="./scripts/stan_models/doubleChangePointFE_gaus.stan",
                      data = data_jags, cores=2,chains=2,
                      pars = stanmonitor,
                      iter=10000, warmup = 4000, thin = 5,control = list(max_treedepth=12))
saveRDS(result_case_DEGAU, file = "./data/brmsModsOut/cp_2_FE_Gaus.RDS")