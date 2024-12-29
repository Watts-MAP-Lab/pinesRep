## This script will be used to run all of the 
## bayesian cp models for the pfactor-nonlinear mauscript

## load library(s)
library(doParallel)
library(rstan)

# Now run through all of these model options in a for parallel loop
mod.dv <- c("cbcl_scr_syn_internal_r", "cbcl_scr_syn_external_r", "cbcl_scr_syn_attention_r", "cbcl_scr_syn_internal_r", "cbcl_scr_syn_thought_r", "cbcl_scr_syn_anxdep_r", "cbcl_scr_syn_withdep_r", "cbcl_scr_syn_somatic_r", "cbcl_scr_syn_aggressive_r", "cbcl_scr_syn_rulebreak_r")
mod.iv <- mod.dv
iter <- 1:4
all.mods <- expand.grid(mod.dv, mod.iv, iter)
all.mods$Var1 <- as.character(all.mods$Var1)
all.mods$Var2 <- as.character(all.mods$Var2)
all.mods <- all.mods[-which(as.character(all.mods$Var1)==as.character(all.mods$Var2)),]
rowID <- as.integer(commandArgs(1))
i <- all.mods[rowID,3]
## Now go through the wave value
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
tmp.dat <- in.dat[complete.cases(tmp.dat[,c(all.mods$Var1[rowID], all.mods$Var2[rowID], "interview_age")]),]
data_jags <- list(
  N = nrow(tmp.dat),
  x = tmp.dat[,all.mods$Var2[rowID]],
  y = tmp.dat[,all.mods$Var1[rowID]],
  x2 = tmp.dat$interview_age,
  N_group = as.numeric(factor(tmp.dat$site_id_l)),
  count_group = length(unique(as.numeric(factor(tmp.dat$site_id_l)))),
  N_group2 = as.numeric(factor(tmp.dat$rel_family_id)),
  count_group2 = length(unique(as.numeric(factor(tmp.dat$rel_family_id)))),
  beta = -5
)
all.dat <- data_jags
stanmonitor <- c("a1","a2","b1","b2","b3","alpha","sigma_p", "sigma_p2")
file.out <- paste("./data/brmsModsOut/model_InvLogitRandSite_allmods_", rowID, ".RDS", sep='')
if(!file.exists(file.out)){
  result_case = stan(file="./scripts/stan_models/updatedMultilevelInvLogit.stan", 
                     data = all.dat, cores=2,chains=2,
                     pars = stanmonitor, 
                     iter=17000, warmup = 5000, thin = 5, control = list(max_treedepth=11))
  saveRDS(result_case, file.out)
}else{
  print("Done")
}

summary(do.call(rbind, 
                args = get_sampler_params(result_case, inc_warmup = FALSE)),
        digits = 2)
