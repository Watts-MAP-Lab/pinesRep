## This script will be used to run all of the 
## Bayesian cp models for the pfactor-nonlinear manuscript

## load library(s)
library(rstan)
library(bayesplot)
library(loo)

# Now run through all of these model options in a for parallel loop
mod.dv <- c("cbcl_scr_syn_internal_r", "cbcl_scr_syn_external_r", "cbcl_scr_syn_attention_r", "cbcl_scr_syn_thought_r")
#mod.dv <- c("cbcl_scr_syn_internal_r", "cbcl_scr_syn_external_r")#, "cbcl_scr_syn_attention_r", "cbcl_scr_syn_thought_r")

mod.iv <- mod.dv
iter <- 1:4
all.mods <- expand.grid(mod.dv, mod.iv, iter)
all.mods$Var1 <- as.character(all.mods$Var1)
all.mods$Var2 <- as.character(all.mods$Var2)
all.mods <- all.mods[-which(as.character(all.mods$Var1)==as.character(all.mods$Var2)),]
all.mods <- unique(all.mods)
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
  beta=-8
)

## Now make the scaled values here
all.dat <- data_jags
file.out <- paste("./data/brmsModsOut/model_InvLogit_rawX_GAUS_2CP_allmods_", rowID, ".RDS", sep='')
stanmonitor = c("alpha", "beta", "phi", "r", "sigma_p", "sigma_p2", "log_lik", "mu")
stanmonitor <- c("a1","a2","b1","b2","b3","alpha", "phi", "mu", "log_lik")
stanmonitor <- c("a1","a2","a3","b1","b2","b3","b4","alpha", "phi", "mu", "log_lik")
if(!file.exists(file.out)){
  result_case = stan(file="./scripts/stan_models/doubleChangePoint_gaus.stan", 
                     data = all.dat, cores=3,chains=3, refresh = 1000, 
                     pars = stanmonitor, 
                     iter=15000, warmup = 5000, thin = 2,control = list(max_treedepth=10))
  #saveRDS(result_case, file.out)
}else{
  print("Done")
}

summary.vals <- rstan::summary(result_case)$summary
pred.vals <- summary.vals[grep(x = rownames(summary.vals), pattern = "mu"),]
out.cor <- cor(pred.vals[,"mean"], data_jags$y, method="s")
log_lik6 <- extract_log_lik(result_case)
out.looic <- loo::loo(log_lik6, moment_match = TRUE)
## Now create all of the figure values
# iter.vals <- c("alpha[1]", "alpha[2]", "beta[1]", "beta[2]", "beta[3]", "r[1]", "r[2]")
# file.out2 <- paste("./data/outPlot/tracePlot_GAUS_2CP_", rowID, ".pdf", sep='')
# pdf(file.out2)
# for(i in iter.vals){
#   print(bayesplot::mcmc_trace(result_case, i))
# }
# dev.off()
out.list <- list(mod.cor = out.cor,
                 out.sum = summary.vals,
                 out.looic = out.looic)
saveRDS(out.list, file.out)