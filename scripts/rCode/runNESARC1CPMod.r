## This script will be used to run all of the 
## Bayesian cp models for the pfactor-nonlinear manuscript

## WATTS message here:

# If this is tenable to pull off in a week/if these models will run in enough time:
#   I want to make a point about substance use disorder in my HiTOP talk. In the NESARC data (a large epidemiologic study), I’m interested in whether there is a nonlinear association between alcohol use disorder and externalizing. I’ve worked with these data a ton and have preexisting factor scores for externalizing. I also have an alcohol use disorder symptom count. So, only AUD is a count variable.
# There are a bunch of other variables in this data file, but here are the ones of interest:
#   AUDsxct = past-year AUD symptom count
# externalizing = externalizing factor score from a paper of mine published in 2023
# For any model, we’ll have:
#   Nested random effects for two clustering variables.
# PSU is nested in STRATUM
# Fixed effects for age (AGE) and sex (SEX.x)
# Finally, I want to look at the results for two sets of data:
#   The entire dataset, which autopopulates folks with 0s for the AUD symptom count if they didn’t drink in the past year
# Conditional on consumption in the past year (py12drk). Where py12drk = 1 (as opposed to 0)
# py12drk = whether or not the person drank 12+ drinks in the past year
# For now, I think we could check out 1cp models first. As you’ll see, the data are extremely zero inflated. The zeros for externalizing are obviously below 0 because they’re factor scores, but the skew is still nuts.  File is here: https://vanderbilt.box.com/s/479q0c3kz86zu80rlu5osztmcm2tlaj8
# What do you think? Sorry to throw this on you so last minute. If it's not possible, that's no problem. I just think HiTOP would be very interested in the SUD link, as we're actively pitching revising the model to move SUD out of externalizing.

## Looks like a total of four possible models:
## 1. 1-cp model 

## load library(s)
library(rstan)
library(bayesplot)
library(loo)
library(mgcv)
library(brms)
library(visreg)
library(mcp)
library(pscl)

mod.execute <- as.integer(commandArgs(1))


in.dat <- read.csv("~/Downloads/NESARCw1_nonlinear.csv")
tmp.dat <- in.dat
tmp.dat <- in.dat[complete.cases(tmp.dat[,c("AUDsxct", "externalizing", "AGE")]),]
tmp.dat1 <- tmp.dat
data_jags <- list(
  N = nrow(tmp.dat),
  y = tmp.dat[,"AUDsxct"],
  x = tmp.dat[,"externalizing"],
  x2 = tmp.dat$AGE,
  x3 = tmp.dat$SEX.x,
  N_group = as.numeric(factor(tmp.dat$PSU)),
  count_group = length(unique(as.numeric(factor(tmp.dat$PSU)))),
  N_group2 = as.numeric(factor(tmp.dat$STRATUM)),
  count_group2 = length(unique(as.numeric(factor(tmp.dat$STRATUM))))
)

## Now do folks who have had atleast 12 drinks in the past 12 months?
tmp.dat <- in.dat[which(in.dat$py12drk==1),]
tmp.dat <- tmp.dat[complete.cases(tmp.dat[,c("AUDsxct", "externalizing", "AGE")]),]
tmp.dat2 <- tmp.dat
data_jags2 <- list(
  N = nrow(tmp.dat),
  y = tmp.dat[,"AUDsxct"],
  x = tmp.dat[,"externalizing"],
  x2 = tmp.dat$AGE,
  N_group = as.numeric(factor(tmp.dat$PSU)),
  count_group = length(unique(as.numeric(factor(tmp.dat$PSU)))),
  N_group2 = as.numeric(factor(tmp.dat$STRATUM)),
  count_group2 = length(unique(as.numeric(factor(tmp.dat$STRATUM))))
)



## Now make the scaled values here
#all.dat <- data_jags
#file.out <- paste("./data/brmsModsOut/model_rawX_NB_1CP_allmods_", rowID, ".RDS", sep='')
stanmonitor = c("alpha", "beta", "phi", "r", "sigma_p", "sigma_p2", "log_lik","mu")
if(mod.execute == 1){
result_case =  stan(file="./scripts/stan_models/quick_cp_test_NESARC.stan", 
                  data = data_jags, cores=3,chains=3, refresh = 100, 
                  pars = stanmonitor, 
                  iter=16000, warmup = 10000, thin = 3,control = list(max_treedepth=9))
saveRDS(object = result_case, file = "./data/brmsModsOut/NESARC_analyses_ALL.RDS")
}
if(mod.execute==2){
result_case2 = stan(file="./scripts/stan_models/quick_cp_test_NESARC.stan", 
                   data = data_jags2, cores=3,chains=3, refresh = 100, 
                   pars = stanmonitor, 
                   iter=16000, warmup = 10000, thin = 3,control = list(max_treedepth=9))
saveRDS(object = result_case2, file = "./data/brmsModsOut/NESARC_analyses_12months.RDS")
}


## Now plot thee results
in.mod <- readRDS("./data/brmsModsOut/NESARC_analyses_ALL.RDS")
## Grab our coefficients
in.coef <- rstan::summary(in.mod)

## Now run the fiexed effects appraoch
model = list(
  AUDsxct ~ 1 + externalizing,  # plateau (int_1)
  ~ 1 + externalizing    # joined slope (time_2) at cp_1
)
fit1 = mcp(model, data = tmp.dat, family=poisson(), cores=2, chains = 2, adapt = 1000, iter=2000)


## Now examine the randef modles
in.mod1 <- readRDS("./data/brmsModsOut/NESARC_analyses_ALL.RDS")
in.mod2 <- readRDS("./data/brmsModsOut/NESARC_analyses_12months.RDS")

## Plot the in.mod1 results here
coef.1 <- summary(in.mod1)$summary
coef.2 <- summary(in.mod2)$summary

## Now prep the data
pred.vals.one <- coef.1["alpha[1]","mean"] + coef.1["beta[1]","mean"] * data_jags$x
pred.vals.two <- coef.1["alpha[2]","mean"] + coef.1["beta[2]","mean"] * data_jags$x
cp.val <- coef.1["r","mean"]
index.lt <- which(data_jags$x<cp.val)
tmp.dat1$predOne <- pred.vals.one
tmp.dat1$predOne[-index.lt] <- pred.vals.two[-index.lt]
tmp.dat1$group <- "LT"
tmp.dat1$group[-index.lt] <- "GT"

pred.vals.one <- coef.2["alpha[1]","mean"] + coef.2["beta[1]","mean"] * data_jags2$x
pred.vals.two <- coef.2["alpha[2]","mean"] + coef.2["beta[2]","mean"] * data_jags2$x
cp.val <- coef.2["r","mean"]
index.lt <- which(data_jags2$x<cp.val)
tmp.dat2$predOne <- pred.vals.one
tmp.dat2$predOne[-index.lt] <- pred.vals.two[-index.lt]
tmp.dat2$group <- "LT"
tmp.dat2$group[-index.lt] <- "GT"

## Now create the figures here
library(ggplot2)
library(ggpubr)
p1 <- ggplot(tmp.dat1, aes(x=externalizing, y = log(AUDsxct + 1))) +
  geom_point(alpha=0) +
  theme_bw() +
  geom_hex(aes(fill=log(..count..)), bins=15) +
  geom_line(data = tmp.dat1, aes(x=externalizing, y=predOne, group=group)) +
  xlab("External") + ylab("log(AUDsxct)") +
  theme(legend.position = "NULL") +
  coord_cartesian(ylim = c(0,5)) +
  ggtitle("Full data") +
  annotate("text", x = 1.5, y=4.5, label = paste("LT slope: ", round(coef.1["beta[1]", "mean"], 2))) +
  annotate("text", x = 1.5, y=4, label = paste("GT slope: ", round(coef.1["beta[2]", "mean"], 2)))


p2 <- ggplot(tmp.dat2, aes(x=externalizing, y = log(AUDsxct + 1))) +
  geom_point(alpha=0) +
  theme_bw() +
  geom_hex(aes(fill=log(..count..)), bins=15) +
  geom_line(data = tmp.dat2, aes(x=externalizing, y=predOne, group=group)) +
  xlab("External") + ylab("log(AUDsxct)") +
  theme(legend.position = "NULL") +
  coord_cartesian(ylim = c(0,5))+
  ggtitle("Reduced data") +
  annotate("text", x = 1.5, y=4.5, label = paste("LT slope: ", round(coef.2["beta[1]", "mean"], 2))) +
  annotate("text", x = 1.5, y=4, label = paste("GT slope: ", round(coef.2["beta[2]", "mean"], 2)))

ggarrange(p1, p2)
