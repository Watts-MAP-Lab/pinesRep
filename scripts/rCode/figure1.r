## Figure 1 will build the case for the need of the changepoint model
## IT will compare three different models
## One with a basic Bayesian Poisson model
## One with the 1 cp model
## One with a 2 cp model

## I will need to compare the convergence of the CP across the 1 and 2 cp models
## and then compare how the prediction varies across the three models
## I am going to do this in the overall 

## First thing is to load the data
in.dat <- read.csv("./data/ABCD_CBCL_BL.csv")

## Library(s)
library(brms)
library(bayesplot)
library(rstan)
#library(mcp)
library(ggplot2)
library(ggpubr)

## Train model with 0 cp
# mod.base <- brms::brm(formula = cbcl_scr_syn_external_r ~ cbcl_scr_syn_internal_r + interview_age + (1|site_id_l/rel_family_id), family = "Poisson", init = 10000, iter = 20000, thin=5, cores=2, chains=2,data = in.dat)
#mod.base.quick <- lme4::glmer.nb(formula = cbcl_scr_syn_external_r ~ cbcl_scr_syn_internal_r + interview_age + (1|rel_family_id/site_id_l), data = in.dat)
# saveRDS(mod.base, file = "./data/brmsModsOut/singlePois.RDS")
mod.base <- readRDS("./data/brmsModsOut/singlePois.RDS")
sum.none <- summary(mod.base)

## Load the 1 cp model 
mod.1cp <- readRDS("./data/brmsModsOut/model_rawX_NB_allmods_1.RDS")

## Now do the 2 cp model
tmp.dat <- in.dat[complete.cases(in.dat),]
data_jags <- list(
  N = nrow(tmp.dat),
  x = tmp.dat[,"cbcl_scr_syn_external_r"],
  y = tmp.dat[,"cbcl_scr_syn_internal_r"],
  x2 = tmp.dat$interview_age,
  N_group = as.numeric(factor(tmp.dat$site_id_l)),
  count_group = length(unique(as.numeric(factor(tmp.dat$site_id_l)))),
  N_group2 = as.numeric(factor(tmp.dat$rel_family_id)),
  count_group2 = length(unique(as.numeric(factor(tmp.dat$rel_family_id)))),
  beta = -6
)
# stanmonitor <- c("a1","a2","a3","b1","b2","b3","b4","alpha", "lp", "aa1", "aa2")
#
# result_case = stan(file="./scripts/stan_models/doubleChangePoint.stan",
#                    data = data_jags, cores=2,chains=2,
#                    pars = stanmonitor,
#                    iter=50000, warmup = 25000, thin = 5,control = list(max_treedepth=12))
# saveRDS(result_case, file="./data/brmsModsOut/two_cp_converge.RDS")
# q()
result_case = readRDS("./data/brmsModsOut/two_cp_converge.RDS")

## Grab all of the aa1 & aa2 values
# aa1.ints <- tmp[grep(pattern = "aa1", x = rownames(tmp)),]
# plot(data_jags$x, aa1.ints[data_jags$N_group])
# cor(data_jags$x, aa1.ints[data_jags$N_group])
# 
# aa2.ints <- tmp[grep(pattern = "aa2", x = rownames(tmp)),]
# plot(data_jags$x, aa2.ints[data_jags$N_group2])
# cor(data_jags$x, aa2.ints[data_jags$N_group2])
# 
# ## Now do these sums
# plot(data_jags$x, aa2.ints[data_jags$N_group2] +  aa1.ints[data_jags$N_group])
# cor(data_jags$x, aa2.ints[data_jags$N_group2])
## Now remove result_case from mem

## Now run the model w/o any random effects here
mod.baseNR <- brms::brm(formula = cbcl_scr_syn_external_r ~ cbcl_scr_syn_internal_r + interview_age, family = "Poisson", init = 10000, iter = 20000, thin=5, cores=2, chains=2,data = in.dat)
sum.none1 <- summary(mod.baseNR)

## Now do linear models
mod.baseLin <- brms::brm(formula = cbcl_scr_syn_external_r ~ cbcl_scr_syn_internal_r + interview_age + (1|site_id_l/rel_family_id), init = 10000, iter = 20000, thin=5, cores=2, chains=2,data = in.dat)
mod.baseLinNR <- brms::brm(formula = cbcl_scr_syn_external_r ~ cbcl_scr_syn_internal_r + interview_age, init = 10000, iter = 20000, thin=5, cores=2, chains=2,data = in.dat)
mod.baseLinNR2 <- lm(formula = cbcl_scr_syn_external_r ~ cbcl_scr_syn_internal_r + interview_age,data = in.dat)

sum.noneLin <- summary(mod.baseLin)
sum.noneLin1 <- summary(mod.baseLinNR)


## Now organize all of these model predictions in the original data frame
tmp.dat$predNone <- sum.none$fixed["Intercept","Estimate"] + tmp.dat$cbcl_scr_syn_internal_r * sum.none$fixed["cbcl_scr_syn_internal_r","Estimate"]
tmp.dat$predNoneFE <- sum.none1$fixed["Intercept","Estimate"] + tmp.dat$cbcl_scr_syn_internal_r * sum.none1$fixed["cbcl_scr_syn_internal_r","Estimate"]
## Now do exp version of this
tmp.dat$predNone_exp <- exp(sum.none$fixed["Intercept","Estimate"] + tmp.dat$cbcl_scr_syn_internal_r * sum.none$fixed["cbcl_scr_syn_internal_r","Estimate"])
tmp.dat$predNoneFE_exp <- exp(sum.none1$fixed["Intercept","Estimate"] + tmp.dat$cbcl_scr_syn_internal_r * sum.none1$fixed["cbcl_scr_syn_internal_r","Estimate"])
tmp.dat$predNone_lin <- sum.noneLin$fixed["Intercept","Estimate"] + tmp.dat$cbcl_scr_syn_internal_r * sum.noneLin$fixed["cbcl_scr_syn_internal_r","Estimate"]
tmp.dat$predNoneFE_lin <- sum.noneLin1$fixed["Intercept","Estimate"] + tmp.dat$cbcl_scr_syn_internal_r * sum.noneLin1$fixed["cbcl_scr_syn_internal_r","Estimate"]


## Now do the pred vals for the 1-cp model
sum.vals <- summary(mod.1cp)$summary
inv.logit.vals <- LaplacesDemon::invlogit(sum.vals["alpha","mean"]+ tmp.dat$cbcl_scr_syn_internal_r * -5)
## Now grab the predicted values
pred.vals.one <- sum.vals["a1","mean"] + sum.vals["b1","mean"] * tmp.dat$cbcl_scr_syn_internal_r
pred.vals.two <- sum.vals["a2","mean"] + sum.vals["b2","mean"] * tmp.dat$cbcl_scr_syn_internal_r
pred.vals.combo <- exp(pred.vals.one * (inv.logit.vals) + (1 - inv.logit.vals) * pred.vals.two)
pred.vals.lin <- pred.vals.one * (inv.logit.vals) + (1 - inv.logit.vals) * pred.vals.two
tmp.dat$predOne <- pred.vals.lin
tmp.dat$predOneExp <- exp(pred.vals.lin)

## Now repeat this process for the FE pois, FE gauss, & RE Gauss models
# Start with the fix effects poisson model 
mod.1cpFE <- readRDS("./data/brmsModsOut/cp_1_FE.RDS")
sum.vals <- summary(mod.1cpFE)$summary
inv.logit.vals <- LaplacesDemon::invlogit(sum.vals["alpha","mean"]+ tmp.dat$cbcl_scr_syn_internal_r * -6)
## Now grab the predicted values
pred.vals.one <- sum.vals["a1","mean"] + sum.vals["b1","mean"] * tmp.dat$cbcl_scr_syn_internal_r
pred.vals.two <- sum.vals["a2","mean"] + sum.vals["b2","mean"] * tmp.dat$cbcl_scr_syn_internal_r
pred.vals.combo <- exp(pred.vals.one * (inv.logit.vals) + (1 - inv.logit.vals) * pred.vals.two)
pred.vals.lin <- pred.vals.one * (inv.logit.vals) + (1 - inv.logit.vals) * pred.vals.two
tmp.dat$predOneFE <- pred.vals.lin
tmp.dat$predOneFEExp <- pred.vals.combo

## Now do the gaussian 1cp models
#FE here 
mod.1cpFEG <- readRDS("./data/brmsModsOut/cp_1_FE_Gaus.RDS")
sum.vals <- summary(mod.1cpFEG)$summary
inv.logit.vals <- LaplacesDemon::invlogit(sum.vals["alpha","mean"]+ tmp.dat$cbcl_scr_syn_internal_r * -6)
## Now grab the predicted values
pred.vals.one <- sum.vals["a1","mean"] + sum.vals["b1","mean"] * tmp.dat$cbcl_scr_syn_internal_r
pred.vals.two <- sum.vals["a2","mean"] + sum.vals["b2","mean"] * tmp.dat$cbcl_scr_syn_internal_r
pred.vals.lin <- pred.vals.one * (inv.logit.vals) + (1 - inv.logit.vals) * pred.vals.two
tmp.dat$predOneFEGaus <- pred.vals.lin

#RE here
mod.1cpREG <- readRDS("./data/brmsModsOut/cp_1_RE_Gaus.RDS")
sum.vals <- summary(mod.1cpREG)$summary
inv.logit.vals <- LaplacesDemon::invlogit(sum.vals["alpha","mean"]+ tmp.dat$cbcl_scr_syn_internal_r * -6)
## Now grab the predicted values
pred.vals.one <- sum.vals["a1","mean"] + sum.vals["b1","mean"] * tmp.dat$cbcl_scr_syn_internal_r
pred.vals.two <- sum.vals["a2","mean"] + sum.vals["b2","mean"] * tmp.dat$cbcl_scr_syn_internal_r
pred.vals.lin <- pred.vals.one * (inv.logit.vals) + (1 - inv.logit.vals) * pred.vals.two
tmp.dat$predOneREGaus <- pred.vals.lin

## Do the same for the 2-cp model
## Create a matrix for the ordered inv logit vals
mod.2 <- readRDS("./data/brmsModsOut/two_cp_converge.RDS")
sum.vals <- summary(mod.2)$summary
inv.logit <- matrix(NA, nrow = nrow(tmp.dat), ncol=3)
inv.logit[,1] <- LaplacesDemon::invlogit(sum.vals["alpha[1]","mean"]+tmp.dat$cbcl_scr_syn_internal_r * -6)
inv.logit[,2] <- 1 - LaplacesDemon::invlogit(sum.vals["alpha[1]","mean"]+tmp.dat$cbcl_scr_syn_internal_r * -6)
inv.logit[,3] <- 1 - LaplacesDemon::invlogit(sum.vals["alpha[2]","mean"]+tmp.dat$cbcl_scr_syn_internal_r * -6) 
pred.vals.one <- sum.vals["a1","mean"] + sum.vals["b1","mean"] * tmp.dat$cbcl_scr_syn_internal_r
pred.vals.two <- sum.vals["a2","mean"] + sum.vals["b2","mean"] * tmp.dat$cbcl_scr_syn_internal_r
pred.vals.thr <- sum.vals["a3","mean"] + sum.vals["b4","mean"] * tmp.dat$cbcl_scr_syn_internal_r
pred.vals.lin2 <- pred.vals.one * inv.logit[,1] + pred.vals.two * inv.logit[,2] + pred.vals.thr * inv.logit[,3] 
pred.vals.combo2 <- exp(pred.vals.lin2)
tmp.dat$predTwo <- pred.vals.lin2
tmp.dat$predTwoExp <- pred.vals.combo2

## Now do the same for the random effect pois link
mod.2FE <- readRDS("./data/brmsModsOut/cp_2_FE.RDS")
sum.vals <- summary(mod.2FE)$summary
inv.logit <- matrix(NA, nrow = nrow(tmp.dat), ncol=3)
inv.logit[,1] <- LaplacesDemon::invlogit(sum.vals["alpha[1]","mean"]+tmp.dat$cbcl_scr_syn_internal_r * -6)
inv.logit[,2] <- 1 - LaplacesDemon::invlogit(sum.vals["alpha[1]","mean"]+tmp.dat$cbcl_scr_syn_internal_r * -6)
inv.logit[,3] <- 1 - LaplacesDemon::invlogit(sum.vals["alpha[2]","mean"]+tmp.dat$cbcl_scr_syn_internal_r * -6) 
pred.vals.one <- sum.vals["a1","mean"] + sum.vals["b1","mean"] * tmp.dat$cbcl_scr_syn_internal_r
pred.vals.two <- sum.vals["a2","mean"] + sum.vals["b2","mean"] * tmp.dat$cbcl_scr_syn_internal_r
pred.vals.thr <- sum.vals["a3","mean"] + sum.vals["b4","mean"] * tmp.dat$cbcl_scr_syn_internal_r
pred.vals.lin2 <- pred.vals.one * inv.logit[,1] + pred.vals.two * inv.logit[,2] + pred.vals.thr * inv.logit[,3] 
pred.vals.combo2 <- exp(pred.vals.lin2)
tmp.dat$predTwoFE <- pred.vals.lin2
tmp.dat$predTwoFEExp <- pred.vals.combo2

## Now do the Gaussian links
mod.2FEGaus <- readRDS("./data/brmsModsOut/cp_2_FE_Gaus.RDS")
sum.vals <- summary(mod.2FEGaus)$summary
inv.logit <- matrix(NA, nrow = nrow(tmp.dat), ncol=3)
inv.logit[,1] <- LaplacesDemon::invlogit(sum.vals["alpha[1]","mean"]+tmp.dat$cbcl_scr_syn_internal_r * -6)
inv.logit[,2] <- 1 - LaplacesDemon::invlogit(sum.vals["alpha[1]","mean"]+tmp.dat$cbcl_scr_syn_internal_r * -6)
inv.logit[,3] <- 1 - LaplacesDemon::invlogit(sum.vals["alpha[2]","mean"]+tmp.dat$cbcl_scr_syn_internal_r * -6) 
pred.vals.one <- sum.vals["a1","mean"] + sum.vals["b1","mean"] * tmp.dat$cbcl_scr_syn_internal_r
pred.vals.two <- sum.vals["a2","mean"] + sum.vals["b2","mean"] * tmp.dat$cbcl_scr_syn_internal_r
pred.vals.thr <- sum.vals["a3","mean"] + sum.vals["b4","mean"] * tmp.dat$cbcl_scr_syn_internal_r
pred.vals.lin2 <- pred.vals.one * inv.logit[,1] + pred.vals.two * inv.logit[,2] + pred.vals.thr * inv.logit[,3] 
tmp.dat$predTwoFEGaus <- pred.vals.lin2

## Now do the Gaussian links
mod.2REGaus <- readRDS("./data/brmsModsOut/cp_2_RE_Gaus.RDS")
sum.vals <- summary(mod.2REGaus)$summary
inv.logit <- matrix(NA, nrow = nrow(tmp.dat), ncol=3)
inv.logit[,1] <- LaplacesDemon::invlogit(sum.vals["alpha[1]","mean"]+tmp.dat$cbcl_scr_syn_internal_r * -6)
inv.logit[,2] <- 1 - LaplacesDemon::invlogit(sum.vals["alpha[1]","mean"]+tmp.dat$cbcl_scr_syn_internal_r * -6)
inv.logit[,3] <- 1 - LaplacesDemon::invlogit(sum.vals["alpha[2]","mean"]+tmp.dat$cbcl_scr_syn_internal_r * -6) 
pred.vals.one <- sum.vals["a1","mean"] + sum.vals["b1","mean"] * tmp.dat$cbcl_scr_syn_internal_r
pred.vals.two <- sum.vals["a2","mean"] + sum.vals["b2","mean"] * tmp.dat$cbcl_scr_syn_internal_r
pred.vals.thr <- sum.vals["a3","mean"] + sum.vals["b4","mean"] * tmp.dat$cbcl_scr_syn_internal_r
pred.vals.lin2 <- pred.vals.one * inv.logit[,1] + pred.vals.two * inv.logit[,2] + pred.vals.thr * inv.logit[,3] 
tmp.dat$predTwoREGaus <- pred.vals.lin2

## Now make a corr matrix with the pred vals
for.cormat <- tmp.dat[,c("cbcl_scr_syn_internal_r", "cbcl_scr_syn_external_r")]
## Now add the cut offs to these data
for.cormat$internalThreshold <- cut(for.cormat$cbcl_scr_syn_internal_r, c(-1, 10, 14, Inf), labels = c("NonclinicalI", "At-riskI", "ClinicalI"))
for.cormat$externalThreshold <- cut(for.cormat$cbcl_scr_syn_external_r, c(-1, 12, 16, Inf), labels = c("NonclinicalE", "At-riskE", "ClinicalE"))
for.cormat$collapse <- paste(for.cormat$internalThreshold, for.cormat$externalThreshold)

## Now do the ggpairs plot
library(GGally)
cor.mat <- ggpairs(for.cormat, columns=1:2,ggplot2::aes(colour = collapse))
## Now do the correlation with the group vars

## Now do a confusion matrix
for.conf.mat <- reshape2::melt(table(for.cormat$internalThreshold, for.cormat$externalThreshold))
conf.mat <- ggplot(data =  for.conf.mat, mapping = aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = value), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", value)), vjust = 1) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() + theme(legend.position = "none") +
  xlab("Internalizing") +
  ylab("Externalizing")

## Now put these together
#ggarrange(cor.mat, conf.mat)
source("~/GitHub/adroseHelperScripts/R/afgrHelpFunc.R")
multiplot(cor.mat, conf.mat, cols=2)

# in.dat.mcp$predOne <- log(predOne$predict)
# predTwo <- predict(fit2)
# in.dat.mcp$predTwo <- log(predTwo$predict)

## Now create the figures here
p1 <- ggplot(tmp.dat, aes(x=cbcl_scr_syn_internal_r, y = log(cbcl_scr_syn_external_r + 1))) +
  geom_point(alpha=0) +
  theme_bw() +
  geom_hex(aes(fill=log(..count..)), bins=30) +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predNone)) +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predNone), color="green") +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predNoneFE), color="red") +
  xlab("Internal") + ylab("log(External)") +
  theme(legend.position = "NULL") +
  coord_cartesian(ylim = c(0,5))
p1
## Now do the exp version of p1 with the four models
p12 <- ggplot(tmp.dat, aes(x=cbcl_scr_syn_internal_r, y = cbcl_scr_syn_external_r)) +
  geom_point(alpha=0) +
  theme_bw() +
  geom_hex(aes(fill=log(..count..)), bins = 50, color="white") +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predNone_exp), color="green") +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predNoneFE_exp), color="red") +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predNone_lin), color="yellow") +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predNoneFE_lin), color="purple") +
  coord_cartesian(ylim=c(0,50), xlim=c(0,50)) +
  xlab("Internal") + ylab("External") +
  theme(legend.position = "NULL")
p12

p13 <- ggplot(tmp.dat, aes(x=cbcl_scr_syn_internal_r, y = log(cbcl_scr_syn_external_r + 1))) +
  geom_point(alpha=0) +
  theme_bw() +
  geom_hex(aes(fill=log(..count..)), bins=30) +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predNone)) +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predNone), color="green") +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predNoneFE), color="red") +
  xlab("Internal") + ylab("log(External)") +
  theme(legend.position = "NULL") +
  coord_cartesian(ylim = c(0,5))

p2 <- ggplot(tmp.dat, aes(x=cbcl_scr_syn_internal_r, y = log(cbcl_scr_syn_external_r + 1))) +
  geom_point(alpha=0) +
  theme_bw() +
  geom_hex(aes(fill=log(..count..))) +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predOne)) +
  xlab("Internal") + ylab("log(External)") +
  theme(legend.position = "NULL") +
  coord_cartesian(ylim = c(0,5))
p22 <- ggplot(tmp.dat, aes(x=cbcl_scr_syn_internal_r, y = cbcl_scr_syn_external_r)) +
  geom_point(alpha=0) +
  theme_bw() +
  geom_hex(aes(fill=log(..count..)), bins = 50, color="white") +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predOneExp), color="green") +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predOneFEExp), color="red") +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predOneREGaus), color="yellow") +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predOneFEGaus), color="purple") +
  coord_cartesian(ylim=c(0,50), xlim=c(0,50)) +
  xlab("Internal") + ylab("External") +
  theme(legend.position = "NULL")
p22
p23 <- ggplot(tmp.dat, aes(x=cbcl_scr_syn_internal_r, y = log(cbcl_scr_syn_external_r + 1))) +
  geom_point(alpha=0) +
  theme_bw() +
  geom_hex(aes(fill=log(..count..)), bins=30) +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predOne)) +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predOne), color="green") +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predOneFE), color="red") +
  xlab("Internal") + ylab("log(External)") +
  theme(legend.position = "NULL") +
  coord_cartesian(ylim = c(0,5))
p3 <- ggplot(tmp.dat, aes(x=cbcl_scr_syn_internal_r, y = log(cbcl_scr_syn_external_r + 1))) +
  geom_point(alpha=0) +
  theme_bw() +
  geom_hex(aes(fill=log(..count..))) +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predTwo)) +
  xlab("Internal") + ylab("log(External)") +
  theme(legend.position = "NULL") +
  coord_cartesian(ylim = c(0,5))
p32 <- ggplot(tmp.dat, aes(x=cbcl_scr_syn_internal_r, y = cbcl_scr_syn_external_r)) +
  geom_point(alpha=0) +
  theme_bw() +
  geom_hex(aes(fill=log(..count..)), bins = 40, color="white") +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predTwoExp), color="green") +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predTwoFEExp), color="red") +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predTwoREGaus), color="yellow") +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predTwoFEGaus), color="purple") +
  coord_cartesian(ylim=c(0,50), xlim=c(0,50)) +
  xlab("Internal") + ylab("External") +
  theme(legend.position = "NULL")
p32
p33 <- ggplot(tmp.dat, aes(x=cbcl_scr_syn_internal_r, y = log(cbcl_scr_syn_external_r + 1))) +
  geom_point(alpha=0) +
  theme_bw() +
  geom_hex(aes(fill=log(..count..)), bins=30) +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predTwo)) +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predTwo), color="green") +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predTwoFE), color="red") +
  xlab("Internal") + ylab("log(External)") +
  theme(legend.position = "NULL") +
  coord_cartesian(ylim = c(0,5))

ggarrange(p12, p22, p32, labels = c("0", "1", "2"))
ggarrange(p13, p23, p33, labels = c("0", "1", "2"))


## Now grab the trace plots
t1 <- bayesplot::mcmc_trace(mod.1cp, "alpha")
t2 <- bayesplot::mcmc_trace(result_case, "alpha[1]")
t3 <- bayesplot::mcmc_trace(result_case, "alpha[2]")

## Make a blank image
blank1 <- ggplot() + theme_void()

## Now plot these
r1 <- ggarrange(p1, blank1, blank1, nrow = 3)
r2 <- ggarrange(p2, t1, blank1, nrow=3)
r3 <- ggarrange(p3, t2, t3, nrow=3)
out.plot <- ggarrange(r1, r2, r3, ncol = 3)
ggsave(filename = "./reports/figure1.png", dpi=300, height = 8, width = 9, plot=out.plot)
ggsave(filename = "~/Downloads/figure1.png", dpi=300, height = 8, width = 9, plot=out.plot)


## NOw plot all of the parameters from the four various esimation techniques
## Start with the linear models
m1 <- summary(mod.base)$fixed[,"Estimate"]
m2 <- summary(mod.baseNR)$fixed[,"Estimate"]
m3 <- summary(mod.baseLin)$fixed[,"Estimate"]
m4 <- summary(mod.baseLinNR)$fixed[,"Estimate"]
plot.dat <- data.frame(parameters=c("Intercept", "Internal", "interview_age"), randNB=m1, fixedNB=m2, randGaus=m3, fixedGaus=m4)
plot.datBase <- reshape2::melt(plot.dat)
p4 <- ggplot(plot.datBase, aes(x=variable, y=value)) +
  geom_bar(stat="identity") +
  facet_wrap(parameters ~ .,scales="free")


## NOw do the 1-cp models
m1 <- summary(mod.1cp)$summary[c("a1", "a2", "b1", "b2", "b3", "alpha"),"mean"]
m2 <- summary(mod.1cpFE)$summary[c("a1", "a2", "b1", "b2", "b3", "alpha"),"mean"]
m3 <- summary(mod.1cpREG)$summary[c("a1", "a2", "b1", "b2", "b3", "alpha"),"mean"]
m4 <- summary(mod.1cpFEG)$summary[c("a1", "a2", "b1", "b2", "b3", "alpha"),"mean"]
plot.dat <- data.frame(parameters=c("Intercept1", "Intercept2","Slope1", "Slope2", "SlopeAge", "CP"), randNB=m1, fixedNB=m2, randGaus=m3, fixedGaus=m4)
plot.dat[6,2:5] <- plot.dat[6,2:5] / 6
plot.dat1CP <- reshape2::melt(plot.dat)
p5 <- ggplot(plot.dat1CP, aes(x=variable, y=value)) +
  geom_bar(stat="identity") +
  facet_wrap(parameters ~ .,scales="free")

## And finally 2-cp models
m1 <- summary(mod.2)$summary[c("a1", "a2", "a3","b1", "b2", "b3", "b4", "alpha[1]", "alpha[2]"),"mean"]
m2 <- summary(mod.2FE)$summary[c("a1", "a2", "a3","b1", "b2", "b3", "b4", "alpha[1]", "alpha[2]"),"mean"]
m3 <- summary(mod.2REGaus)$summary[c("a1", "a2", "a3","b1", "b2", "b3", "b4", "alpha[1]", "alpha[2]"),"mean"]
m4 <- summary(mod.2FEGaus)$summary[c("a1", "a2", "a3","b1", "b2", "b3", "b4", "alpha[1]", "alpha[2]"),"mean"]
plot.dat <- data.frame(parameters=c("Intercept1", "Intercept2","Intercept3","Slope1", "Slope2", "SlopeAge","Slope3", "CP1", "CP2"), randNB=m1, fixedNB=m2, randGaus=m3, fixedGaus=m4)
plot.dat[8:9,2:5] <- plot.dat[8:9,2:5] / 6
plot.dat[7,2:5] <- plot.dat[7,2:5] + plot.dat[5,2:5]
plot.dat2CP <- reshape2::melt(plot.dat)
p6 <- ggplot(plot.dat2CP, aes(x=variable, y=value)) +
  geom_bar(stat="identity") +
  facet_wrap(parameters ~ .,scales="free")


## Now put all of these together
ggarrange(p4, p5, p6, nrow=3)
