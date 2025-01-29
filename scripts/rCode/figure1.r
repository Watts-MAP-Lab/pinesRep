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
stanmonitor <- c("a1","a2","a3","b1","b2","b3","b4","alpha","sigma_p","sigma_p2","phi")

result_case = stan(file="./scripts/stan_models/doubleChangePoint.stan",
                   data = data_jags, cores=2,chains=2,
                   pars = stanmonitor,
                   iter=20000, warmup = 10000, thin = 5,control = list(max_treedepth=12))

saveRDS(result_case, file="./data/brmsModsOut/two_cp_converge.RDS")
result_case = readRDS("./data/brmsModsOut/two_cp_converge.RDS")
tmp <- summary(result_case)
tmp <- tmp$summary

## Now organize all of these model predictions in the original data frame
tmp.dat$predNone <- sum.none$fixed["Intercept","Estimate"] + tmp.dat$cbcl_scr_syn_internal_r * sum.none$fixed["cbcl_scr_syn_internal_r","Estimate"]
## Now do the pred vals for the 1-cp model
sum.vals <- summary(mod.1cp)$summary
inv.logit.vals <- LaplacesDemon::invlogit(sum.vals["alpha","mean"]+ tmp.dat$cbcl_scr_syn_internal_r * -5)
## Now grab the predicted values
pred.vals.one <- sum.vals["a1","mean"] + sum.vals["b1","mean"] * tmp.dat$cbcl_scr_syn_internal_r
pred.vals.two <- sum.vals["a2","mean"] + sum.vals["b2","mean"] * tmp.dat$cbcl_scr_syn_internal_r
pred.vals.combo <- exp(pred.vals.one * (inv.logit.vals) + (1 - inv.logit.vals) * pred.vals.two)
pred.vals.lin <- pred.vals.one * (inv.logit.vals) + (1 - inv.logit.vals) * pred.vals.two
tmp.dat$predOne <- pred.vals.lin

## Do the same for the 2-cp model
## Create a matrix for the ordered inv logit vals
sum.vals <- summary(result_case)$summary
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

# in.dat.mcp$predOne <- log(predOne$predict)
# predTwo <- predict(fit2)
# in.dat.mcp$predTwo <- log(predTwo$predict)

## Now create the figures here
p1 <- ggplot(tmp.dat, aes(x=cbcl_scr_syn_internal_r, y = log(cbcl_scr_syn_external_r + 1))) +
  geom_point(alpha=0) +
  theme_bw() +
  geom_hex(aes(fill=log(..count..)), bins=15) +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predNone)) +
  xlab("Internal") + ylab("log(External)") +
  theme(legend.position = "NULL") +
  coord_cartesian(ylim = c(0,5))
p1
p2 <- ggplot(tmp.dat, aes(x=cbcl_scr_syn_internal_r, y = log(cbcl_scr_syn_external_r + 1))) +
  geom_point(alpha=0) +
  theme_bw() +
  geom_hex(aes(fill=log(..count..)), bins=15) +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predOne)) +
  xlab("Internal") + ylab("log(External)") +
  theme(legend.position = "NULL") +
  coord_cartesian(ylim = c(0,5))
p3 <- ggplot(tmp.dat, aes(x=cbcl_scr_syn_internal_r, y = log(cbcl_scr_syn_external_r + 1))) +
  geom_point(alpha=0) +
  theme_bw() +
  geom_hex(aes(fill=log(..count..)), bins=15) +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predTwo)) +
  xlab("Internal") + ylab("log(External)") +
  theme(legend.position = "NULL") +
  coord_cartesian(ylim = c(0,5))

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
