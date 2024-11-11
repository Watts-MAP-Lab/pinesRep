## Load data
in.dat <- read.csv("./data/ABCD_CBCL_BL.csv")
in.dat <- in.dat[complete.cases(in.dat$cbcl_scr_syn_external_r),]
in.dat <- in.dat[complete.cases(in.dat$cbcl_scr_syn_internal_r),]

## Load library(s)
library(mgcv)
library(brms)
library(visreg)

## run spline model bayes
mod.one <- brm(cbcl_scr_syn_external_r ~ s(cbcl_scr_syn_internal_r), data = in.dat, cores = 5, thin = 10, iter = 6000, warmup = 2000)
mod.one2 <- brm(cbcl_scr_syn_external_r ~ s(cbcl_scr_syn_internal_r), data = in.dat, cores = 5, thin = 10, iter = 6000, warmup = 2000, family="poisson")

## run spline mgcv REML
mod.two <- mgcv::gam(cbcl_scr_syn_external_r ~ s(cbcl_scr_syn_internal_r), data = in.dat)
mod.two2 <- mgcv::gam(cbcl_scr_syn_external_r ~ s(cbcl_scr_syn_internal_r), data = in.dat, family="poisson")

## Plot these gams
visreg(mod.two2)
in.dat$modTwoPred <- predict(mod.two)
in.dat$modTwoPred2 <- exp(predict(mod.two2))
ggplot(in.dat, aes(x=cbcl_scr_syn_internal_r, y = cbcl_scr_syn_external_r)) +
  geom_point() +
  geom_line(aes(x=cbcl_scr_syn_internal_r, y = modTwoPred), color="blue") +
  geom_line(aes(x=cbcl_scr_syn_internal_r, y = modTwoPred2), color="red") +
  theme_bw()
  

## Now look into lin mod factors
in.dat$decileGroup <- cut(in.dat$cbcl_scr_syn_internal_r, quantile(in.dat$cbcl_scr_syn_internal_r, prob=seq(0,1,length=4), na.rm=TRUE), include.lowest = TRUE)
mod.thr <- glm(cbcl_scr_syn_external_r ~ cbcl_scr_syn_internal_r*decileGroup, data = in.dat, family="poisson")
visreg(mod.thr, "cbcl_scr_syn_internal_r", "decileGroup")
mod.thr2 <- glm(cbcl_scr_syn_external_r ~ cbcl_scr_syn_internal_r*decileGroup, data = in.dat)
visreg(mod.thr, "cbcl_scr_syn_internal_r", "decileGroup")
in.dat$modThrPred <- predict(mod.thr,type = "response")
in.dat$modThrPred2 <- predict(mod.thr2)

## Now plot these decile group fit vals
ggplot(in.dat, aes(x=cbcl_scr_syn_internal_r, y = cbcl_scr_syn_external_r)) +
  geom_point() +
  geom_line(aes(x = cbcl_scr_syn_internal_r, y = modThrPred, group=decileGroup))


## Now try quant reg
library(quantreg)
mod.fou <- rq(cbcl_scr_syn_external_r ~ cbcl_scr_syn_internal_r, tau = .2, data = in.dat)
mod.fou <- rq(cbcl_scr_syn_external_r ~ cbcl_scr_syn_internal_r, tau = .4, data = in.dat)
mod.fou <- rq(cbcl_scr_syn_external_r ~ cbcl_scr_syn_internal_r, tau=c(0.1,.25, .5, .7,.95), data = in.dat)
mod.fou <- rq(cbcl_scr_syn_external_r ~ cbcl_scr_syn_internal_r, tau=c(0.1,.2, .3, .7,.95), data = in.dat)
plot(mod.fou)

## Now do count quant reg
library(lqmm)
m1 <- lqm.counts(cbcl_scr_syn_external_r ~ cbcl_scr_syn_internal_r, data = in.dat, tau=c(.1),control = list(iterations=10000))
m2 <- lqm.counts(cbcl_scr_syn_external_r ~ cbcl_scr_syn_internal_r, data = in.dat, tau=c(.2),control = list(iterations=10000))
m3 <- lqm.counts(cbcl_scr_syn_external_r ~ cbcl_scr_syn_internal_r, data = in.dat, tau=c(.3))
m4 <- lqm.counts(cbcl_scr_syn_external_r ~ cbcl_scr_syn_internal_r, data = in.dat, tau=c(.4))
m5 <- lqm.counts(cbcl_scr_syn_external_r ~ cbcl_scr_syn_internal_r, data = in.dat, tau=c(.5))
m6 <- lqm.counts(cbcl_scr_syn_external_r ~ cbcl_scr_syn_internal_r, data = in.dat, tau=c(.6))
m7 <- lqm.counts(cbcl_scr_syn_external_r ~ cbcl_scr_syn_internal_r, data = in.dat, tau=c(.7))
m8 <- lqm.counts(cbcl_scr_syn_external_r ~ cbcl_scr_syn_internal_r, data = in.dat, tau=c(.8))
m9 <- lqm.counts(cbcl_scr_syn_external_r ~ cbcl_scr_syn_internal_r, data = in.dat, tau=c(.9))

## Mod overall
mod <- glm(cbcl_scr_syn_external_r ~ cbcl_scr_syn_internal_r, data = in.dat, family="poisson")
visreg(mod)
## Now make a list for easier operation
mod.list <- list(m1, m2, m3, m4, m5, m6, m7, m8, m9)
out.coef <- lapply(mod.list, coef)
pred.vals <- lapply(mod.list, function(x) exp(predict(x)))
pred.vals <- dplyr::bind_cols(pred.vals)
colnames(pred.vals) <- paste("predVal_", 1:9, sep='')
## Now plot
for.plot <- dplyr::bind_rows(out.coef)
for.plot$percentile <- 1:9
for.plot <- reshape2::melt(for.plot, id.vars="percentile")
library(ggplot2)
ggplot(for.plot, aes(x=percentile, y=value)) +
  geom_point() +
  facet_grid(variable~ .)
p1 <- ggplot(for.plot[which(for.plot$variable=="cbcl_scr_syn_internal_r"),], aes(x=percentile, y=value)) +
  geom_point() + ylab("Slope") +theme_bw()
p2 <- ggplot(for.plot[which(for.plot$variable!="cbcl_scr_syn_internal_r"),], aes(x=percentile, y=value)) +
  geom_point() + ylab("Intercept")+theme_bw()

## Plot these quant reg values
ggplot(in.dat, aes(x=cbcl_scr_syn_internal_r, y=cbcl_scr_syn_external_r)) + geom_point() +
  geom_quantile(quantiles = 0.9, color="red") +
  geom_quantile(quantiles = 0.5, color="blue") +
  geom_quantile(quantiles = 0.1, color="green") +
  theme_bw()
  
ggplot(in.dat, aes(x=cbcl_scr_syn_internal_r, y=cbcl_scr_syn_external_r)) + geom_point() +
  geom_quantile(quantiles = 0.9, color="red") +
  geom_quantile(quantiles = 0.5, color="blue") +
  geom_quantile(quantiles = 0.1, color="green") +
  theme_bw()
ggplot(in.dat, aes(x=cbcl_scr_syn_internal_r, y=cbcl_scr_syn_external_r)) + geom_point() +
  geom_quantile(quantiles = seq(0,.9,.1), color="red") +
  theme_bw()
## Now add the predVals to these data
in.dat2 <- dplyr::bind_cols(in.dat, pred.vals)
p3 <- ggplot(in.dat2, aes(x=cbcl_scr_syn_internal_r, y=cbcl_scr_syn_external_r)) + geom_point() +
  geom_line(aes(y=predVal_1)) +
  geom_line(aes(y=predVal_3)) +
  geom_line(aes(y=predVal_5)) +
  geom_line(aes(y=predVal_7)) +
  geom_line(aes(y=predVal_9)) +
  theme_bw() +
  coord_cartesian(xlim=c(0, 40), ylim=c(0, 60))
library(ggpubr)
p.poisson <- ggarrange(ggarrange(p1, p2, nrow = 1), p3, nrow=2)

## Now do the same for the gaussian rergession
mod.fou <- rq(cbcl_scr_syn_external_r ~ cbcl_scr_syn_internal_r, tau=c(0.1,.2, .3,.4,.5, .6, .7, .8, .9), data = in.dat)
for.plot.gaus <- reshape2::melt(t(coef(mod.fou)))
for.plot.gaus$percentile <- unlist(lapply(strsplit(as.character(for.plot.gaus$Var1),split = "="), function(x) as.numeric(x[2])))
p11 <- ggplot(for.plot.gaus[which(for.plot$variable=="cbcl_scr_syn_internal_r"),], aes(x=percentile, y=value)) +
  geom_point() + ylab("Slope") +theme_bw()
p22 <- ggplot(for.plot.gaus[which(for.plot$variable!="cbcl_scr_syn_internal_r"),], aes(x=percentile, y=value)) +
  geom_point() + ylab("Intercept")+theme_bw()
p33 <- ggplot(in.dat, aes(x=cbcl_scr_syn_internal_r, y=cbcl_scr_syn_external_r)) + geom_point() +
  geom_quantile(quantiles = seq(0,.9,.1), color="black") +
  theme_bw()
p.gaus <- ggarrange(ggarrange(p11, p22, nrow = 1), p33, nrow=2)
