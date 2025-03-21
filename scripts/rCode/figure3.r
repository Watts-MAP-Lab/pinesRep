## Figure 3 production script here
## This script will be used to produce the figure 3 results

## Load library(s)
library(ggplot2)
library(bayesplot)
library(ggpubr)
library(ggExtra)

## Now expand this to the all_mods runs
## First thing is to load the data
in.dat <- read.csv("./data/ABCD_CBCL_BL.csv")
in.dat <- in.dat[complete.cases(in.dat[,c("cbcl_scr_syn_internal_r", "cbcl_scr_syn_external_r", "interview_age")]),]

## Library(s)
library(brms)
library(bayesplot)
library(rstan)
#library(mcp)
library(ggplot2)
library(ggpubr)

## Grab all model permutations
mod.dv <- c("cbcl_scr_syn_internal_r", "cbcl_scr_syn_external_r", "cbcl_scr_syn_attention_r", "cbcl_scr_syn_internal_r", "cbcl_scr_syn_thought_r")
mod.iv <- mod.dv
iter <- 1:4
all.mods <- expand.grid(mod.dv, mod.iv, iter)
all.mods$Var1 <- as.character(all.mods$Var1)
all.mods$Var2 <- as.character(all.mods$Var2)
all.mods <- all.mods[-which(as.character(all.mods$Var1)==as.character(all.mods$Var2)),]
all.mods <- unique(all.mods)

## Load our models
in.ols <- readRDS("./data/brmsModsOut/model_rawX_GAUS_NOCP_allmods_1.RDS")
in.1cpNB <- readRDS("./data/brmsModsOut/model_rawX_NB_1CP_allmods_1.RDS")
in.2cpNB <- readRDS("./data/brmsModsOut/model_rawX_NB_2CP_allmods_1.RDS")

## Start with the OLS
sum.valsEI <- in.ols$out.sum
in.dat$basicLMEREst <- sum.valsEI["alpha[1]","mean"] + sum.valsEI["beta[1]","mean"] * in.dat$cbcl_scr_syn_internal_r
## Now create the model residuals
in.dat$modelResidLMER <- in.dat$cbcl_scr_syn_external_r - in.ols$out.sum[grep(pattern = "mu", x = rownames(in.ols$out.sum)),"mean"]


## Now prep our data with the estimates from each of these models
sum.valsEI <- in.1cpNB$out.sum
pred.vals.one <- sum.valsEI["alpha[1]","mean"] + sum.valsEI["beta[1]","mean"] * in.dat$cbcl_scr_syn_internal_r
pred.vals.two <- sum.valsEI["alpha[2]","mean"] + sum.valsEI["beta[2]","mean"] * in.dat$cbcl_scr_syn_internal_r
cp.val <- sum.valsEI["r","mean"]
index.lt <- which(tmp.dat$cbcl_scr_syn_internal_r<cp.val)
in.dat$predOne <- ifelse(in.dat$cbcl_scr_syn_internal_r < cp.val, pred.vals.one, pred.vals.two)
#in.dat$predOne[-index.lt] <- pred.vals.two[-index.lt]
in.dat$predOneExp <- exp(in.dat$predOne)
## Now gorup these values
in.dat$predOneFit <- in.1cpNB$out.sum[grep(pattern = "mu", x = rownames(in.1cpNB$out.sum)),"mean"]
in.dat$predOneGroup <- ifelse(in.dat$cbcl_scr_syn_internal_r < cp.val, "LT", "GT")
in.dat$modelResid1CP <- in.dat$cbcl_scr_syn_external_r - exp(in.1cpNB$out.sum[grep(pattern = "mu", x = rownames(in.1cpNB$out.sum)),"mean"])

## Now do the 2cp model
sum.valsEI <- in.2cpNB$out.sum
pred.vals.one <- sum.valsEI["alpha[1]","mean"] + sum.valsEI["beta[1]","mean"] * in.dat$cbcl_scr_syn_internal_r
pred.vals.two <- sum.valsEI["alpha[2]","mean"] + sum.valsEI["beta[2]","mean"] * in.dat$cbcl_scr_syn_internal_r
pred.vals.thr <- sum.valsEI["alpha[3]","mean"] + sum.valsEI["beta[3]","mean"] * in.dat$cbcl_scr_syn_internal_r
cp1.val <- sum.valsEI["r[1]","mean"]
cp2.val <- sum.valsEI["r[2]","mean"]
index.lt2 <- which(tmp.dat$cbcl_scr_syn_internal_r<cp2.val)
index.lt1 <- which(tmp.dat$cbcl_scr_syn_internal_r<cp1.val)
in.dat$predTwo <- ifelse(in.dat$cbcl_scr_syn_internal_r < cp1.val, pred.vals.one, ifelse(in.dat$cbcl_scr_syn_internal_r < cp2.val, pred.vals.two, pred.vals.thr))
in.dat$predTwoExp <- exp(in.dat$predTwo)
in.dat$predTwoGroup <- ifelse(in.dat$cbcl_scr_syn_internal_r < cp1.val, "LT", ifelse(in.dat$cbcl_scr_syn_internal_r < cp2.val, "MD", "GT"))


## Now plot all of these
out.plot <- ggplot(in.dat, aes(x=cbcl_scr_syn_internal_r, y=cbcl_scr_syn_external_r)) +
  geom_hex(aes(fill=log(..count..)),bins=40) +
  geom_line(data = in.dat, aes(x=cbcl_scr_syn_internal_r, y=basicLMEREst), color="blue") +
  geom_line(data = in.dat, aes(x=cbcl_scr_syn_internal_r, y=predOneExp), color="orange") +
  geom_line(data = in.dat, aes(x=cbcl_scr_syn_internal_r, y=predTwoExp), color="purple") +
  coord_cartesian(ylim=c(0,50)) +
  theme_bw() +
  xlab("Internalizing") + ylab("Externalizing") + theme(legend.position = "NULL")

## Now save this figure
ggsave(filename = "./reports/figure2ModelComparison.png", plot = out.plot, dpi = 300, width=9, height = 6, units="in")
