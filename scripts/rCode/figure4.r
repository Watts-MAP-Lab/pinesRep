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
## Identify our models of int here
which(all.mods$Var1=="cbcl_scr_syn_external_r" & all.mods$Var2=="cbcl_scr_syn_internal_r")

## Load our models
in.1cpNB <- readRDS("./data/brmsModsOut/model_rawX_NB_1CP_allmods_1.RDS")
in.2cpNB <- readRDS("./data/brmsModsOut/model_rawX_NB_1CP_allmods_13.RDS")
in.3cpNB <- readRDS("./data/brmsModsOut/model_rawX_NB_1CP_allmods_25.RDS")
in.4cpNB <- readRDS("./data/brmsModsOut/model_rawX_NB_1CP_allmods_37.RDS")

## Now prep our data with the estimates from each of these models
sum.valsEI <- in.1cpNB$out.sum
pred.vals.one <- sum.valsEI["alpha[1]","mean"] + sum.valsEI["beta[1]","mean"] * in.dat$cbcl_scr_syn_internal_r
pred.vals.two <- sum.valsEI["alpha[2]","mean"] + sum.valsEI["beta[2]","mean"] * in.dat$cbcl_scr_syn_internal_r
cp.val <- sum.valsEI["r","mean"]
index.lt <- which(tmp.dat$cbcl_scr_syn_internal_r<cp.val)
in.dat$predOne <- ifelse(in.dat$cbcl_scr_syn_internal_r < cp.val, pred.vals.one, pred.vals.two)
#in.dat$predOne[-index.lt] <- pred.vals.two[-index.lt]
in.dat$predOneExp <- exp(in.dat$predOne)
in.dat$predOneGroup <- ifelse(in.dat$cbcl_scr_syn_internal_r < cp.val, "LT", "GT")

## Now gorup these values
in.dat$predOneFit <- in.1cpNB$out.sum[grep(pattern = "mu", x = rownames(in.1cpNB$out.sum)),"mean"]
in.dat$predOneGroup <- ifelse(in.dat$cbcl_scr_syn_internal_r < cp.val, "LT", "GT")
in.dat$modelResid1CP <- in.dat$cbcl_scr_syn_external_r - exp(in.1cpNB$out.sum[grep(pattern = "mu", x = rownames(in.1cpNB$out.sum)),"mean"])

## Now do second timepoint
sum.valsEI <- in.2cpNB$out.sum
pred.vals.one <- sum.valsEI["alpha[1]","mean"] + sum.valsEI["beta[1]","mean"] * in.dat$cbcl_scr_syn_internal_r
pred.vals.two <- sum.valsEI["alpha[2]","mean"] + sum.valsEI["beta[2]","mean"] * in.dat$cbcl_scr_syn_internal_r
cp.val <- sum.valsEI["r","mean"]
index.lt <- which(tmp.dat$cbcl_scr_syn_internal_r<cp.val)
in.dat$predTwo <- ifelse(in.dat$cbcl_scr_syn_internal_r < cp.val, pred.vals.one, pred.vals.two)
#in.dat$predOne[-index.lt] <- pred.vals.two[-index.lt]
in.dat$predTwoExp <- exp(in.dat$predTwo)
in.dat$predTwoGroup <- ifelse(in.dat$cbcl_scr_syn_internal_r < cp.val, "LT", "GT")


sum.valsEI <- in.3cpNB$out.sum
pred.vals.one <- sum.valsEI["alpha[1]","mean"] + sum.valsEI["beta[1]","mean"] * in.dat$cbcl_scr_syn_internal_r
pred.vals.two <- sum.valsEI["alpha[2]","mean"] + sum.valsEI["beta[2]","mean"] * in.dat$cbcl_scr_syn_internal_r
cp.val <- sum.valsEI["r","mean"]
index.lt <- which(tmp.dat$cbcl_scr_syn_internal_r<cp.val)
in.dat$predThr <- ifelse(in.dat$cbcl_scr_syn_internal_r < cp.val, pred.vals.one, pred.vals.two)
#in.dat$predOne[-index.lt] <- pred.vals.two[-index.lt]
in.dat$predTheExp <- exp(in.dat$predThr)
in.dat$predThrGroup <- ifelse(in.dat$cbcl_scr_syn_internal_r < cp.val, "LT", "GT")


sum.valsEI <- in.4cpNB$out.sum
pred.vals.one <- sum.valsEI["alpha[1]","mean"] + sum.valsEI["beta[1]","mean"] * in.dat$cbcl_scr_syn_internal_r
pred.vals.two <- sum.valsEI["alpha[2]","mean"] + sum.valsEI["beta[2]","mean"] * in.dat$cbcl_scr_syn_internal_r
cp.val <- sum.valsEI["r","mean"]
index.lt <- which(tmp.dat$cbcl_scr_syn_internal_r<cp.val)
in.dat$predFou <- ifelse(in.dat$cbcl_scr_syn_internal_r < cp.val, pred.vals.one, pred.vals.two)
#in.dat$predOne[-index.lt] <- pred.vals.two[-index.lt]
in.dat$predFouExp <- exp(in.dat$predFou)
in.dat$predFouGroup <- ifelse(in.dat$cbcl_scr_syn_internal_r < cp.val, "LT", "GT")


## Now get all of the cp val locations
sum.valsEI <- in.1cpNB$out.sum
cp.val1 <- sum.valsEI["r","mean"]
sum.valsEI <- in.2cpNB$out.sum
cp.val2 <- sum.valsEI["r","mean"]
sum.valsEI <- in.3cpNB$out.sum
cp.val3 <- sum.valsEI["r","mean"]
sum.valsEI <- in.4cpNB$out.sum
cp.val4 <- sum.valsEI["r","mean"]

## pred df
plot.df <- data.frame(cp.vals = c(cp.val1, cp.val2, cp.val3, cp.val4), wave=c("BP", "W1", "W2", "W3"))


## Now plot all of these
out.plot <- ggplot(in.dat, aes(x=cbcl_scr_syn_internal_r, y=log(cbcl_scr_syn_external_r+1))) +
  #geom_hex(aes(fill=log(..count..)),bins=20) +
  geom_line(data = in.dat, aes(x=cbcl_scr_syn_internal_r, y=predOne, group = predOneGroup), linetype="solid") +
  geom_vline(data = plot.df, aes(xintercept = cp.vals, x=NULL, y=NULL, linetype=wave)) +
  scale_linetype_manual(values=c("solid","longdash","dotdash","dotted")) +
  geom_line(data = in.dat, aes(x=cbcl_scr_syn_internal_r, y=predTwo, group = predTwoGroup), linetype="longdash") +
  geom_line(data = in.dat, aes(x=cbcl_scr_syn_internal_r, y=predThr, group = predThrGroup), linetype="dotdash") +
  geom_line(data = in.dat, aes(x=cbcl_scr_syn_internal_r, y=predFou, group = predFouGroup), linetype="dotted") +
  coord_cartesian(ylim=c(0,5), xlim=c(0,7)) +
  theme_bw() +
  xlab("Internal") + ylab("log(External)") +
  scale_x_continuous(breaks = 1:10)
ggsave(filename = "./reports/figure4WaveCompare.png", plot = out.plot, width = 6, height=5, units = "in", dpi = 300)
