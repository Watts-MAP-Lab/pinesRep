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

## Grab all model permutations
mod.dv <- c("cbcl_scr_syn_internal_r", "cbcl_scr_syn_external_r", "cbcl_scr_syn_attention_r", "cbcl_scr_syn_internal_r", "cbcl_scr_syn_thought_r")
mod.iv <- mod.dv
iter <- 1:4
all.mods <- expand.grid(mod.dv, mod.iv, iter)
all.mods$Var1 <- as.character(all.mods$Var1)
all.mods$Var2 <- as.character(all.mods$Var2)
all.mods <- all.mods[-which(as.character(all.mods$Var1)==as.character(all.mods$Var2)),]
all.mods <- unique(all.mods)

## Find our models of interest
mod.int <- c(1,4)

## Load the 1 cp model 
mod.1cpEI <- readRDS("./data/brmsModsOut/model_rawX_NB_1CP_allmods_1.RDS")
mod.1cpIE <- readRDS("./data/brmsModsOut/model_rawX_NB_1CP_allmods_4.RDS")

## Now prep all of our plot vals
tmp.dat <- in.dat[complete.cases(in.dat[,c("cbcl_scr_syn_external_r", "cbcl_scr_syn_internal_r", "interview_age")]),]
sum.valsEI <- mod.1cpEI$out.sum
sum.valsIE <- mod.1cpIE$out.sum

## Now grab the predicted values
pred.vals.one <- sum.valsEI["alpha[1]","mean"] + sum.valsEI["beta[1]","mean"] * tmp.dat$cbcl_scr_syn_internal_r
pred.vals.two <- sum.valsEI["alpha[2]","mean"] + sum.valsEI["beta[2]","mean"] * tmp.dat$cbcl_scr_syn_internal_r
cp.val <- sum.valsEI["r","mean"]
index.lt <- which(tmp.dat$cbcl_scr_syn_internal_r<cp.val)
tmp.dat$predOne <- pred.vals.one
tmp.dat$predOne[-index.lt] <- pred.vals.two[-index.lt]

## Now do the same for IE
pred.vals.one <- sum.valsIE["alpha[1]","mean"] + sum.valsIE["beta[1]","mean"] * tmp.dat$cbcl_scr_syn_external_r
pred.vals.two <- sum.valsIE["alpha[2]","mean"] + sum.valsIE["beta[2]","mean"] * tmp.dat$cbcl_scr_syn_external_r
cp.val <- sum.valsIE["r","mean"]
index.lt <- which(tmp.dat$cbcl_scr_syn_external_r<cp.val)
tmp.dat$predOne2 <- pred.vals.one
tmp.dat$predOne2[-index.lt] <- pred.vals.two[-index.lt]

## Now create the figures here
p1 <- ggplot(tmp.dat, aes(x=cbcl_scr_syn_internal_r, y = log(cbcl_scr_syn_external_r + 1))) +
  geom_point(alpha=0) +
  theme_bw() +
  geom_hex(aes(fill=log(..count..)), bins=15) +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_internal_r, y=predOne)) +
  xlab("Internal") + ylab("log(External)") +
  theme(legend.position = "NULL") +
  coord_cartesian(ylim = c(0,5))

p2 <- ggplot(tmp.dat, aes(x=cbcl_scr_syn_external_r, y = log(cbcl_scr_syn_internal_r + 1))) +
  geom_point(alpha=0) +
  theme_bw() +
  geom_hex(aes(fill=log(..count..)), bins=15) +
  geom_line(data = tmp.dat, aes(x=cbcl_scr_syn_external_r, y=predOne2)) +
  ylab("log(Internal)") + xlab("External") +
  theme(legend.position = "NULL") +
  coord_cartesian(ylim = c(0,5))

out.plot <- ggarrange(p1, p2, labels = c("A","B"))
ggsave(filename = "./reports/figure3EIRevCompare.png", plot = out.plot, dpi = 300, width=7, height = 4, units="in")
