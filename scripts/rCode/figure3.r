## Figure 3 production script here
## This script will be used to produce the figure 3 results

## Load library(s)
library(ggplot2)
library(bayesplot)
library(ggpubr)
library(ggExtra)

## Now expand this to the all_mods runs
mod.dv <- c("cbcl_scr_syn_internal_r", "cbcl_scr_syn_external_r", "cbcl_scr_syn_attention_r", "cbcl_scr_syn_internal_r", "cbcl_scr_syn_thought_r", "cbcl_scr_syn_anxdep_r", "cbcl_scr_syn_withdep_r", "cbcl_scr_syn_somatic_r", "cbcl_scr_syn_aggressive_r", "cbcl_scr_syn_rulebreak_r")
mod.iv <- mod.dv
iter <- 1:4
all.mods <- expand.grid(mod.dv, mod.iv, iter)
all.mods$Var1 <- as.character(all.mods$Var1)
all.mods$Var2 <- as.character(all.mods$Var2)
all.mods <- all.mods[-which(as.character(all.mods$Var1)==as.character(all.mods$Var2)),]
all.mods <- unique(all.mods)
iter.val <- 1
index.val <- 1
out.dat <- NULL
no.mods <- NULL
## NOw isolate only the models of interest
iso.vals <- c("cbcl_scr_syn_internal_r", "cbcl_scr_syn_external_r", "cbcl_scr_syn_attention_r", "cbcl_scr_syn_internal_r", "cbcl_scr_syn_thought_r")
iter.pats <- which(all.mods$Var1 %in% iso.vals & all.mods$Var2 %in% iso.vals)
for(z in iter.pats){
  ## create the data
  i = all.mods[z,3]
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
  ## This line of code finds the smallest common multiple
  ## In case we want to consider resclaing these data
  cheapr::scm(as.numeric(apply(tmp.dat, 2, function(x) range(x, na.rm=TRUE))[2,130:135]))
  tmp.dat <- in.dat[complete.cases(tmp.dat[,c(all.mods$Var1[z], all.mods$Var2[z], "interview_age")]),]
  data_jags <- list(
    N = nrow(tmp.dat),
    x = tmp.dat[,all.mods$Var2[z]],
    y = tmp.dat[,all.mods$Var1[z]],
    x2 = tmp.dat$interview_age,
    #N_group = as.numeric(factor(tmp.dat$site_id_l)),
    #count_group = length(unique(as.numeric(factor(tmp.dat$site_id_l)))),
    #N_group2 = as.numeric(factor(tmp.dat$rel_family_id)),
    #count_group2 = length(unique(as.numeric(factor(tmp.dat$rel_family_id)))),
    beta = -5
  )
  ## Now read the brms model
  in.mod <- paste("./data/brmsModsOut/model_rawX_NB_allmods_", z, ".RDS", sep='')
  if(file.exists(in.mod)){
    in.mod <- readRDS(in.mod)
    ## Check if Rhat are lrage here -- just print a warning
    if(sum(rstan::summary(in.mod)$summary[,"Rhat"] > 1.2)>0){
      print("Warning")
      print(c(all.mods$Var2[z], all.mods$Var1[z], i))
    }
    ## First obtain the inv logit vals
    sum.vals <- rstan::summary(in.mod)$summary[c("a1", "a2", "b1", "b2", "b3", "alpha","sigma_p", "sigma_p2", "phi"),c("mean")]
    sum.valsL <- rstan::summary(in.mod)$summary[c("a1", "a2", "b1", "b2", "b3", "alpha","sigma_p", "sigma_p2", "phi"),c("2.5%")]
    sum.valsU <- rstan::summary(in.mod)$summary[c("a1", "a2", "b1", "b2", "b3", "alpha","sigma_p", "sigma_p2", "phi"),c("97.5%")]
    sum.vals["alpha"] <- sum.vals["alpha"]/5
    sum.valsL["alpha"] <- sum.valsL["alpha"]/5
    sum.valsU["alpha"] <- sum.valsU["alpha"]/5
    
    sum.vals <- c(sum.vals,sum.valsL, sum.valsU,all.mods$Var2[z], all.mods$Var1[z], i)
    out.dat <- rbind(sum.vals, out.dat)
    index.val <- index.val + 1
    
  }else{
    print(paste("No Mod: ", z))
    no.mods <- c(no.mods, z)
  }
}

out.dat <- data.frame(out.dat)
colnames(out.dat) <- c("Pre_Intercept", "Post_Intercept", "Pre_Slope", "Post_Slope", "Age_Slope", "ChangePoint","var_Site","var_family","phi",
                       "Pre_Intercept_Lower", "Post_Intercept_Lower", "Pre_Slope_Lower", "Post_Slope_Lower", "Age_Slope_Lower", "ChangePoint_Lower","var_Site_Lower","var_family_Lower","phi_lower",
                       "Pre_Intercept_Upper", "Post_Intercept_Upper", "Pre_Slope_Upper", "Post_Slope_Upper", "Age_Slope_Upper", "ChangePoint_Upper","var_Site_Upper","var_family_Upper","phi_Upper",
                       "Predictor","Outcome", "Wave")
write.csv(out.dat, file = "./data/outputSummaryStatsCpModels_AllExtend_RAWDAT.csv", quote=FALSE, row.names=FALSE)

## Now do the scaled values here
mod.dv <- c("cbcl_scr_syn_internal_r", "cbcl_scr_syn_external_r", "cbcl_scr_syn_attention_r", "cbcl_scr_syn_internal_r", "cbcl_scr_syn_thought_r", "cbcl_scr_syn_anxdep_r", "cbcl_scr_syn_withdep_r", "cbcl_scr_syn_somatic_r", "cbcl_scr_syn_aggressive_r", "cbcl_scr_syn_rulebreak_r")
mod.iv <- mod.dv
iter <- 1:4
all.mods <- expand.grid(mod.dv, mod.iv, iter)
all.mods$Var1 <- as.character(all.mods$Var1)
all.mods$Var2 <- as.character(all.mods$Var2)
all.mods <- all.mods[-which(as.character(all.mods$Var1)==as.character(all.mods$Var2)),]
all.mods <- unique(all.mods)
iter.val <- 1
index.val <- 1
out.dat <- NULL
no.mods <- NULL
## NOw isolate only the models of interest
iso.vals <- c("cbcl_scr_syn_internal_r", "cbcl_scr_syn_external_r", "cbcl_scr_syn_attention_r", "cbcl_scr_syn_internal_r", "cbcl_scr_syn_thought_r")
iter.pats <- which(all.mods$Var1 %in% iso.vals & all.mods$Var2 %in% iso.vals)
for(z in iter.pats){
  ## create the data
  i = all.mods[z,3]
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
  tmp.dat <- in.dat[complete.cases(tmp.dat[,c(all.mods$Var1[z], all.mods$Var2[z], "interview_age")]),]
  data_jags <- list(
    N = nrow(tmp.dat),
    x = tmp.dat[,all.mods$Var2[z]],
    y = tmp.dat[,all.mods$Var1[z]],
    x2 = tmp.dat$interview_age,
    #N_group = as.numeric(factor(tmp.dat$site_id_l)),
    #count_group = length(unique(as.numeric(factor(tmp.dat$site_id_l)))),
    #N_group2 = as.numeric(factor(tmp.dat$rel_family_id)),
    #count_group2 = length(unique(as.numeric(factor(tmp.dat$rel_family_id)))),
    beta = -5
  )
  ## Now read the brms model
  in.mod <- paste("./data/brmsModsOut/model_scaleXY_NB_allmods_", z, ".RDS", sep='')
  if(file.exists(in.mod)){
    in.mod <- readRDS(in.mod)
    ## Check if Rhat are lrage here -- just print a warning
    if(sum(rstan::summary(in.mod)$summary[,"Rhat"] > 1.2)>0){
      print("Warning")
      print(c(all.mods$Var2[z], all.mods$Var1[z], i))
    }
    ## First obtain the inv logit vals
    sum.vals <- rstan::summary(in.mod)$summary[c("a1", "a2", "b1", "b2", "b3", "alpha","sigma_p", "sigma_p2", "phi"),c("mean")]
    sum.valsL <- rstan::summary(in.mod)$summary[c("a1", "a2", "b1", "b2", "b3", "alpha","sigma_p", "sigma_p2", "phi"),c("2.5%")]
    sum.valsU <- rstan::summary(in.mod)$summary[c("a1", "a2", "b1", "b2", "b3", "alpha","sigma_p", "sigma_p2", "phi"),c("97.5%")]
    
    sum.vals["alpha"] <- sum.vals["alpha"]/5
    sum.valsL["alpha"] <- sum.valsL["alpha"]/5
    sum.valsU["alpha"] <- sum.valsU["alpha"]/5
    
    ## Now add the range of the x predictor
    min.x <- min(data_jags$x)
    max.x <- max(data_jags$x)
    range.x <- max.x-min.x
    ## Now look into other higher order moments
    all.db <- psych::describe(data_jags$x)
    all.db2 <- psych::describe(data_jags$y)
    ##Now do MV skew & Kurt
    all.mv <- psych::mardia(x = cbind(data_jags$x, data_jags$y), plot=FALSE)
    sum.vals <- c(sum.vals,sum.valsL, sum.valsU,all.mods$Var2[z], all.mods$Var1[z], i)
    out.dat <- rbind(sum.vals, out.dat)
    index.val <- index.val + 1
    
  }else{
    print(paste("No Mod: ", z))
    no.mods <- c(no.mods, z)
  }
}

out.dat <- data.frame(out.dat)
colnames(out.dat) <- c("Pre_Intercept", "Post_Intercept", "Pre_Slope", "Post_Slope", "Age_Slope", "ChangePoint","var_Site","var_family","phi",
                       "Pre_Intercept_Lower", "Post_Intercept_Lower", "Pre_Slope_Lower", "Post_Slope_Lower", "Age_Slope_Lower", "ChangePoint_Lower","var_Site_Lower","var_family_Lower","phi_lower",
                       "Pre_Intercept_Upper", "Post_Intercept_Upper", "Pre_Slope_Upper", "Post_Slope_Upper", "Age_Slope_Upper", "ChangePoint_Upper","var_Site_Upper","var_family_Upper","phi_Upper",
                       "Predictor","Outcome", "Wave")
write.csv(out.dat, file = "./data/outputSummaryStatsCpModels_AllExtend_SCALEDAT.csv", quote=FALSE, row.names=FALSE)

## Now load all model summary stats
sum.stats.raw <- read.csv("./data/outputSummaryStatsCpModels_AllExtend_RAWDAT.csv")
sum.stats.scale <- read.csv("./data/outputSummaryStatsCpModels_AllExtend_SCALEDAT.csv")
sum.stats.raw$predPat <- paste(sum.stats.raw$Predictor, sum.stats.raw$Outcome)
sum.stats.scale$predPat <- paste(sum.stats.scale$Predictor, sum.stats.scale$Outcome)

## Now merge these
sum.stats.merge <- merge(sum.stats.raw, sum.stats.scale, by=c("Predictor", "Outcome", "predPat", "Wave"), suffixes = c("_Raw", "_Scale"))

## Now isolate vals of interest
grep.vals <- c(1:4)
for(i in c("Pre_Intercept", "Post_Intercept", "Pre_Slope", "Post_Slope", "ChangePoint")){
  grep.vals <- append(grep.vals, grep(pattern=i, x = names(sum.stats.merge)))
}
## Now isolate data
iso.dat <- sum.stats.merge[,grep.vals]


## Now plot the location of the cp versus time?
p1 <- ggplot(iso.dat, aes(x=Wave, y= ChangePoint_Raw)) +
  geom_point() +
  geom_smooth(method="lm") +
  facet_grid(Predictor ~ Outcome, scales="free") +
  theme_bw()
ggplot(iso.dat, aes(x=Wave, y= Pre_Slope_Raw)) +
  geom_point() +
  geom_smooth(method="lm") +
  facet_grid(Predictor ~ Outcome, scales="free")
ggplot(iso.dat, aes(x=Wave, y= Post_Slope_Raw)) +
  geom_point() +
  geom_smooth(method="lm") +
  facet_grid(Predictor ~ Outcome, scales="free")

## Now melt the data so these can all be put on the same scale?
melt.dat <- iso.dat[,c("Wave", "Predictor", "Outcome", "Post_Slope_Raw", "Pre_Slope_Raw")]
melt.dat <- reshape2::melt(melt.dat, id.vars=c("Wave", "Predictor", "Outcome"))
p2 <- ggplot(melt.dat, aes(x=Wave, y= value, group = variable, fill=variable, colour = variable)) + 
  geom_point() + geom_smooth(method="lm") +
  facet_grid(Predictor ~ Outcome, scales="free") +
  theme_bw()
out.plot <- ggarrange(p1, p2, nrow = 2)
ggsave(filename = "./reports/figure3.png", plot = out.plot, dpi = 300, height = 11, width=9)
ggsave(filename = "~/Downloads/figure3.png", plot = out.plot, dpi = 300, height = 11, width=9)
