## This script will be used to produce figure 2 for the cp paper
## figure 2 will examine all of the cp models across all higher order constructs
## Each figure will show the linear portion of the logistic regression model
## and this will be isolated to all wave 1 data
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
iter.pats <- which(all.mods$Var1 %in% iso.vals & all.mods$Var2 %in% iso.vals & all.mods$Var3==1)
out.plots <- list()
out.plots2 <- list()
list.count <- 1
## Now go through and plot all of these data
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
  in.mod <- paste("./data/brmsModsOut/model_rawX_NB_allmods_", z, ".RDS", sep='')
  if(file.exists(in.mod)){
    in.mod <- readRDS(in.mod)
    ## Check if Rhat are lrage here -- just print a warning
    sum.vals <- rstan::summary(in.mod)$summary[c("a1", "a2", "b1", "b2", "b3", "alpha","sigma_p", "sigma_p2", "phi"),c("mean")]
    
    ## First obtain the inv logit vals
    inv.logit.vals <- LaplacesDemon::invlogit(sum.vals["alpha"] + data_jags$x * -5)
    ## Now plot the predicted values
    pred.vals.one <- sum.vals["a1"] + sum.vals["b1"] * data_jags$x
    pred.vals.two <- sum.vals["a2"] + sum.vals["b2"] * data_jags$x
    pred.vals.combo <- exp(pred.vals.one * (inv.logit.vals) + (1 - inv.logit.vals) * pred.vals.two)
    pred.vals.lin <- pred.vals.one * (inv.logit.vals) + (1 - inv.logit.vals) * pred.vals.two
    plot.df <- data.frame(x=data_jags$x, y=data_jags$y, pred = pred.vals.combo, pred.lin = pred.vals.lin)
    plot.df$slopeVals = "LT"
    plot.df$slopeVals[which(inv.logit.vals>.5)] <- "GT"
    ## Now make the plot of these data
    for.heatmap <- reshape2::melt(table(plot.df$x, plot.df$y))
    for.heatmap$valueLog <- log(for.heatmap$value)
    for.heatmap <- for.heatmap[which(for.heatmap$valueLog!=-Inf),]
    ggplot(for.heatmap, aes(x=Var1,y=Var2, fill=valueLog)) +
      geom_tile()
    max.val <- max(range(for.heatmap$Var1))
    max.valx <- max.val
    max.val <- max(range(for.heatmap$Var2))
    max.valy <- max.val
    
    x.lab.val <- all.mods$Var2[z]
    y.lab.val <- all.mods$Var1[z]
    title.lab <- paste("Wave: ", all.mods$Var3[z], sep='')
    p1 <- ggplot(plot.df, aes(x=x, y=y)) +
      #geom_violin(aes(x=xFactor, y),scale = "count") +
      #geom_tile(data = for.heatmap, aes(x=Var1, y=Var2, fill=valueLog), fill = NA, color = "black", size = 1) +
      geom_tile(data = for.heatmap, aes(x=Var1, y=Var2, fill=valueLog)) +
      #geom_jitter(alpha=.01) +
      theme_bw() +
      geom_point(alpha = 0)+
      geom_line(data = plot.df,aes(y=pred, group = slopeVals), linewidth=2.75, alpha=.75,color="grey") +
      geom_line(data = plot.df,aes(y=pred, group = slopeVals), linewidth=1.75, color="grey") +
      geom_line(data = plot.df,aes(y=pred, group = slopeVals), linewidth=1, color="black") +
      coord_cartesian(xlim=c(0,max.valx), ylim=c(0,max.valy)) +
      xlab(x.lab.val) + ylab(y.lab.val) +
      #geom_smooth(method="gam", method.args = list(family=poisson)) +
      theme(legend.position = "NULL") +
      ggtitle(title.lab)
    #p1 <- ggMarginal(p1, type="histogram",xparams = list(binwidth = 1), yparams = list(binwidth = 1))
    #p1 <- ggMarginal(p1, adjust=1.75)
    p1lin <- ggplot(plot.df, aes(x=x, y=log(y+1))) +
      #geom_point() +
      #geom_violin(aes(group = cut_width(x, 1))) +
      geom_point(alpha=0) +
      theme_bw() +
      geom_hex(aes(fill=log(..count..)), bins=15) +
      geom_line(data = plot.df,aes(y=pred.lin, group = slopeVals), linewidth=2.75, alpha=.75,color="grey") +
      geom_line(data = plot.df,aes(y=pred.lin, group = slopeVals), linewidth=1.75, color="grey") +
      geom_line(data = plot.df,aes(y=pred.lin, group = slopeVals), linewidth=1, color="black") +
      theme(legend.position = "NULL") +
      xlab(x.lab.val) + ylab(paste("Log(",y.lab.val,")", sep='')) +
      coord_cartesian(xlim=c(0,max.valx)) +
      ggtitle(title.lab)
    #p1lin <- ggMarginal(p1lin, adjust=1.75)
    
    out.plots[[list.count]] <- p1
    out.plots2[[list.count]] <- p1lin
    list.count <- list.count + 1
    
  }
}
## Now go ahead and make the ggarrange for all of these
all.pats <- all.mods[iter.pats,]
fig.one <- out.plots2[[which(all.pats[,1] == "cbcl_scr_syn_external_r" & all.pats[,2] == "cbcl_scr_syn_internal_r" & all.pats[,3] == 1)]]
fig.two <- out.plots2[[which(all.pats[,1] == "cbcl_scr_syn_external_r" & all.pats[,2] == "cbcl_scr_syn_thought_r" & all.pats[,3] == 1)]]
fig.thr <- out.plots2[[which(all.pats[,1] == "cbcl_scr_syn_external_r" & all.pats[,2] == "cbcl_scr_syn_attention_r" & all.pats[,3] == 1)]]
fig.fou <- out.plots2[[which(all.pats[,1] == "cbcl_scr_syn_internal_r" & all.pats[,2] == "cbcl_scr_syn_thought_r" & all.pats[,3] == 1)]]
fig.fiv <- out.plots2[[which(all.pats[,1] == "cbcl_scr_syn_internal_r" & all.pats[,2] == "cbcl_scr_syn_attention_r" & all.pats[,3] == 1)]]

## Now modify the figures ranges here
fig.one <- fig.one +theme_bw() + coord_cartesian(ylim=c(0,5)) + xlab("Internal") + ylab("log(External)") + theme(legend.position = "NULL") + ggtitle(NULL)
fig.two <- fig.two +theme_bw() + coord_cartesian(ylim=c(0,5)) + xlab("Thought") + ylab("log(External)") + theme(legend.position = "NULL") + ggtitle(NULL)
fig.thr <- fig.thr + theme_bw() + coord_cartesian(ylim = c(0,5), xlim = c(0,20)) + xlab("Attention") + ylab("log(External)") + theme(legend.position = "NULL") + ggtitle(NULL)
fig.fou <- fig.fou + theme_bw() + coord_cartesian(ylim = c(0,5), xlim=c(0,20)) + xlab("Thought") + ylab("log(Internal)") + theme(legend.position = "NULL") + ggtitle(NULL)
fig.fiv <- fig.fiv + theme_bw() + coord_cartesian(ylim = c(0,5)) + xlab("Attention") + ylab("log(Internal)") + theme(legend.position = "NULL") + ggtitle(NULL)
## Now add the clinical threshold delins here
fig.one <- fig.one + geom_rect(aes(xmin=12, xmax=Inf, ymin=-Inf, ymax=+Inf),fill='yellow', alpha=0.01)
fig.one <- fig.one + geom_rect(aes(xmin=16.5, xmax=Inf, ymin=-Inf, ymax=+Inf),fill='red', alpha=0.01)
fig.one <- fig.one + geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=log(12), ymax=+Inf),fill='yellow', alpha=0.01)
fig.one <- fig.one + geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=log(16.5), ymax=+Inf),fill='red', alpha=0.01)



## Now do the same for the raw units
fig.oner <- out.plots[[which(all.pats[,1] == "cbcl_scr_syn_external_r" & all.pats[,2] == "cbcl_scr_syn_internal_r" & all.pats[,3] == 1)]]
fig.twor <- out.plots[[which(all.pats[,1] == "cbcl_scr_syn_external_r" & all.pats[,2] == "cbcl_scr_syn_thought_r" & all.pats[,3] == 1)]]
fig.thrr <- out.plots[[which(all.pats[,1] == "cbcl_scr_syn_external_r" & all.pats[,2] == "cbcl_scr_syn_attention_r" & all.pats[,3] == 1)]]
fig.four <- out.plots[[which(all.pats[,1] == "cbcl_scr_syn_internal_r" & all.pats[,2] == "cbcl_scr_syn_thought_r" & all.pats[,3] == 1)]]
fig.fivr <- out.plots[[which(all.pats[,1] == "cbcl_scr_syn_internal_r" & all.pats[,2] == "cbcl_scr_syn_attention_r" & all.pats[,3] == 1)]]

## Now modify the figures ranges here
fig.oner <- fig.oner +theme_bw() + xlab("Internal") + ylab("External") + theme(legend.position = "NULL") + ggtitle(NULL)
fig.twor <- fig.twor +theme_bw() + xlab("Thought") + ylab("External") + theme(legend.position = "NULL") + ggtitle(NULL)
fig.thrr <- fig.thrr + theme_bw() + xlab("Attention") + ylab("External") + theme(legend.position = "NULL") + ggtitle(NULL)
fig.four <- fig.four + theme_bw() +  xlab("Thought") + ylab("Internal") + theme(legend.position = "NULL") + ggtitle(NULL)
fig.fivr <- fig.fivr + theme_bw() +  xlab("Attention") + ylab("Internal") + theme(legend.position = "NULL") + ggtitle(NULL)


## Now position these next to each other
out.plot <- ggarrange(
ggarrange(fig.one, fig.oner, labels = c("A","B")),
ggarrange(fig.two, fig.twor, labels = c("C","D")),
ggarrange(fig.thr, fig.thrr, labels = c("E","F")),
ggarrange(fig.fou, fig.four, labels = c("G","H")),
ggarrange(fig.fiv, fig.fivr, labels = c("I","J")),
nrow = 5)
ggsave(filename ="./reports/figure2.png", height=10, width=6, dpi=300)
ggsave(filename ="~/Downloads/figure2.png", height=10, width=6, dpi=300)


## Now do this same figure with all of the marginals
fig.oner <- ggMarginal(fig.oner, xparams = list(adjust=2), yparams = list(adjust=5))
fig.twor <- ggMarginal(fig.twor, xparams = list(adjust=3), yparams = list(adjust=2)) 
fig.thrr <- ggMarginal(fig.thrr, xparams = list(adjust=3), yparams = list(adjust=2)) 
fig.four <- ggMarginal(fig.four, xparams = list(adjust=3), yparams = list(adjust=2))
fig.fivr <- ggMarginal(fig.fivr, xparams = list(adjust=3), yparams = list(adjust=2))

out.plot <- ggarrange(
  ggarrange(fig.one, fig.oner, labels = c("A","B")),
  ggarrange(fig.two, fig.twor, labels = c("C","D")),
  ggarrange(fig.thr, fig.thrr, labels = c("E","F")),
  ggarrange(fig.fou, fig.four, labels = c("G","H")),
  ggarrange(fig.fiv, fig.fivr, labels = c("I","J")),
  nrow = 5)
ggsave(filename ="./reports/figure2b.png", height=10, width=6, dpi=300)
ggsave(filename ="~/Downloads/figure2b.png", height=10, width=6, dpi=300)


out.plot1 <- ggarrange(fig.one, fig.two, fig.thr, fig.fou, fig.fiv)

fig.oneM <- ggMarginal(fig.one, xparams=list(adjust=3), yparams = list(adjust=1))
fig.twoM <- ggMarginal(fig.two, xparams=list(adjust=3), yparams = list(adjust=1.75))
fig.thrM <- ggMarginal(fig.thr, xparams=list(adjust=3), yparams = list(adjust=1.75))
fig.fouM <- ggMarginal(fig.fou, xparams=list(adjust=3), yparams = list(adjust=1.75))
fig.fivM <- ggMarginal(fig.fiv, xparams=list(adjust=3), yparams = list(adjust=1.75))
plotg <- get_legend(fig.one + theme(legend.position = "right"))
out.plot1M <- ggarrange(fig.oneM, fig.twoM, fig.thrM, fig.fouM, fig.fivM, plotg)

ggsave(filename ="./reports/figure2c.png", plot = out.plot1M,height=5, width=6, dpi=300)
ggsave(filename ="~/Downloads/figure2c.png", plot = out.plot1M,height=5, width=6, dpi=300)
