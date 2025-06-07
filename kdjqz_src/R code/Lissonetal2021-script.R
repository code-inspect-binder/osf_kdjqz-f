###############################################################
##### Paula Lissón - lissonh@gmail.com ########################
###############################################################
##### Outline of the script:
##### - 1. Load packages
##### - 2. Data preprocessing and descriptive stats
##### - 3. Data analysis (Bayesian mixed effects models)
##### - 4. Activation-based model (fit, summary, plots)
##### - 5. Direct-access model (fit, summary, plots)
##### - 6. Posterior predictive checks (PPCs) for both models
##### - 7. Crossvalidation for both models, elpd plots & table 
###############################################################

###############################################################
## 1. Load packages
###############################################################
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(plyr)
library(ggplot2)
library(reshape2)
library(MASS)
library(gtools)
library(scales)
library(tidyverse)
library(gridExtra)
library(brms)
library(bayesplot)
library(purrr)
theme_set(theme_bw())
options(scipen=999)

###############################################################
## 2. Preprocessing of the data                            ####
###############################################################
data <- read.table("data/RC_LT.csv", sep=",", header = TRUE, stringsAsFactors = FALSE)

data <- filter(data, subj != "54002" & subj !=  "54015" & subj != "54065"
               & subj != "54047" & subj !="54050" & subj != "54051" & subj != "54054"
               & subj != "54056" & subj !="54059" & subj != "54060" & subj !="54061"
               & subj != "54062" & subj !="54044" & subj != "54023" & subj !="54024") %>% 
  filter(subj != "54041" & subj !=  "54007" & subj != "54045"  & subj != "54052"
         & subj != "54048"& subj != "54049" & subj !="54046") 

# Sum listening times (LT) per item, subj, and sentence type
data$rt <- as.integer(data$rt)
data <- data %>% 
  group_by(item, sent, subj, grp, acc) %>% 
  summarise_at(vars(rt), sum) %>% 
  filter(rt > 200) # outliers
  
# generate relevant columns 
data$rctype_meff <- ifelse(data$sent=="SS",-1,1) # Subject relatives are -1
data$group_meff <- ifelse(data$grp=="EC",-1,1) # Controls are -1
data$acc_lab <- ifelse(data$acc==0,"Incorrect","Correct")
data$rctype <- ifelse(data$sent=="SS","SR","OR")
data$group <- ifelse(data$group_meff==-1,"Controls", "IWAs")

# The actual answer given (1 for SR, 2 for OR)
data$answer <- ifelse( (data$rctype_meff==-1 & data$acc==1) | (data$rctype_meff==1 & data$acc==0), 1, 2)


###############################################################
##  Descriptive stats                 
###############################################################

# SE
se <- function(x) {
  return(sqrt(var(x)/length((x))))
}

# binomial SE
se.bin <- function(x){
  n.success = sum(x) 
  n         = length(x)
  p         = n.success / n
  return(sqrt(p*(1-p)/n))
}

# compute accuracy per condition and group
acc <- ddply(data, .(rctype,group), summarize,
             mean.acc = 100*mean(acc),
             se.lower = mean.acc - 100*(se.bin(acc)),
             se.upper = mean.acc + 100*(se.bin(acc))) 

# rename factors
acc$rctype_lab <- factor(ifelse(acc$rctype=="SR","SR", "OR"),levels=c("SR", "OR"))
acc$group_lab  <- factor(ifelse(acc$group=="Controls","Controls", "IWAs"),levels=c("Controls", "IWAs"))

# compute mean LT per condition, group and accuracy
LT <- ddply(data, .(rctype,group,acc_lab), summarize,
            mean.rt = mean(rt),
            se.lower = mean.rt - se(rt),
            se.upper = mean.rt + se(rt))
# rename factors
LT$acc_lab    <- factor(ifelse(LT$acc_lab=="Incorrect","Incorrect", "Correct"),levels=c("Correct","Incorrect"))
LT$rctype_lab <- factor(ifelse(LT$rctype=="SR","SR", "OR"),levels=c("SR", "OR"))
LT$group_lab  <- factor(ifelse(LT$group=="Controls","Controls", "IWAs"),levels=c("Controls", "IWAs"))

# figure 2
#pdf("figures/rawACC.pdf")
ggplot(data =acc,
       aes(x=rctype_lab, y=mean.acc, group=group, color=group)) +
  geom_errorbar(aes(ymin=se.lower, ymax=se.upper),
                width=0.1,size=.4) + ylim(60, 100) +
  geom_line(aes(linetype=group),size=.4) + geom_point(size=2) +
  scale_color_manual(values=c("darkgrey","black")) + 
  labs(y = "%", title = "Mean accuracy")  + xlab("Condition")
labs(group="Relative clause") +
  theme(legend.position="bottom")

# figure 3
#pdf("figures/rawLT.pdf", width = 8, height=6)
ggplot(data = LT.noacc,
       aes(x = rctype_lab, y= mean.rt, group = group, color = group)) +
  geom_errorbar(aes(ymin = se.lower, ymax = se.upper),
                width=.1,size=.4) + ylim(1000, 2600) +
  geom_line(aes(linetype =group), size=.4) + geom_point(size=2) +
  geom_line(aes(linetype = group), size=.4) +
  scale_color_manual(values=c("darkgrey","black")) + 
  labs(y = "ms", title = "Listening times")+  xlab("Condition")

###############################################################
##   3 Data Analysis                               ############
###############################################################
# priors for accuracies
priors_acc <- c(set_prior("normal(0, 1)", class = "Intercept"),
                set_prior("normal(0, .5)", class = "b", coef="rctype_meff"),
                set_prior("normal(0, .5)", class = "b", coef="group_meff"),
                set_prior("normal(0, .1)", class = "sd"),
                set_prior("lkj(2)", class = "cor"))
# priors for listening times
priors_lt <- c(set_prior("normal(7, .6)", class = "Intercept"),
               set_prior("normal(0, .5)", class = "b", coef="rctype_meff"),
               set_prior("normal(0, .5)", class = "b", coef="group_meff"),
               set_prior("normal(0, .1)", class = "sd"),
               set_prior("normal(0, .5)", class = "sigma"),
               set_prior("lkj(2)", class = "cor"))

# logit model for accuracies 
SROR.acc <- brm(acc ~ rctype_meff+group_meff+rctype_meff*group_meff +
                  (1+rctype_meff|subj) + (1+group_meff|item), data, 
                family = bernoulli(), prior = priors_acc, control = list(adapt_delta = 0.9))

# extract posterior samples 
SROR.acc.posterior <- posterior_samples(SROR.acc, "^b") %>%
  mutate(RC = (inv.logit(b_Intercept + b_rctype_meff) - inv.logit(b_Intercept - b_rctype_meff))*100,
         group = (inv.logit(b_Intercept + b_group_meff) - inv.logit(b_Intercept - b_group_meff))*100,
         `RC × group` = (inv.logit(b_Intercept + `b_rctype_meff:group_meff`) - inv.logit(b_Intercept - `b_rctype_meff:group_meff`))*100)  

# plot posterior samples 
# Figure 4a
posterior.acc.plot <- mcmc_areas(SROR.acc.posterior[,5:7], prob = 0.80, prob_outer = 0.95,point_est = "mean") +
  geom_vline(xintercept=0, linetype="dashed", size=.5) +
  labs(x = "Estimates (%)",title = "A. Accuracy") +
  theme_default(base_family = getOption("bayesplot.base_family", "sans")) 

# lognormal model for LT
SROR.LT <-brm(rt ~ rctype_meff+group_meff + rctype_meff*group_meff+
                (1+rctype_meff|subj) + (1+group_meff|item), 
              data,family=lognormal(),prior=priors_lt)
# extract posterior samples
SROR.LT.posterior <-posterior_samples(SROR.LT, "^b") %>%
  mutate(RC = exp(b_Intercept + b_rctype_meff) - exp(b_Intercept - b_rctype_meff),
         group = exp(b_Intercept + b_group_meff) - exp(b_Intercept - b_group_meff),
         `RC × group` = exp(b_Intercept + `b_rctype_meff:group_meff`) - exp(b_Intercept - `b_rctype_meff:group_meff`))

# plot posterior
# Figure 4b
posterior.LT.plot <- mcmc_areas(SROR.LT.posterior[,5:7], prob = 0.80, prob_outer = 0.95,point_est = "mean") +
  geom_vline(xintercept=0, linetype="dashed", color="black", size=.5) +
  labs(x = "Estimates (ms)",title = "B. Listening times") +
  bayesplot::theme_default(base_family = getOption("bayesplot.base_family", "sans")) 

# Figure 4a & 4b 
pdf("figures/estimates.pdf", width=8, height=6)
grid.arrange(posterior.acc.plot, posterior.LT.plot, nrow=1)
dev.off()

# Accuracy model: extracting and backtransforming the effects to percentages 
SROR.acc.posterior <-posterior_samples(SROR.acc, "^b")

rctype.acc <- SROR.acc.posterior  %>%
  mutate(pred = inv.logit(b_Intercept + b_rctype_meff) - inv.logit(b_Intercept - b_rctype_meff))  %>%
  summarize(Estimate = (round(mean(pred),2)*100) %>%  paste("%"),
            `95% CrI` = paste((round(quantile(pred,c(0.025,.975)),2)*100),
                              collapse=", ") %>%
              paste0("CrI: [",.,"]"))

group.acc <- SROR.acc.posterior  %>%
  mutate(pred = inv.logit(b_Intercept + b_group_meff) - inv.logit(b_Intercept - b_group_meff))  %>%
  summarize(Estimate = (round(mean(pred),2)*100) %>%  paste("%"),
            `95% CrI` = paste((round(quantile(pred,c(0.025,.975)),2)*100),
                              collapse=", ") %>%
              paste0("CrI: [",.,"]"))

inter.acc <- SROR.acc.posterior  %>%
  mutate(pred = inv.logit(b_Intercept + `b_rctype_meff:group_meff`) - inv.logit(b_Intercept -`b_rctype_meff:group_meff`))  %>%
  summarize(Estimate = (round(mean(pred),2)*100) %>%  paste("%"),
            `95% CrI` = paste((round(quantile(pred,c(0.025,.975)),2)*100),
                              collapse=", ") %>%
              paste0("CrI: [",.,"]"))

## LT model: extracting and backtransforming the effects to ms
SROR.LT.posterior <-posterior_samples(SROR.LT, "^b")
rctype.LT <- SROR.LT.posterior  %>%
  mutate(pred = exp(b_Intercept + b_rctype_meff) - exp(b_Intercept - b_rctype_meff))  %>%
  summarize(Estimate = mean(pred) %>%
              round %>% paste("ms"),
            `95% CrI` = paste(round(quantile(pred,c(0.025,.975)),0),
                              collapse=", ") %>%
              paste0("CrI: [",.,"]"))

group.LT <- SROR.LT.posterior  %>%
  mutate(pred = exp(b_Intercept + b_group_meff) - exp(b_Intercept - b_group_meff))  %>%
  summarize(Estimate = mean(pred) %>%
              round %>% paste("ms"),
            `95% CrI` = paste(round(quantile(pred,c(0.025,.975)),0),
                              collapse=", ") %>%
              paste0("CrI: [",.,"]"))

inter.LT <- SROR.LT.posterior  %>%
  mutate(pred = exp(b_Intercept + `b_rctype_meff:group_meff`) - exp(b_Intercept -`b_rctype_meff:group_meff`))  %>%
  summarize(Estimate = mean(pred) %>%
              round %>% paste("ms"),
            `95% CrI` = paste(round(quantile(pred,c(0.025,.975)),0),
                              collapse=", ") %>%
              paste0("CrI: [",.,"]"))

###############################################################
###############################################################
##   4. Activation-based model                   ##############
###############################################################
###############################################################
### load the model (already saved)
ACT.fit.real <- readRDS("./fits/ACT/ACT_fit_real.rds")
###############################################################
ACT.params.short <- c("alpha", "beta", "sigma_e", "tau_u", "tau_w", "L_u", "L_w")
ACT.params.long <- c(
  "alpha[1]",
  "alpha[2]",
  "beta[1]",
  "beta[2]",
  "beta[3]",
  "beta[4]",
  "beta[5]",
  "beta[6]",
  "beta[7]",
  "tau_u[1]",
  "tau_u[2]",
  "tau_u[3]",
  "tau_u[4]",
  "tau_w[1]",
  "tau_w[2]",
  "tau_w[3]",
  "tau_w[4]",
  "L_u[1,1]",
  "L_u[2,1]",
  "L_u[2,2]",
  "L_u[3,1]",
  "L_u[3,2]",
  "L_u[3,3]",
  "L_u[4,1]",
  "L_u[4,2]",
  "L_u[4,3]",
  "L_u[4,4]",
  "L_w[1,1]",
  "L_w[2,1]",
  "L_w[2,2]",
  "L_w[3,1]",
  "L_w[3,2]",
  "L_w[3,3]",
  "L_w[4,1]",
  "L_w[4,2]",
  "L_w[4,3]",
  "L_w[4,4]",
  "sigma_e")

###############################################################
## Activation - fit real data 
###############################################################
###############################################################
# create list of data for Stan
ACT.data <- within(list(),
                   {
                     N_obs<-nrow(data)
                     RT <- data$rt
                     N_choices <- length(unique(data$answer))
                     winner <- data$answer
                     rctype <- data$rctype_meff
                     group <- data$group_meff
                     group_v <- data$group_meff
                     
                     subj <- as.integer(factor(data$subj))
                     item <- as.integer(factor(data$item))
                     
                     N_subj <- length(unique(data$subj))
                     N_item <- length(unique(data$item))
                     
                     n_u <- 4 # no. subj reff
                     n_w <- 4 # no. item reff
                   }
)

# fit the model
ACT.fit.real <- stan(data=ACT.data,file="./models/activation_based.stan",
                      iter=6000,chains=3,control = list(adapt_delta = 0.9))

# summary of the model (can be used with ACT.params.long or with ACT.params.short)
ACT.fit.summary <-print(summary(ACT.fit.real,pars=ACT.params.long, probs=c(0.025,0.975))$summary)

# Save the fit 
# filename <- paste("./fits/act_fit_real_", Sys.time(), ".rds", sep="")
# filename <- gsub(":","",gsub(" ","_",filename))
# saveRDS(ACT.fit.shift, filename)

###############################################################
## Activation - plot accumulators #############################
###############################################################
ACT.accumulators <- data.frame(
  rt_1 = colMeans(rstan::extract(ACT.fit.real)$rt_1),
  rt_2 = colMeans(rstan::extract(ACT.fit.real)$rt_2),
  gen_RT = colMeans(rstan::extract(ACT.fit.real)$gen_RT),
  group = colMeans(rstan::extract(ACT.fit.real)$gen_group),
  rctype = colMeans(rstan::extract(ACT.fit.real)$gen_rctype)
)

ACT.accumulators$group <- ifelse(ACT.accumulators$group=="-1","controls","IWAs")
ACT.accumulators$rctype <- ifelse(ACT.accumulators$rctype=="-1","SR","OR")

# gather is used to get the dataframe in proper format for plots 
controls_SR <- filter(ACT.accumulators, group=="controls"&rctype=="SR") %>% 
  gather(accumulator,rt,rt_1:rt_2,factor_key=TRUE)  
controls_OR <- filter(ACT.accumulators, group=="controls"&rctype=="OR") %>% 
  gather(accumulator,rt,rt_1:rt_2,factor_key=TRUE)
IWA_SR <- filter(ACT.accumulators, group=="IWAs"&rctype=="SR") %>% 
  gather(accumulator,rt,rt_1:rt_2,factor_key=TRUE)
IWA_OR <- filter(ACT.accumulators, group=="IWAs"&rctype=="OR") %>% 
  gather(accumulator,rt,rt_1:rt_2,factor_key=TRUE)


## plot accumulators, first accumulator is SR interpretation, second is OR interpretation
### controls SR 
controls_SR$accumulator <- ifelse(controls_SR$accumulator=="rt_1","SR","OR")
mean_controls_SR <- ddply(controls_SR, "accumulator", summarise, rt.mean=round(mean(rt)))

p1 <- ggplot(controls_SR, aes(rt, fill = accumulator)) + 
  geom_density(adjust=2,alpha=.7) + labs(title="a) SR - controls")  + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(breaks = c(1000,2000,4000,6000,8000,12000)) +
  xlab("finishing times (ms)") + theme(axis.text=element_text(size=5), 
                                       legend.key.size = unit(0.3, "cm"), legend.position=c(0.80, 0.75)) +
  scale_fill_manual(values=c("lightgrey","darkgrey")) +
  geom_vline(data=mean_controls_SR, aes(xintercept=rt.mean, colour=accumulator),
             linetype="dashed", size=0.8) + scale_color_manual(values=c("lightgrey","darkgrey"))

# controls OR
controls_OR$accumulator <- ifelse(controls_OR$accumulator=="rt_1","SR","OR")
mean_controls_OR <- ddply(controls_OR, "accumulator", summarise, rt.mean=round(mean(rt)))

p2 <- ggplot(controls_OR, aes(rt, fill = accumulator)) + 
  geom_density(adjust=2,alpha=.7) + labs(title="b) OR - controls")  + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(breaks = c(1000,2000,4000,6000,8000,12000)) +
  xlab("finishing times (ms)") + theme(axis.text=element_text(size=5), 
                                       legend.key.size = unit(0.3, "cm"), legend.position=c(0.80, 0.75)) +
  scale_fill_manual(values=c("lightgrey","darkgrey")) +
  geom_vline(data=mean_controls_OR, aes(xintercept=rt.mean, colour=accumulator),
             linetype="dashed", size=0.8) + scale_color_manual(values=c("lightgrey","darkgrey"))

# iwa SR
IWA_SR$accumulator <- ifelse(IWA_SR$accumulator=="rt_1","SR","OR")
mean_IWA_SR <- ddply(IWA_SR, "accumulator", summarise, rt.mean=round(mean(rt)))

p3 <- ggplot(IWA_SR, aes(rt, fill = accumulator)) + 
  geom_density(adjust=2,alpha=.7) + labs(title="c) SR - IWA")  + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = c(1000,2000,4000,6000,8000,12000)) +
  xlab("finishing times (ms)") + theme(axis.text=element_text(size=5), 
                                       legend.key.size = unit(0.3, "cm"), legend.position=c(0.80, 0.75)) +
  scale_fill_manual(values=c("lightgrey","darkgrey")) +
  geom_vline(data=mean_IWA_SR, aes(xintercept=rt.mean, colour=accumulator),
             linetype="dashed", size=0.8) + scale_color_manual(values=c("lightgrey","darkgrey"))

# iwa OR
IWA_OR$accumulator <- ifelse(IWA_OR$accumulator=="rt_1","SR","OR")
mean_IWA_OR <- ddply(IWA_OR, "accumulator", summarise, rt.mean=round(mean(rt)))

p4 <- ggplot(IWA_OR, aes(rt, fill = accumulator)) + 
  geom_density(adjust=2,alpha=.7) + labs(title="d) OR - IWA")  + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(breaks = c(1000,2000,4000,6000,8000,12000)) +
  xlab("finishing times (ms)") + theme(axis.text=element_text(size=5), 
                                       legend.key.size = unit(0.3, "cm"), legend.position=c(0.80, 0.75)) +
  scale_fill_manual(values=c("lightgrey","darkgrey")) +
  geom_vline(data=mean_IWA_OR, aes(xintercept=rt.mean, colour=accumulator),
             linetype="dashed", size=0.8) + scale_color_manual(values=c("lightgrey","darkgrey"))

# Figure 6
png(filename="accumulators.png", res=300, width=1800, height=1100) 
grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()


# extract sigma
sigma_i <- rstan::extract(ACT.fit.real)$sigma_i
sigma_c <- rstan::extract(ACT.fit.real)$sigma_c

sigma_act <- data.frame(probability = c(sigma_i,sigma_c),
                        group = c(rep("IWA",length(sigma_i)),
                                  rep("control",length(sigma_c))))
sigma_act_mean <- ddply(sigma_act, .(group), summarize,
                        mean = mean(probability),2)

# Figure 7
pdf("figures/sigma_act.pdf", width= 8, height=6)
ggplot(sigma_act, aes(probability, fill = group)) +
  geom_density(adjust=2,alpha=.5)  + xlab("parameter estimate")+
  geom_vline(data=sigma_act_mean, aes(xintercept=mean, colour=group), linetype="dashed", size=0.8) + 
  scale_color_manual(values=c("grey","black")) +
  theme_bw() + scale_fill_brewer(palette="Greys") + labs(title=expression("Posterior distribution of"~sigma~"(activation-based)"))
dev.off()

###############################################################
###############################################################
##   5. Direct-access model                      ##############
###############################################################
###############################################################
####### load the model (previously saved)
DA.fit.real <- readRDS("./fits/DA/DA_fit_real.rds")
###############################################################
DA.params.short <- c("mu_0","alpha", "beta", "tau_u","tau_w","sigma_e_0","delta_0","gamma")
DA.params.long <- c(
  "mu_0",
  "alpha",
  "beta[1]",
  "beta[2]",
  "beta[3]",
  "beta[4]",
  "beta[5]",
  "beta[6]",
  "beta[7]",
  "delta_0",
  "sigma_e_0",
  "gamma",
  "tau_u[1]",
  "tau_u[2]",
  "tau_u[3]",
  "tau_w[1]",
  "tau_w[2]",
  "tau_w[3]",
  "L_u[1,1]",
  "L_u[2,1]",
  "L_u[2,2]",
  "L_u[3,1]",
  "L_u[3,2]",
  "L_u[3,3]",
  "L_w[1,1]",
  "L_w[2,1]",
  "L_w[2,2]",
  "L_w[3,1]",
  "L_w[3,2]",
  "L_w[3,3]"
  )

###############################################################
## direct access - fit real data
###############################################################
# list of data for Stan
DA.standata <- within(list(),
                      {
                        N_obs<-nrow(data)
                        RT <- data$rt
                        N_choices <- length(unique(data$answer))
                        rctype <- data$rctype_meff
                        accuracy <- data$acc
                        group <- data$group_meff
                        group_v <- data$group_meff
                        
                        N_subj <- length(unique(data$subj))
                        N_item <- length(unique(data$item))
                        
                        subj <- as.integer(factor(data$subj))
                        item <- as.integer(factor(data$item))
                        
                        n_u <- 4  # no. subj reff 
                        n_w <- 3  # no. item reff 
                      }
)

# fit the model 
DA.fit.real <- stan(data=DA.standata,file="./models/direct_access.stan",iter=3000, chains=3)
                                                              
# Save the fit to file so we dont have to rerun it every time
# filename <- paste("./fits/DA_fit_real_", Sys.time(), ".rds", sep="")
# filename <- gsub(":","",gsub(" ","_",filename))
# saveRDS(DA.fit.real, filename)

# print the model (can be also used with DA.params.long)
summaryDA <- print(summary(DA.fit.real,pars=DA.params.short, probs=c(0.025,0.975))$summary)

###############################################################
## direct access  - plot thetas and P_b
###############################################################
# Extract thetas by condition and group
prob_sr_i <- rstan::extract(DA.fit.real)$prob_sr_i
prob_sr_c <- rstan::extract(DA.fit.real)$prob_sr_c
prob_or_i <- rstan::extract(DA.fit.real)$prob_or_i
prob_or_c <- rstan::extract(DA.fit.real)$prob_or_c

# dataframe of thetas
prb_df_theta <- data.frame(probability = c(prob_sr_i,prob_sr_c,prob_or_i,prob_or_c),
                     rctype =  c(rep("SR",length(prob_sr_i)),
                                 rep("SR",length(prob_sr_c)),
                                 rep("OR",length(prob_or_i)),
                                 rep("OR",length(prob_or_c))),
                     group = c(rep("IWA",length(prob_sr_i)),
                               rep("control",length(prob_sr_c)),
                               rep("IWA",length(prob_or_i)),
                               rep("control",length(prob_or_c))))
# rename factor
prb_df_theta$rctype <- factor(prb_df_theta$rctype, levels=c("SR","OR"), labels=c("SR","OR")) 
# summarize means by group and condition
theta <- ddply(prb_df_theta, .(group, rctype), summarize,
                       mean = mean(probability))

# plot thetas
#png(filename="figures/theta_DA_real.png", res=300, width=1800, height=1100) # use this to save it separately
theta_plot <- ggplot(prb_df_theta, aes(probability, fill = group)) +
  geom_density(adjust=2,alpha=.5) + facet_grid(. ~ rctype) + xlim(c(0,1)) + xlab("parameter estimate")+
  geom_vline(data=theta, aes(xintercept=mean, colour=group), linetype="dashed", size=0.8) + 
  scale_color_manual(values=c("grey","black")) +
  theme_bw() + scale_fill_brewer(palette="Greys") + labs(title=" a) Probability of correct initial retrieval")

## Extract probability of backtracking
pb_i <- rstan::extract(DA.fit.real)$P_b_i
pb_c <- rstan::extract(DA.fit.real)$P_b_c

## dataframe
p_b_df <- data.frame(probability = c(pb_i,pb_c),
                     group = c(rep("IWA",length(pb_i)),
                               rep("control",length(pb_c))))
# summarize mean by group
p_b_df_mean <- ddply(p_b_df, .(group), summarize,
                     mean = round(mean(probability),2))
# percentages
p_b_df_mean_p <- ddply(p_b_df, .(group), summarize,
                       mean = round(mean(probability),2)*100)
# plot prob of backtracking
backtrack <- ggplot(p_b_df, aes(probability, fill = group)) +
  geom_density(alpha=.7) + xlim(c(0,1)) + xlab("parameter estimate") +
  geom_vline(data=p_b_df_mean, aes(xintercept=mean, color=group), linetype="dashed", size=0.8) +  ylab(" ") +
  scale_color_manual(values=c("gray","black")) +
  theme_bw() + ggtitle("b) Probability of backtracking") + scale_fill_brewer(palette="Greys")

# Figure 11 (plot theta and backtracking together)
grid.arrange(theta_plot,backtrack, nrow=2)

# extract delta
delta_i <- rstan::extract(DA.fit.real)$delta_i
delta_c<- rstan::extract(DA.fit.real)$delta_c
# dataframe
delta_df <- data.frame(probability = c(delta_i,delta_c),
                       group = c(rep("IWA",length(delta_i)),
                                 rep("control",length(delta_c))))
#summarize by group
delta_df_mean <- ddply(delta_df, .(group), summarize,
                       mean = round(mean(probability),0))

# plot delta (Figure 11)
ggplot(delta_df, aes(probability, fill = group)) +
  geom_density(alpha=.5) + xlim(c(0,3000)) + ylab(" ") + xlab("Parameter estimate") +
  geom_vline(data=delta_df_mean, aes(xintercept=mean, color=group), linetype="dashed", size=0.8) + 
  scale_color_manual(values=c("gray","black")) +
  theme_bw() + ggtitle("Time needed for backtracking (in ms)") + scale_fill_brewer(palette="Greys")


# extract sigma
sigma_e_i <- rstan::extract(DA.fit.real)$sigma_e_i
sigma_e_c <- rstan::extract(DA.fit.real)$sigma_e_c

#dataframe
sigma_DA <- data.frame(probability = c(sigma_e_i,sigma_e_c),
                       group = c(rep("IWA",length(sigma_e_i)),
                                 rep("control",length(sigma_e_c))))
# summarize by group
sigma_DA_mean <- ddply(sigma_DA, .(group), summarize,
                       mean = mean(probability))

# plot sigma 
sigma_DA_plot <- ggplot(sigma_DA, aes(probability, fill = group)) +
  geom_density(adjust=2,alpha=.5) + xlab("Parameter estimate")+
  geom_vline(data=sigma_DA_mean, aes(xintercept=mean, colour=group), linetype="dashed", size=0.8) + ylab(" ") +
  scale_color_manual(values=c("grey","black")) +
  theme_bw() + scale_fill_brewer(palette="Greys") + labs(title=expression("Posterior distribution of"~sigma))

### extract mu
mu_i <- rstan::extract(DA.fit.real)$mu_i
mu_c <- rstan::extract(DA.fit.real)$mu_c

# dataframe
mu_df <- data.frame(probability = c(mu_i,mu_c),
                     group = c(rep("IWA",length(mu_i)),
                               rep("control",length(mu_c))))
# summarize mean by group 
mu_mean <- ddply(mu_df, .(group), summarize,
                 mean = round(mean(probability),2))

# plot mu 
mu_plot <- ggplot(mu_df, aes(probability, fill = group)) +
  geom_density(adjust=2,alpha=.5)  + xlab("Parameter estimate")+
  geom_vline(data=mu_mean, aes(xintercept=mean, colour=group), linetype="dashed", size=0.8) +  ylab(" ") +
  scale_color_manual(values=c("grey","black")) +
  theme_bw() + scale_fill_brewer(palette="Greys") + labs(title=expression("Posterior distribution of"~mu))

# plot mu and sigma together (Figure 12)
pdf("figures/mu_sigma_DA.pdf")
grid.arrange(mu_plot,sigma_DA_plot, nrow=2)
dev.off()

#############################
#############################
####### 6. PPCs #############
#############################
#############################

#############################
##### ACTIVATION-BASED ######
#############################
#############################
###### ACCURACY        
#############################
# extract posterior and convert to proper format 
ACT.ppc.acc <- rstan::extract(ACT.fit.real, pars = c("gen_acc","gen_rctype","gen_RT","gen_group")) %>%
  purrr::map( ~  purrr::array_tree(.x, 1)) %>% purrr::transpose() %>%
  purrr::map_dfr(~ as_tibble(.x), .id = "sim")

# summarize by mean acc
ACT.ppc.acc.df <- ddply(ACT.ppc.acc , .(sim,gen_rctype,gen_group), summarize,
                        acc = mean(gen_acc))
# convert to percentages
ACT.ppc.acc.df$acc <-  ACT.ppc.acc.df$acc*100

# labels
ACT.ppc.acc.df$rctype_lab <- factor(ifelse(ACT.ppc.acc.df$gen_rctype==-1,"SR", "OR"),levels=c("SR", "OR"))
ACT.ppc.acc.df$group_lab  <- factor(ifelse(ACT.ppc.acc.df$gen_group==-1,"Controls", "IWAs"),levels=c("Controls", "IWAs"))

# mean accuracy per group and condition
ACT.acc.mean <- ddply(ACT.ppc.acc.df , .(rctype_lab, group_lab), summarize, acc = mean(acc))

# plot PPCs for accuracy 
# Figure 8
ACT.acc <- ggplot(ACT.ppc.acc.df, aes(x = rctype_lab, y = acc)) +
  geom_violin() + facet_grid(group_lab ~ .) +
  geom_pointrange(data=acc, mapping=aes(x=rctype_lab, y=mean.acc, ymin=se.lower, ymax=se.upper), size=.3) + theme_bw() + 
  ylab("Accuracy (%)") + 
  xlab("Condition") +
  ggtitle("Activation-based race model") 


#############################
###### LISTENING TIMES 
#############################
# extract posterior and convert to proper format 
ACT.ppc.lt <- rstan::extract(ACT.fit.real, pars = c("gen_acc","gen_rctype","gen_RT","gen_group")) %>%
purrr::map( ~  purrr::array_tree(.x, 1)) %>% purrr::transpose() %>%
  purrr::map_dfr(~ as_tibble(.x), .id = "sim")

# summarize, mean LT
ACT.ppc.lt.df <- ddply(ACT.ppc.lt, .(sim,gen_rctype,gen_group, gen_acc), summarize, rt = mean(gen_RT))

# labels
ACT.ppc.lt.df$acc_lab    <- factor(ifelse(ACT.ppc.lt.df$gen_acc==0,"Incorrect", "Correct"),levels=c("Correct","Incorrect"))
ACT.ppc.lt.df$rctype_lab <- factor(ifelse(ACT.ppc.lt.df$gen_rctype==-1,"SR", "OR"),levels=c("SR", "OR"))
ACT.ppc.lt.df$group_lab  <- factor(ifelse(ACT.ppc.lt.df$gen_group==-1,"Controls", "IWAs"),levels=c("Controls", "IWAs"))

# mean LT per condition, group, and accuracy 
ACT.LT.mean <- ddply(ACT.ppc.lt.df , .(acc_lab, rctype_lab, group_lab), summarize, rt = mean(rt))

# plot PPC for LT
# Figure 9
ACT.LT.ppc <- ggplot(ACT.ppc.lt.df, aes(x = acc_lab, y = rt)) +
  geom_violin() + facet_grid(group_lab ~ rctype_lab) + ylim(0, 4000) +
  geom_pointrange(data=LT, mapping=aes(x=acc_lab, y=mean.rt, ymin=se.lower, ymax=se.upper), size=1, shape=95) +
  ylab("Listening times (ms)") + 
  xlab("Accuracy in picture selection") +
  ggtitle("Activation-based model") 


#############################
##### DIRECT ACCESS #########
#############################
#############################
###### ACCURACY     
#############################
DA.ppc.acc <- rstan::extract(DA.fit.real, pars = c("gen_acc","gen_rctype","gen_RT","gen_group")) %>%
  purrr::map( ~  purrr::array_tree(.x, 1)) %>% purrr::transpose() %>%
  purrr::map_dfr(~ as_tibble(.x), .id = "sim")

# mean LT per condition, group 
DA.ppc.acc.df <- ddply(DA.ppc.acc , .(sim,gen_rctype,gen_group), summarize,
                       acc = mean(gen_acc))
DA.ppc.acc.df$acc <-  DA.ppc.acc.df$acc*100 # percentages

# labels
DA.ppc.acc.df$rctype_lab <- factor(ifelse(DA.ppc.acc.df$gen_rctype==-1,"SR", "OR"),levels=c("SR", "OR"))
DA.ppc.acc.df$group_lab  <- factor(ifelse(DA.ppc.acc.df$gen_group==-1,"Controls", "IWAs"),levels=c("Controls", "IWAs"))

# plot PPCs for accuracy, Figure 13
DA.acc <- ggplot(DA.ppc.acc.df, aes(x = rctype_lab, y = acc)) +
geom_violin() + facet_grid(group_lab ~ .) +
  geom_pointrange(data=acc, mapping=aes(x=rctype_lab, y=mean.acc, ymin=se.lower, ymax=se.upper), size=.3) + theme_bw() + 
  ylab("Accuracy (%)") + 
  xlab("Condition") +
  ggtitle("Direct-access model") 

#############################
###### LISTENING TIMES 
#############################
# extract posterior and convert to proper format 
DA.ppc.lt <- rstan::extract(DA.fit.real, pars = c("gen_acc","gen_rctype","gen_RT","gen_group")) %>%
  purrr::map( ~  purrr::array_tree(.x, 1)) %>% purrr::transpose() %>%
  purrr::map_dfr(~ as_tibble(.x), .id = "sim")

# mean LT per condition, group, and accuracy 
DA.ppc.lt.df <- ddply(DA.ppc.lt , .(sim,gen_rctype,gen_group, gen_acc), summarize, rt = mean(gen_RT))

# labels
DA.ppc.lt.df$acc_lab    <- factor(ifelse(DA.ppc.lt.df$gen_acc==0,"Incorrect", "Correct"),levels=c("Correct","Incorrect"))
DA.ppc.lt.df$rctype_lab <- factor(ifelse(DA.ppc.lt.df$gen_rctype==-1,"SR", "OR"),levels=c("SR", "OR"))
DA.ppc.lt.df$group_lab  <- factor(ifelse(DA.ppc.lt.df$gen_group==-1,"Controls", "IWAs"),levels=c("Controls", "IWAs"))

# plot PPC for LT, Figure 14 
DA.LT.ppc <- ggplot(DA.ppc.lt.df, aes(x = acc_lab, y = rt)) +
  geom_violin() + facet_grid(group_lab ~ rctype_lab) + ylim(0, 4000) +
  geom_pointrange(data=LT, mapping=aes(x=acc_lab, y=mean.rt, ymin=se.lower, ymax=se.upper), size=1, shape=95) +
  ylab("Listening times (ms)") + 
  xlab("Accuracy in picture selection") +
  ggtitle("Direct-access model") 


# pdf("figures/PPC_acc.pdf", width=8, height=6)
# grid.arrange(ACT.acc,DA.acc, nrow=1)
# dev.off()
# 
# pdf("figures/PPC_LT.pdf", width=10, height=7)
# grid.arrange(ACT.LT.ppc,DA.LT.ppc, nrow=1)
# dev.off()

###############################################################
## 7. Crossvalidation 
###############################################################

## prepare dataset for crossvalidation
data <- as.data.frame(data)
row.names(data) <- 1:dim(data)[1]
data$row <- row.names(data)

K <- 10 #folds
d <- data
G <- list()
for (i in 1:K) {
  G[[i]] <- sample_frac(group_by(d, subj),
                        (1/(K + 1 - i)))
  G[[i]]$k <- i
  d <<- anti_join(d, G[[i]],
                  by = c("subj", "item",
                         "rctype", "rt","group_meff","answer","acc"))
}
head(G[[1]]) 


dK <- bind_rows(G)
dK <- dK[order(dK$row), ]


############################
## Activation-based model ## 
############################
ACTdata <- plyr::llply(1:K, function(i) {
  list(
    RT = dK$rt,
    N_obs = nrow(dK),
    N_choices = length(unique(dK$answer)),
    winner = dK$answer,
    rctype = dK$rctype_meff,
    accuracy = dK$acc,
    group = dK$group_meff,
    
    N_subj = length(unique(dK$subj)),
    N_item = length(unique(dK$item)),
    
    subj = as.integer(factor(dK$subj)),
    item = as.integer(factor(dK$item)),
    
    n_u =  4, # no. subj reff
    n_w = 4, # no. item reff
    
    heldout = ifelse(dK$k == i, 1, 0))
})

str(ACTdata)
pointwise_ACT <- list()
for(i in 1:10){
  hnxval.dat<-ACTdata[[i]]
  ACTxval <- stan(file = "./models/activation_based.stan",
                  data = hnxval.dat,
                  iter = 6000, chains = 3,
                  refresh=0)
  
  ACTxval<-rstan::extract(ACTxval,pars="log_lik") 
  loglik<-ACTxval$log_lik 
  hldout<-which(hnxval.dat$heldout==1) 
  logmeans<-rep(NA,length(hldout))
  for(j in 1:length(hldout)){
    logmeans[j]<-log(mean(exp(loglik[,hldout[j]]))) 
  }
  pointwise_ACT[[i]]<-logmeans 
}
saveRDS(pointwise_ACT, "./crossvalidation/pointwise_ACT.rds")

## k-fold pointwise log mean likelihood:
pointwiseACT_flat<-Reduce(c,pointwise_ACT) 
## sum it up: 
(elpd_ACT<-sum(pointwiseACT_flat))
# -12514.99 
## compute SE using formula from Vehtari paper:
(elpd_ACT_se<-sqrt(length(pointwiseACT_flat)*
                     var(pointwiseACT_flat)))
# 49

###########################
## Direct-access model   ## 
###########################
DAdata <- plyr::llply(1:K,function(i) {
  list(
    RT = dK$rt,
    N_obs = nrow(dK),
    N_choices = length(unique(dK$answer)),
    winner = dK$answer,
    rctype = dK$rctype_meff,
    accuracy = dK$acc,
    group = dK$group_meff,
    
    N_subj = length(unique(dK$subj)),
    N_item = length(unique(dK$item)),
    
    subj = as.integer(factor(dK$subj)),
    item = as.integer(factor(dK$item)),
    
    n_u =  4, # no. subj reff
    n_w = 3, # no. item reff
    
    heldout = ifelse(dK$k == i, 1, 0))
})

pointwise_DA <- list()
for(i in 1:10){
  hnxval.dat<-DAdata[[i]]
  DAxval <- stan(file="./models/direct_access.stan",
                 data = hnxval.dat,
                 iter = 7000, chains = 3, refresh=0, 
                 control = list(adapt_delta = 0.9))
  DAxval<-rstan::extract(DAxval,pars="log_lik")
  loglik<-DAxval$log_lik
  hldout<-which(hnxval.dat$heldout==1)
  logmeans<-rep(NA,length(hldout))
  for(j in 1:length(hldout)){
    logmeans[j]<-log(mean(exp(loglik[,hldout[j]]))) 
  }
  pointwise_DA[[i]]<-logmeans 
}

saveRDS(pointwise_DA, "./crossvalidation/pointwise_DA.rds")
## flatten the k-fold pointwise log mean likelihood:
pointwiseDA_flat<-Reduce(c,pointwise_DA) 
## sum it up: 
(elpd_DA<-sum(pointwiseDA_flat))
# -12630 

## compute SE using formula from Vehtari paper:
(elpd_DA_se<-sqrt(length(pointwiseDA_flat)*
                    var(pointwiseDA_flat)))
# 51 

# Difference in elpd:
elpd_ACT-elpd_DA
# 115 

# SE of the difference:
SE.ACT.DA.diff  <- round(sqrt(length(pointwiseACT_flat)*var(pointwiseACT_flat-pointwiseDA_flat)))
# 69

# Plot elpd differences
# add a column with the elpd of each model to the dataset, also the differences
data$diff_elpd <- pointwiseACT_flat - pointwiseDA_flat 
data$elpd_act <- pointwiseACT_flat
data$elpd_da <- pointwiseDA_flat

# Label for the y axis:
y_axis_da_ab <- expression(widehat(elpd)[activation-based~~model] - 
                             widehat(elpd)[direct-acces~~model])

# Figure 15
elpd_plot <- ggplot(data,aes(x = rt, y = diff_elpd)) + 
  facet_grid(group_lab ~ acc_lab) +
  scale_x_continuous(name= "log-scaled LT (ms)", 
                     trans = log_trans(), 
                     breaks=c(200,500,1000,2000,3000,4000,5000,7000), 
                     labels=c(200,500,1000,2000,3000,4000,5000,7000)) + 
  theme(axis.text.x = 
          element_text(angle = 90, vjust = 1,hjust = 1),
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8)) + 
  scale_y_continuous(name = y_axis_da_ab, breaks = seq(-7,7,1)) +
  geom_hex(bins = 8) + 
  scale_fill_gradientn(colours = c("grey","black"), 
                       name = "Number of\nobservations") + 
  facet_grid(group_lab~acc_lab)  + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  ggtitle("Comparison of models")

# compute elpd values and corresponding SE by condition and group:
elpd_OR_controls <- filter(data, group_lab=="Controls"&rctype=="OR")
elpd_OR_control_DA <- sum(elpd_OR_controls$elpd_da)
elpd_OR_control_ACT <- sum(elpd_OR_controls$elpd_act)
diff.ACT.DA.controls.OR <- round(elpd_OR_control_ACT-elpd_OR_control_DA)
SE.ACT.DA.controls.OR  <- round(sqrt(length(elpd_OR_controls$elpd_act)*
                                       var(elpd_OR_controls$elpd_act-elpd_OR_controls$elpd_da)))

elpd_SR_controls <- filter(data, group_lab=="Controls"&rctype=="SR")
elpd_SR_control_DA <- sum(elpd_SR_controls$elpd_da)
elpd_SR_control_ACT <- sum(elpd_SR_controls$elpd_act)
diff.ACT.DA.controls.SR <- round(elpd_SR_control_ACT-elpd_SR_control_DA)
SE.ACT.DA.controls.SR  <- round(sqrt(length(elpd_SR_controls$elpd_act)*
                                       var(elpd_SR_controls$elpd_act-elpd_SR_controls$elpd_da)))

elpd_OR_IWA <- filter(data, group_lab=="IWAs"&rctype=="OR")
elpd_OR_IWA_DA <- sum(elpd_OR_IWA$elpd_da)
elpd_OR_IWA_ACT <- sum(elpd_OR_IWA$elpd_act)
diff.ACT.DA.IWA.OR <- round(elpd_OR_IWA_ACT-elpd_OR_IWA_DA)
SE.ACT.DA.IWA.OR <- round(sqrt(length(elpd_OR_IWA$elpd_act)*
                                 var(elpd_OR_IWA$elpd_act-elpd_OR_IWA$elpd_da)))

elpd_SR_IWA <- filter(data, group_lab=="IWAs"&rctype=="SR")
elpd_SR_IWA_DA <- sum(elpd_SR_IWA$elpd_da)
elpd_SR_IWA_ACT <- sum(elpd_SR_IWA$elpd_act)
diff.ACT.DA.IWA.SR <- round(elpd_SR_IWA_ACT-elpd_SR_IWA_DA)
SE.ACT.DA.IWA.SR <- round(sqrt(length(elpd_SR_IWA$elpd_act)*
                                 var(elpd_SR_IWA$elpd_act-elpd_SR_IWA$elpd_da)))

elpds <- t(data.frame(SR.controls = diff.ACT.DA.controls.SR,
                      OR.controls = diff.ACT.DA.controls.OR,
                      SR.IWA = diff.ACT.DA.IWA.SR,
                      OR.IWA = diff.ACT.DA.IWA.OR))
SE <- t(data.frame(SR.controls = SE.ACT.DA.controls.SR,
                   OR.controls = SE.ACT.DA.controls.OR,
                   SR.IWA = SE.ACT.DA.IWA.SR,
                   OR.IWA = SE.ACT.DA.IWA.OR))

# dataframe that gathers all elpds and SE
elpd_table <- cbind(elpds,SE)
colnames(elpd_table) <- c("$\\widehat{elpd}_{diff}$", "SE")

# print table with elpds and SE
# Table 1
apa_table(elpd_table
          , caption = "$\\widehat{Elpd}$ differences between the activation-based and the direct-access model across conditions and groups"
          , escape=FALSE
          , digits=0
)

