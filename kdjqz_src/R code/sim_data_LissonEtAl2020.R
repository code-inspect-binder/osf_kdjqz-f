###############################################################
##### Paula Lissón - lissonh@gmail.com ########################
############################################################### 
##### This script is used to generate simulated data for     ##
#####  the recovery of the parameters in Lisson et al. 2020  ##
########## Outline of the script: #############################
##### - 0. Load packages
##### - 1. Activation-based model: extract parameters
#####   1.1. Generate data from posterior means of the model
#####   1.2. Fit simulated data
#####   1.3. Recovery of the parameters (plots)
##### - 2. Direct-access model: extract parameters
#####   2.1. Generate data from posterior means of the model
#####   2.2. Fit simulated data
#####   2.3. Recovery of the parameters (plots)
 
###############################################################
## 0. Load packages -----
###############################################################
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(ggplot2)
library(reshape2)
library(MASS)
library(gtools)
library(loo)
library(scales)

# user-made functions to extract posteriors from stan fit and get them to proper format:
extract_mean_posteriors <- function(data, params){
  df <- data.frame(matrix(ncol = length(params), nrow = 1))
  colnames(df) <- params
  for(p in params){
    df[1,p] <- mean(subset(data, variable==p)$value)
  }
  return(df)
}

super_melt <- function(stanfit, params){
  df <- as.data.frame(stanfit)
  df <- df[,params]
  df <- melt(df)
  return(df)
}

# load the models previously fitted:
# 
# real data 
# ACT.fit.real <- readRDS("./fits/ACT/ACT_fit_real.rds")
# DA.fit.real <- readRDS("./fits/DA/DA_fit_real.rds")
# 
# simulated data 
# ACT.fit.sim <- readRDS("./fits/ACT/ACT_fit_sim.rds")
# ACT.fit.sim <- readRDS("./fits/DA/DA_fit_sim.rds")

###############################################################
###############################################################
##  1. ACTIVATION-BASED MODEL: extract parameters--------------
###############################################################
###############################################################
# parameters of the model:
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
  "sigma_e", 
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
  "L_w[4,4]"
)


###############################################################
## 1.1. Activation-based: generate data from posterior means ----
###############################################################
# Create template. The fields 'acc_meff', 'answer'_meff and 'rt' will be populated later
ACT.data.sim <- data.frame(rctype = rep(c(rep(-1, 1000), rep(1, 1000)),1),
                           item = rep(1:10, 200),
                           subj = rep(1:50, each=10, times=4),
                           rt = rep(0, 2000),
                           answer = rep(0,2000),
                           group = rep(c(rep(-1, 10), rep(1, 10)), 100))


ACT.data.sim.nsubj<- length(unique(ACT.data.sim$subj))
ACT.data.sim.nitem<- length(unique(ACT.data.sim$item))
ACT.data.sim.ntotal<-nrow(ACT.data.sim)

# mean posteriors (always use params.long here)
ACT.mp <- extract_mean_posteriors(super_melt(ACT.fit.real, ACT.params.long), ACT.params.long)

#################################
### Random effects by subject ###
################################# 
# Cholesky 
L_u <- matrix(c(ACT.mp$`L_u[1,1]`,ACT.mp$`L_u[2,1]`,ACT.mp$`L_u[3,1]`,ACT.mp$`L_u[4,1]`,
                0,ACT.mp$`L_u[2,2]`,ACT.mp$`L_u[3,2]`,ACT.mp$`L_u[4,2]`,
                0,0,ACT.mp$`L_u[3,3]`,ACT.mp$`L_u[4,3]`,
                0,0,0,ACT.mp$`L_u[4,4]`
),nrow=4,ncol=4)

## random uncorrelated variables 
z_u1 <- rnorm(ACT.data.sim.nsubj, 0, .5)
z_u2 <- rnorm(ACT.data.sim.nsubj, 0, .5)
z_u3 <- rnorm(ACT.data.sim.nsubj, 0, .5)
z_u4 <- rnorm(ACT.data.sim.nsubj, 0, .5)
z_u <- matrix(c(z_u1, z_u2, z_u3,z_u4), ncol = ACT.data.sim.nsubj, byrow = T)

# Diagonal matrix Tau
Tau_u <- matrix(c(ACT.mp$`tau_u[1]`,0,0,0,
                  0,ACT.mp$`tau_u[2]`,0,0,
                  0,0,ACT.mp$`tau_u[3]`,0,
                  0,0,0,ACT.mp$`tau_u[4]`
),nrow=4,ncol=4)

# final matrix of correlated adjustments for subject
ACT.u <- Tau_u %*% L_u %*% z_u %>% t

#############################
### Random effects by item ###
##############################
# cholesky
L_w <- matrix(c(ACT.mp$`L_w[1,1]`,ACT.mp$`L_w[2,1]`,ACT.mp$`L_w[3,1]`,ACT.mp$`L_w[4,1]`,
                0,ACT.mp$`L_w[2,2]`,ACT.mp$`L_w[3,2]`,ACT.mp$`L_w[4,2]`,
                0,0,ACT.mp$`L_w[3,3]`,ACT.mp$`L_w[4,3]`,
                0,0,0,ACT.mp$`L_w[4,4]`
),nrow=4,ncol=4)

## random uncorrelated variables 
z_w1 <- rnorm(ACT.data.sim.nitem, 0, .5)
z_w2 <- rnorm(ACT.data.sim.nitem, 0, .5)
z_w3 <- rnorm(ACT.data.sim.nitem, 0, .5)
z_w4 <- rnorm(ACT.data.sim.nitem, 0, .5)
z_w <- matrix(c(z_w1, z_w2, z_w3,z_w4), ncol = ACT.data.sim.nitem, byrow = T)

# Diagonal matrix Tau
Tau_w <- matrix(c(ACT.mp$`tau_w[1]`,0,0,0,
                  0,ACT.mp$`tau_w[2]`,0,0,
                  0,0,ACT.mp$`tau_w[3]`,0,
                  0,0,0,ACT.mp$`tau_w[4]`
),nrow=4,ncol=4)

# final matrix of correlated adjustments 
ACT.w <- Tau_w %*% L_w %*% z_w %>% t

for(n in 1:ACT.data.sim.ntotal){
  # sample from both, choose shortest time
  sub<- ACT.data.sim$subj[n]
  item<-ACT.data.sim$item[n]
  grp<- ACT.data.sim$group[n]
  rct<- ACT.data.sim$rctype[n]
  
  mu_1 <- ACT.mp$`alpha[1]`+ACT.u[sub,1]+ACT.w[item,1]+ grp*(ACT.mp$`beta[1]`+ACT.w[item,3]) + rct*(ACT.mp$`beta[3]`+ACT.u[sub,3]) + grp*rct*ACT.mp$`beta[5]`
  mu_2 <- ACT.mp$`alpha[2]`+ACT.u[sub,2]+ACT.w[item,2]+ grp*(ACT.mp$`beta[2]`+ACT.w[item,4]) + rct*(ACT.mp$`beta[4]`+ACT.u[sub,4]) + grp*rct*ACT.mp$`beta[6]`
  
  sig <- ACT.mp$sigma_e + grp*ACT.mp$`beta[7]`
  rt1 <- rlnorm(1, mu_1, sig);
  rt2 <- rlnorm(1, mu_2, sig); 
  
  ACT.data.sim$answer[n]  <- ifelse(rt1 < rt2, 1, 2)
  ACT.data.sim$rt[n]      <- ifelse(rt1 < rt2, rt1, rt2)
  ACT.data.sim$acc_meff[n]<- ifelse((rt1 < rt2 & rct == -1)|(rt1 > rt2 & rct == 1), 1, 0)
}

###############################################################
## 1.2. Activation-based: fit simulated data     -------------
###############################################################
# data for Stan:
ACT.standata.sim <- within(list(),
                           {
                             N_obs<-nrow(ACT.data.sim)
                             RT <- ACT.data.sim$rt
                             N_choices <- length(unique(ACT.data.sim$answer))
                             winner <- ACT.data.sim$answer
                             rctype <- ACT.data.sim$rctype
                             group <- ACT.data.sim$group
                             
                             subj <- as.integer(factor(ACT.data.sim$subj))
                             item <- as.integer(factor(ACT.data.sim$item))
                             
                             N_subj <- length(unique(ACT.data.sim$subj))
                             N_item <- length(unique(ACT.data.sim$item))
                             
                             n_u <- 4 # no. subj reff
                             n_w <- 4 # no. item reff
                           }
)

# fit simulated data
ACT.fit.sim <- stan(data=ACT.standata.sim,chains=3,
                    file="./models/activation_based.stan", iter=6000)

# print summary                    
print(summary(ACT.fit.sim,pars=ACT.params.short, probs=c(0.025,0.975))$summary)

# Save the fit to file so we dont have to rerun it every time
# filename <- paste("./fits/act_fit_sim_", Sys.time(), ".rds", sep="")
# filename <- gsub(":","",gsub(" ","_",filename))
# saveRDS(ACT.fit.sim, filename)

###############################################################
## 1.3. Recovery of the parameters (plots) ------
## Activation - compare posteriors from simulated data with
##              posteriors from real data
###############################################################
ACT.fit.summary.sim <- summary(ACT.fit.sim,pars=ACT.params.long, probs=c(0.025,0.975))$summary
# arrange dataframe for plot
ACT.comparison.df <- data.frame(par = rownames(ACT.fit.summary.sim),
                                mean = ACT.fit.summary.sim[, "mean"],
                                p025 = ACT.fit.summary.sim[, "2.5%"],
                                p975 = ACT.fit.summary.sim[, "97.5%"],
                                # to get from key/value dataframe to proper format
                                gen = as.data.frame(t(ACT.mp[,]))$V1)

ACT.comparison.df$par <-   with(ACT.comparison.df, factor(par, rev(par)))
ACT.comparison.df$lower <- with(ACT.comparison.df, p025 - gen)
ACT.comparison.df$middle <-with(ACT.comparison.df, mean - gen)
ACT.comparison.df$upper <- with(ACT.comparison.df, p975 - gen)

# Figure A1
# Plot the discrepancies
#png(filename="figures/discrepanciesACT.png", res=300, width=1500, height=800)
ggplot(ACT.comparison.df[1:12,]) +
  aes(x = par, y = middle, ymin = lower, ymax = upper) +
  scale_x_discrete() +
  labs(y = "Discrepancy", x = NULL) +
  geom_abline(intercept = 0, slope = 0, color = "darkgrey", lwd=1) +
  geom_linerange(lwd=1.0) +
  geom_point(size = 2.3) +
  coord_flip() +
  theme_light()


###############################################################
###############################################################
##    2. DIRECT ACCESS MODEL: extract parameters ---------------                
###############################################################
###############################################################
# parameters of the model:
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
  "tau_u[4]",
  "tau_w[1]",
  "tau_w[2]",
  "tau_w[3]",
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
  "L_w[3,3]"
)

###############################################################
## 2.1. Direct access  - generate data from posterior means----
###############################################################
# Simulation data
DA.data.sim <- data.frame(rctype = rep(c(rep(-1, 1000), rep(1, 1000)),2),
                          items = rep(1:10, 400),
                          subjects = rep(1:50, each=10, times=8),
                          rt = rep(0, 4000),
                          acc = rep(0,4000),
                          group = rep(c(rep(-1, 10), rep(1, 10)), 200))

DA.data.sim.nsubj<- length(unique(DA.data.sim$subjects))
DA.data.sim.nitem<- length(unique(DA.data.sim$items))
DA.data.sim.ntotal<-nrow(DA.data.sim)

# mean posteriors
DA.mp <- extract_mean_posteriors(super_melt(DA.fit.real,DA.params.long), DA.params.long)

#################################
### Random effects by subject ###
#################################
# Correlation matrix C
L_u <- matrix(c(DA.mp$`L_u[1,1]`,DA.mp$`L_u[2,1]`,DA.mp$`L_u[3,1]`,DA.mp$`L_u[4,1]`,
                0,DA.mp$`L_u[2,2]`,DA.mp$`L_u[3,2]`,DA.mp$`L_u[4,2]`,
                0,0,DA.mp$`L_u[3,3]`,DA.mp$`L_u[4,3]`,
                0,0,0,DA.mp$`L_u[4,4]`
),nrow=4,ncol=4)


## random uncorrelated variables 
z_u1 <- rnorm(DA.data.sim.nsubj, 0, .5)
z_u2 <- rnorm(DA.data.sim.nsubj, 0, .5)
z_u3 <- rnorm(DA.data.sim.nsubj, 0, .5)
z_u4 <- rnorm(DA.data.sim.nsubj, 0, .5)
z_u <- matrix(c(z_u1, z_u2, z_u3, z_u4), ncol = DA.data.sim.nsubj, byrow = T)

# Diagonal matrix Tau
Tau_u <- matrix(c(DA.mp$`tau_u[1]`,0,0,0,
                  0,DA.mp$`tau_u[2]`,0,0,
                  0,0,DA.mp$`tau_u[3]`,0,
                  0,0,0,DA.mp$`tau_u[4]`
),nrow=4,ncol=4)

# final matrix of correlated adjustments 
DA.u <- Tau_u %*% L_u %*% z_u %>% t

##############################
### Random effects by item ###
##############################
##############################
## Cholesky factor 
L_w <- matrix(c(DA.mp$`L_w[1,1]`,DA.mp$`L_w[2,1]`,DA.mp$`L_w[3,1]`,
                0,DA.mp$`L_w[2,2]`,DA.mp$`L_w[3,2]`,
                0,0,DA.mp$`L_w[3,3]`
),nrow=3,ncol=3)

## random uncorrelated variables 
z_w1 <- rnorm(DA.data.sim.nitem, 0, .5)
z_w2 <- rnorm(DA.data.sim.nitem, 0, .5)
z_w3 <- rnorm(DA.data.sim.nitem, 0, .5)
z_w <- matrix(c(z_w1, z_w2, z_w3), ncol = DA.data.sim.nitem, byrow = T)

# Diagonal matrix Tau used to scale the z_w
Tau_w <- matrix(c(DA.mp$`tau_w[1]`,0,0,
                  0,DA.mp$`tau_w[2]`,0,
                  0,0,DA.mp$`tau_w[3]`
),nrow=3,ncol=3)

# final correlated adjustments 
DA.w <- Tau_w %*% L_w %*% z_w %>% t

for(i in 1:DA.data.sim.ntotal){
  mu <- DA.mp$mu_0 + DA.u[DA.data.sim$subjects[i],1] + DA.w[DA.data.sim$items[i],1] + DA.data.sim$group[i]*DA.mp$`beta[1]`;
  theta <- inv.logit(DA.mp$alpha + DA.u[DA.data.sim$subjects[i],2] + DA.w[DA.data.sim$items[i],2] + DA.data.sim$rctype[i]*(DA.mp$`beta[2]`+ DA.u[DA.data.sim$subjects[i],3]) + 
                       DA.data.sim$group[i]*(DA.mp$`beta[3]`+DA.w[DA.data.sim$items[i],3])+ DA.data.sim$rctype[i]*DA.data.sim$group[i]*DA.mp$`beta[4]`);
  P_b <- inv.logit(DA.mp$gamma + DA.u[DA.data.sim$subjects[i],4] +DA.data.sim$group[i]*DA.mp$`beta[5]`);
  delta <- DA.mp$delta_0 + DA.data.sim$group[i]*DA.mp$`beta[6]`;
  sigma_e <- DA.mp$sigma_e_0 + DA.data.sim$group[i]*DA.mp$`beta[7]`;

  # original retrieval
  o_r <- rbinom(1,1,theta)
  if(o_r){
    DA.data.sim$acc[i] <- 1
    DA.data.sim$rt[i] <- rlnorm(1, mu, sigma_e);
  }
  else {
    acc <- rbinom(1,1,P_b)
    DA.data.sim$acc[i] <- acc
    if(acc){
      DA.data.sim$rt[i] <- rlnorm(1, mu + delta, sigma_e);
    }
    else {
      DA.data.sim$rt[i] <- rlnorm(1, mu, sigma_e);
    }
  }
}


###############################################################
## 2.2. Direct access  - fit simulated data -------------------
###############################################################
DA.standata.sim <- within(list(),
                          {
                            N_obs<-nrow(DA.data.sim)
                            RT <- DA.data.sim$rt
                            accuracy <- DA.data.sim$acc
                            rctype <- DA.data.sim$rctype
                            group <- DA.data.sim$group
                            
                            N_subj <- length(unique(DA.data.sim$subjects))
                            N_item <- length(unique(DA.data.sim$items))
                            
                            subj <- as.integer(factor(DA.data.sim$subjects))
                            item <- as.integer(factor(DA.data.sim$items))

                            n_u <- 4  # no. subj reff 
                            n_w <- 3  # no. item reff 
                          }
                          
)

DA.fit.sim <- stan(data=DA.standata.sim,
                   file="./models/direct_access.stan", 
                   iter=7000,chains=3)

# Save the fit to file so we dont have to rerun it every time
# filename <- paste("./fits/da_fit_sim_", Sys.time(), ".rds", sep="")
# filename <- gsub(":","",gsub(" ","_",filename))
# saveRDS(DA.fit.sim, filename)

###############################################################
## 2.3. Recovery of the parameters (plots) --------------------
## Direct access  - compare posteriors from simulated data with
##                  posteriors from real data
###############################################################

DA.fit.summary.sim <- summary(DA.fit.sim,pars=DA.params.long, probs=c(0.025,0.975))$summary
saveRDS(DA.fit.summary.sim, "./fits/DA.fit.summary.sim.rds")

DA.comparison.df <- data.frame(par = rownames(DA.fit.summary.sim),
                               mean = DA.fit.summary.sim[, "mean"],
                               p025 = DA.fit.summary.sim[, "2.5%"],
                               p975 = DA.fit.summary.sim[, "97.5%"],
                               # to get from key/value dataframe to proper format
                               gen = as.data.frame(t(DA.mp[,]))$V1)

DA.comparison.df$par <-   with(DA.comparison.df, factor(par, rev(par)))
DA.comparison.df$lower <- with(DA.comparison.df, p025 - gen)
DA.comparison.df$middle <-with(DA.comparison.df, mean - gen)
DA.comparison.df$upper <- with(DA.comparison.df, p975 - gen)

# Figure A2
# Plot the discrepancies
ggplot(DA.comparison.df[1:14,]) +
  aes(x = par, y = middle, ymin = lower, ymax = upper) +
  scale_x_discrete() +
  labs(y = "Discrepancy", x = NULL) +
  geom_abline(intercept = 0, slope = 0, color = "darkgrey", lwd=1) +
  geom_linerange(lwd=1.0) +
  geom_point(size = 2.3) +
  coord_flip() +
  theme_light()


