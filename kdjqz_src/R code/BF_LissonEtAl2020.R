###############################################################
##### Paula Lissón - lissonh@gmail.com ########################
###############################################################
#### Code to perform the BF analysis of Lisson et al. 2020  ###
###############################################################
########## Outline of the script:
##### - Load packages
##### - 0. Preprocessing of the data
##### - 1. Activation-based model
#####   1.1. Mu 
#####   1.1.1. BF and plots for mu
#####   1.2. Sigma
##### - 2. Direct-access model
#####   2.1. Mu
#####   2.2. Theta
#####   2.2.1. BF and plots for theta
#####   2.3. Pb 
#####   2.4. Delta
#####   2.5. Sigma
 
## 0. Load packages ------------------------------------------------
library(rstan)
options(mc.cores = parallel::detectCores())
library(plyr)
library(tidyverse)
library(bridgesampling)
library(bayesplot)
library(purrr)
options(scipen=999)
set.seed(4)

# priors sd used for the sensitivity analysis
priors <- c(0.1,0.3, 0.5)


## 0. Preprocessing of the data ---------------------------------
data <- read.table("data/segLT_AMW.csv", sep=",", header = TRUE, stringsAsFactors = FALSE) %>% 
  select(grp, subj, sent, seg, item, acc, rt, segtype, contr) %>% # relevant columns
  filter(segtype!= "final") %>%  # rows with listening times
  filter(contr== "SOSS") %>% # relevant contrast
  filter(item>=11) %>% # relevant items 
  filter(seg=="c" | seg=="d") # relevant segments (NP2-VP1)

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


#### 1.Activation-based model --------------------------------------------------
# create list of data for Stan
# this model has 2 reff for items, so only needed when there's no beta for group: 
ACT.data.mu0 <- within(list(),
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
                     n_w <- 2 # no. item reff
                   }
)

# this is the list needed for all the other activation-based models:
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

##### 1.1. Mu -----------------------------
# No group adjustment for mu (mu0)
ACT.BF.mu0 <- stan(data=ACT.data.mu0,file="./BayesFactor/models/ACT/ACT.BF.mu0.stan",
                     iter=40000, warmup =1000, chains=3,control = list(adapt_delta = 0.9))
ACT.BF.mu0.bridge <- bridge_sampler(ACT.BF.mu0)
saveRDS(ACT.BF.mu0.bridge, "./BayesFactor/marginal_lik/ACT.BF.mu0.rds")

# prior sd of the beta for group 01
# and 0,0.1 for the beta interaction
ACT.BF.mu.group01.int01 <- stan(data=ACT.data,file="./BayesFactor/models/ACT/ACT.BF.mu.group01.int01.stan",
                            iter=40000, warmup =1000, chains=3,control = list(adapt_delta = 0.9))
ACT.BF.mu.group01.int01.bridge <- bridge_sampler(ACT.BF.mu.group01.int01)
saveRDS(ACT.BF.mu.group01.int01.bridge, "./BayesFactor/marginal_lik/ACT/ACT.BF.mu.group01.int01.rds")

#  and 0,0.3 for the beta interaction
ACT.BF.mu.group01.int03 <- stan(data=ACT.data,file="./BayesFactor/models/ACT/ACT.BF.mu.group01.int03.stan",
                                iter=40000, warmup =1000, chains=3,control = list(adapt_delta = 0.9))
ACT.BF.mu.group01.int03.bridge <- bridge_sampler(ACT.BF.mu.group01.int03)
saveRDS(ACT.BF.mu.group01.int03.bridge, "./BayesFactor/marginal_lik/ACT/ACT.BF.mu.group01.int03.rds")

# and 0,0.5 for the beta interaction
ACT.BF.mu.group01.int05 <- stan(data=ACT.data,file="./BayesFactor/models/ACT/ACT.BF.mu.group01.int05.stan",
                                iter=40000, warmup =1000, chains=3,control = list(adapt_delta = 0.9))
ACT.BF.mu.group01.int05.bridge <- bridge_sampler(ACT.BF.mu.group01.int05)
saveRDS(ACT.BF.mu.group01.int05.bridge, "./BayesFactor/marginal_lik/ACT/ACT.BF.mu.group01.int05.rds")

# prior sd of the beta for group 03
# and 0,0.1 for the beta interaction
ACT.BF.mu.group03.int01 <- stan(data=ACT.data,file="./BayesFactor/models/ACT/ACT.BF.mu.group03.int01.stan",
                                iter=40000, warmup =1000, chains=3,control = list(adapt_delta = 0.9))
ACT.BF.mu.group03.int01.bridge <- bridge_sampler(ACT.BF.mu.group03.int01)
saveRDS(ACT.BF.mu.group03.int01.bridge, "./BayesFactor/marginal_lik/ACT/ACT.BF.mu.group03.int01.rds")

#  and 0,0.3 for the beta interaction
ACT.BF.mu.group03.int03 <- stan(data=ACT.data,file="./BayesFactor/models/ACT/ACT.BF.mu.group03.int03.stan",
                                iter=40000, warmup =1000, chains=3,control = list(adapt_delta = 0.9))
ACT.BF.mu.group03.int03.bridge <- bridge_sampler(ACT.BF.mu.group03.int03)
saveRDS(ACT.BF.mu.group03.int03.bridge, "./BayesFactor/marginal_lik/ACT/ACT.BF.mu.group03.int03.rds")

# and 0,0.5 for the beta interaction
ACT.BF.mu.group03.int05 <- stan(data=ACT.data,file="./BayesFactor/models/ACT/ACT.BF.mu.group03.int05.stan",
                                iter=40000, warmup =1000, chains=3,control = list(adapt_delta = 0.9))
ACT.BF.mu.group03.int05.bridge <- bridge_sampler(ACT.BF.mu.group03.int05)
saveRDS(ACT.BF.mu.group03.int05.bridge, "./BayesFactor/marginal_lik/ACT/ACT.BF.mu.group03.int05.rds")


# prior sd of the beta for group 05
# and 0,0.1 for the beta interaction
ACT.BF.mu.group05.int01 <- stan(data=ACT.data,file="./BayesFactor/models/ACT/ACT.BF.mu.group05.int01.stan",
                                iter=40000, warmup =1000, chains=3,control = list(adapt_delta = 0.9))
ACT.BF.mu.group05.int01.bridge <- bridge_sampler(ACT.BF.mu.group05.int01)
saveRDS(ACT.BF.mu.group05.int01.bridge, "./BayesFactor/marginal_lik/ACT/ACT.BF.mu.group05.int01.rds")

#  and 0,0.3 for the beta interaction
ACT.BF.mu.group05.int03 <- stan(data=ACT.data,file="./BayesFactor/models/ACT/ACT.BF.mu.group05.int03.stan",
                                iter=40000, warmup =1000, chains=3,control = list(adapt_delta = 0.9))
ACT.BF.mu.group05.int03.bridge <- bridge_sampler(ACT.BF.mu.group05.int03)
saveRDS(ACT.BF.mu.group05.int03.bridge, "./BayesFactor/marginal_lik/ACT/ACT.BF.mu.group05.int03.rds")

# and 0,0.5 for the beta interaction
ACT.BF.mu.group05.int05 <- stan(data=ACT.data,file="./BayesFactor/models/ACT/ACT.BF.mu.group05.int05.stan",
                                iter=40000, warmup =1000, chains=3,control = list(adapt_delta = 0.9))
ACT.BF.mu.group05.int05.bridge <- bridge_sampler(ACT.BF.mu.group05.int01)
saveRDS(ACT.BF.mu.group05.int05.bridge, "./BayesFactor/marginal_lik/ACT/ACT.BF.mu.group05.int05.rds")

#### 1.1.1. BF and plots ----
# load the bridge fits
ACT.BF.mu0 <- readRDS("./BayesFactor/marginal_lik/ACT/ACT.BF.mu0.rds")

ACT.BF.mu.group01.int01 <- readRDS("./BayesFactor/marginal_lik/ACT/ACT.BF.mu.group01.int01.rds")
ACT.BF.mu.group01.int03 <- readRDS("./BayesFactor/marginal_lik/ACT/ACT.BF.mu.group01.int03.rds")
ACT.BF.mu.group01.int05 <- readRDS("./BayesFactor/marginal_lik/ACT/ACT.BF.mu.group01.int05.rds")

ACT.BF.mu.group03.int01 <- readRDS("./BayesFactor/marginal_lik/ACT/ACT.BF.mu.group03.int01.rds")
ACT.BF.mu.group03.int03 <- readRDS("./BayesFactor/marginal_lik/ACT/ACT.BF.mu.group03.int03.rds")
ACT.BF.mu.group03.int05 <- readRDS("./BayesFactor/marginal_lik/ACT/ACT.BF.mu.group03.int05.rds")

ACT.BF.mu.group05.int01 <- readRDS("./BayesFactor/marginal_lik/ACT/ACT.BF.mu.group05.int01.rds")
ACT.BF.mu.group05.int03 <- readRDS("./BayesFactor/marginal_lik/ACT/ACT.BF.mu.group05.int03.rds")
ACT.BF.mu.group05.int05 <- readRDS("./BayesFactor/marginal_lik/ACT/ACT.BF.mu.group05.int05.rds")

BFmu0.g01.i01 <- (bayes_factor(ACT.BF.mu0, ACT.BF.mu.group01.int01))[[1]]
BFmu0.g01.i03 <- (bayes_factor(ACT.BF.mu0, ACT.BF.mu.group01.int03))[[1]]
BFmu0.g01.i05 <- (bayes_factor(ACT.BF.mu0, ACT.BF.mu.group01.int05))[[1]]

BFmu0.g03.i01 <- (bayes_factor(ACT.BF.mu0, ACT.BF.mu.group03.int01))[[1]]
BFmu0.g03.i03 <- (bayes_factor(ACT.BF.mu0, ACT.BF.mu.group03.int03))[[1]]
BFmu0.g03.i05 <- (bayes_factor(ACT.BF.mu0, ACT.BF.mu.group03.int05))[[1]]

BFmu0.g05.i01 <- (bayes_factor(ACT.BF.mu0, ACT.BF.mu.group05.int01))[[1]]
BFmu0.g05.i03 <- (bayes_factor(ACT.BF.mu0, ACT.BF.mu.group05.int03))[[1]]
BFmu0.g05.i05 <- (bayes_factor(ACT.BF.mu0, ACT.BF.mu.group05.int05))[[1]]


BF.ACT.mu <- data.frame(group = c(0.1,0.1,0.1, 0.3, 0.3, 0.3, 0.5, 0.5, 0.5),
                     interaction = rep(c(0.1,0.3,0.5),3), 
                              BF.raw = c(BFmu0.g01.i01, BFmu0.g01.i03, BFmu0.g01.i05,
                                     BFmu0.g03.i01, BFmu0.g03.i03, BFmu0.g03.i05,
                                     BFmu0.g05.i01, BFmu0.g05.i03, BFmu0.g05.i05))
                              
# 3D plot
# scatterplot3d(x=BF.ACT.mu$group, y=as.numeric(BF.ACT.mu$BF), z=BF.ACT.mu$interaction,
#               grid=TRUE, box=FALSE, color="black", pch = 16, 
#               xlab = "Prior SD for the group" ~beta, ylab="BF",
#               zlab="Prior SD for the interaction" ~beta,
#               title(main="BF with different prior SD for the group and the 
#           interaction adjustments in the parameter" ~mu))

##### 1.2. Sigma ------------------------------------
# no group adjustment for sigma 
ACT.BF.sigma0 <- stan(data=ACT.data,file="./BayesFactor/models/ACT/ACT.BF.sigma0.stan",
                               iter=40000, warmup =1000, chains=3,control = list(adapt_delta = 0.9))
ACT.BF.sigma0.bridge <- bridge_sampler(ACT.BF.sigma0) 
saveRDS(ACT.BF.sigma0.bridge, "./BayesFactor/marginal_lik/ACT/ACT.BF.sigma0.rds")

# group adjustment for sigma with prior 0,0.1 for the beta
ACT.BF.sigma.prior01 <- stan(data=ACT.data,file="./BayesFactor/models/ACT/ACT.BF.sigma.prior01.stan",
                            iter=40000, warmup =1000, chains=3,control = list(adapt_delta = 0.9))
ACT.BF.sigma.prior01.bridge <- bridge_sampler(ACT.BF.sigma.prior01) 
saveRDS(ACT.BF.sigma.prior01.bridge, "./BayesFactor/marginal_lik/ACT/ACT.BF.sigma.prior01.rds")

# group adjustment for sigma with prior 0,0.3 for the beta
ACT.BF.sigma.prior03 <- stan(data=ACT.data,file="./BayesFactor/models/ACT/ACT.BF.sigma.prior03.stan",
                              iter=40000, warmup =1000, chains=3,control = list(adapt_delta = 0.9))
ACT.BF.sigma.prior03.bridge <- bridge_sampler(ACT.BF.sigma.prior03) 
saveRDS(ACT.BF.sigma.prior03.bridge, "./BayesFactor/marginal_lik/ACT/ACT.BF.sigma.prior03.rds")

# group adjustment for sigma with prior 0,0.5 for the beta
ACT.BF.sigma.prior05 <- stan(data=ACT.data,file="./BayesFactor/models/ACT/ACT.BF.sigma.prior05.stan",
                               iter=40000, warmup =1000, chains=3,control = list(adapt_delta = 0.9))
ACT.BF.sigma.prior05.bridge <- bridge_sampler(ACT.BF.sigma.prior05) 
saveRDS(ACT.BF.sigma.prior05.bridge, "./BayesFactor/marginal_lik/ACT/ACT.BF.sigma.prior05.rds")


# load the fits and compute BF
ACT.BF.sigma0 <- readRDS("./BayesFactor/marginal_lik/ACT/ACT.BF.sigma0.rds")
ACT.BF.sigma.prior01 <- readRDS("./BayesFactor/marginal_lik/ACT/ACT.BF.sigma01.rds")
ACT.BF.sigma.prior03 <- readRDS("./BayesFactor/marginal_lik/ACT/ACT.BF.sigma03.rds")
ACT.BF.sigma.prior05 <- readRDS("./BayesFactor/marginal_lik/ACT/ACT.BF.sigma05.rds")

BF.sigma01.act <- (bayes_factor(ACT.BF.sigma0, ACT.BF.sigma.prior01))[[1]]
BF.sigma03.act <- (bayes_factor(ACT.BF.sigma0, ACT.BF.sigma.prior03))[[1]]
BF.sigma05.act <- (bayes_factor(ACT.BF.sigma0, ACT.BF.sigma.prior05))[[1]]

BF.sigma.act <- data.frame(priors= priors, param=rep("sigma",length(priors)),
                       BF = c(BF.sigma01.act, BF.sigma03.act, BF.sigma05.act))

BF.sigma.act.plot <- BF.sigma.act %>% 
  ggplot( aes(x=priors, y=BF)) +
  geom_line() + geom_point() + theme_minimal() +
  xlab("Prior SD") + scale_color_grey() +
  ggtitle(~sigma) + theme(plot.title = element_text(hjust = 0.5, size=22))

#### 2.Direct-access model   ---------------------------------------------------
# Create list of stan data:
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

##### 2.1. Mu ----------------------------------
# No group adjustment for mu 
DA.BF.mu0 <- stan(data=DA.standata,file="./BayesFactor/models/DA/DA.BF.mu0.stan",
                   iter=40000, warmup =1000, chains=3)
DA.BF.mu0.bridge <- bridge_sampler(DA.BF.mu0) 
saveRDS(DA.BF.mu0.bridge, "./BayesFactor/marginal_lik/DA/DA.BF.mu0.rds")
# sanity check: this model was run 3 times:
# first run -12906.39, second run -12906.34, third run

# group adjustment for mu with prior 0,0.1 for the beta 
DA.BF.mu.prior01 <- stan(data=DA.standata,file="./BayesFactor/models/DA/DA.BF.mu01.stan",
                  iter=40000, warmup =1000, chains=3)
DA.BF.mu.prior01.bridge <- bridge_sampler(DA.BF.mu.prior01) 
saveRDS(DA.BF.mu.prior01.bridge, "./BayesFactor/marginal_lik/DA/DA.BF.mu.prior01.rds")

# group adjustment for mu with prior 0,0.3 for the beta 
DA.BF.mu.prior03 <- stan(data=DA.standata,file="./BayesFactor/models/DA/DA.BF.mu03.stan",
                           iter=40000, warmup =1000, chains=3)
DA.BF.mu.prior03.bridge <- bridge_sampler(DA.BF.mu.prior03) 
saveRDS(DA.BF.mu.prior03.bridge, "./BayesFactor/marginal_lik/DA/DA.BF.mu.prior03.rds")

# group adjustment for mu with prior 0,0.5 for the beta 
DA.BF.mu.prior05 <- stan(data=DA.standata,file="./BayesFactor/models/DA/DA.BF.mu05.stan",
                           iter=40000, warmup =1000, chains=3)
DA.BF.mu.prior05.bridge <- bridge_sampler(DA.BF.mu.prior05) 
saveRDS(DA.BF.mu.prior05.bridge, "./BayesFactor/marginal_lik/DA/DA.BF.mu.prior05.rds")

#### load the fits and compute BF
DA.BF.mu0 <- readRDS("./Bayesfactor/marginal_lik/DA/DA.BF.mu0.rds")
DA.BF.mu.prior01 <- readRDS("./BayesFactor/marginal_lik/DA/DA.BF.mu.group01.rds")
DA.BF.mu.prior03 <- readRDS("./BayesFactor/marginal_lik/DA/DA.BF.mu.group03.rds")
DA.BF.mu.prior05 <- readRDS("./BayesFactor/marginal_lik/DA/DA.BF.mu.group05.rds")

DA.BF.mu.01 <- (bayes_factor(DA.BF.mu0, DA.BF.mu.prior01))[[1]]
DA.BF.mu.03 <- (bayes_factor(DA.BF.mu0, DA.BF.mu.prior03))[[1]]
DA.BF.mu.05 <- (bayes_factor(DA.BF.mu0, DA.BF.mu.prior05))[[1]]

BF.mu <- data.frame(priors= priors, param=rep("mu",length(priors)),
                       BF = c(DA.BF.mu.01, DA.BF.mu.03, DA.BF.mu.05))

##### 2.2. Theta -----------------------
# list of data for theta (one less reff for item)
DA.standata.theta0 <- within(list(),
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
                        n_w <- 2  # no. item reff 
                      }
)

# null model:
# Theta: No group adjustment for group , no interaction
DA.BF.theta0 <- stan(data=DA.standata.theta0,file="./BayesFactor/models/DA/DA.BF.theta0.stan",
                     iter=40000, warmup =1000, chains=3)
DA.BF.theta0.bridge <- bridge_sampler(DA.BF.theta0) 
saveRDS(DA.BF.theta0.bridge, "./BayesFactor/marginal_lik/DA/DA.BF.theta0.rds")


## models with group beta with prior 0,0.1
# Theta: group adjustment prior 0,0.1 , interaction prior 0,0.1
DA.BF.theta.group01.int01 <- stan(data=DA.standata,file="./BayesFactor/models/DA/DA.BF.theta.group01.int01.stan",
                              iter=40000, warmup =1000, chains=3)
DA.BF.theta.group01.int01.bridge <- bridge_sampler(DA.BF.theta.group01.int01) 
saveRDS(DA.BF.theta.group01.int01.bridge, "./BayesFactor/marginal_lik/DA.BF.theta.group01.int01.rds")

# Theta: group adjustment prior 0,0.1 , interaction prior 0,0.3
DA.BF.theta.group01.int03 <- stan(data=DA.standata,file="./BayesFactor/models/DA/DA.BF.theta.group01.int03.stan",
                                  iter=40000, warmup =1000, chains=3)
DA.BF.theta.group01.int03.bridge <- bridge_sampler(DA.BF.theta.group01.int03) 
saveRDS(DA.BF.theta.group01.int03.bridge, "./BayesFactor/marginal_lik/DA/DA.BF.theta.group01.int03.rds")

# Theta: group adjustment prior 0,0.1 , interaction prior 0,0.5
DA.BF.theta.group01.int05 <- stan(data=DA.standata,file="./BayesFactor/models/DA/DA.BF.theta.group01.int05.stan",
                                  iter=40000, warmup =1000, chains=3)
DA.BF.theta.group01.int05.bridge <- bridge_sampler(DA.BF.theta.group01.int05) 
saveRDS(DA.BF.theta.group01.int05.bridge, "./BayesFactor/marginal_lik/DA/DA.BF.theta.group01.int05.rds")


## models with group beta with prior 0,0.3
# Theta: group adjustment prior 0,0.3 , interaction prior 0,0.1
DA.BF.theta.group03.int01 <- stan(data=DA.standata,file="./BayesFactor/models/DA/DA.BF.theta.group03.int01.stan",
                                  iter=40000, warmup =1000, chains=3)
DA.BF.theta.group03.int01.bridge <- bridge_sampler(DA.BF.theta.group03.int01) 
saveRDS(DA.BF.theta.group03.int01.bridge, "./BayesFactor/marginal_lik/DA.BF.theta.group03.int01.rds")

# Theta: group adjustment prior 0,0.3 , interaction prior 0,0.3
DA.BF.theta.group03.int03 <- stan(data=DA.standata,file="./BayesFactor/models/DA/DA.BF.theta.group03.int03.stan",
                                  iter=40000, warmup =1000, chains=3)
DA.BF.theta.group03.int03.bridge <- bridge_sampler(DA.BF.theta.group03.int03) 
saveRDS(DA.BF.theta.group03.int03.bridge, "./BayesFactor/marginal_lik/DA/DA.BF.theta.group03.int03.rds")

# Theta: group adjustment prior 0,0.3 , interaction prior 0,0.5
DA.BF.theta.group03.int05 <- stan(data=DA.standata,file="./BayesFactor/models/DA/DA.BF.theta.group03.int05.stan",
                                  iter=40000, warmup =1000, chains=3)
DA.BF.theta.group03.int05.bridge <- bridge_sampler(DA.BF.theta.group03.int05) 
saveRDS(DA.BF.theta.group03.int05.bridge, "./BayesFactor/marginal_lik/DA/DA.BF.theta.group03.int05.rds")

## models with group beta with prior 0,0.5
# Theta: group adjustment prior 0,0.5 , interaction prior 0,0.1
DA.BF.theta.group05.int01 <- stan(data=DA.standata,file="./BayesFactor/models/DA/DA.BF.theta.group05.int01.stan",
                                  iter=40000, warmup =1000, chains=3)
DA.BF.theta.group05.int01.bridge <- bridge_sampler(DA.BF.theta.group05.int01) 
saveRDS(DA.BF.theta.group05.int01.bridge, "./BayesFactor/marginal_lik/DA/DA.BF.theta.group05.int01.rds")

# Theta: group adjustment prior 0,0.5 , interaction prior 0,0.3
DA.BF.theta.group05.int03 <- stan(data=DA.standata,file="./BayesFactor/models/DA/DA.BF.theta.group05.int03.stan",
                                  iter=40000, warmup =1000, chains=3)
DA.BF.theta.group05.int03.bridge <- bridge_sampler(DA.BF.theta.group05.int03) 
saveRDS(DA.BF.theta.group05.int03.bridge, "./BayesFactor/marginal_lik/DA/DA.BF.theta.group05.int03.rds")

# Theta: group adjustment prior 0,0.5 , interaction prior 0,0.5
DA.BF.theta.group05.int05 <- stan(data=DA.standata,file="./BayesFactor/models/DA/DA.BF.theta.group05.int05.stan",
                                  iter=40000, warmup =1000, chains=3)
DA.BF.theta.group05.int05.bridge <- bridge_sampler(DA.BF.theta.group05.int05) 
saveRDS(DA.BF.theta.group05.int05.bridge, "./BayesFactor/marginal_lik/DA/DA.BF.theta.group05.int05.rds")

## 2.2.1 BF and plots ----
# load the fits
DA.BF.theta0 <- readRDS("./BayesFactor/marginal_lik/DA/DA.BF.theta0.rds")

DA.BF.theta.group01.int01 <- readRDS("./BayesFactor/marginal_lik/DA/DA.BF.theta.group01.int01.rds")
DA.BF.theta.group01.int03 <- readRDS("./BayesFactor/marginal_lik/DA/DA.BF.theta.group01.int01.rds")
DA.BF.theta.group01.int05 <- readRDS("./BayesFactor/marginal_lik/DA/DA.BF.theta.group01.int05.rds")

DA.BF.theta.group03.int01 <- readRDS("./BayesFactor/marginal_lik/DA/DA.BF.theta.group03.int01.rds")
DA.BF.theta.group03.int03 <- readRDS("./BayesFactor/marginal_lik/DA/DA.BF.theta.group03.int03.rds")
DA.BF.theta.group03.int05 <- readRDS("./BayesFactor/marginal_lik/DA/DA.BF.theta.group03.int05.rds")

DA.BF.theta.group05.int01 <- readRDS("./BayesFactor/marginal_lik/DA/DA.BF.theta.group05.int01.rds")
DA.BF.theta.group05.int03 <- readRDS("./BayesFactor/marginal_lik/DA/DA.BF.theta.group05.int03.rds")
DA.BF.theta.group05.int05 <- readRDS("./BayesFactor/marginal_lik/DA/DA.BF.theta.group05.int05.rds")

# prior sd 0.01 for group beta and varying prior for interaction
DA.BF.theta.g01.i01 <- (bayes_factor(DA.BF.theta0, DA.BF.theta.group01.int01))[[1]]
DA.BF.theta.g01.i03 <- (bayes_factor(DA.BF.theta0, DA.BF.theta.group01.int03))[[1]]
DA.BF.theta.g01.i05 <- (bayes_factor(DA.BF.theta0, DA.BF.theta.group01.int05))[[1]]

# prior sd 0.03 for group beta and varying prior for interaction
DA.BF.theta.g03.i01 <- (bayes_factor(DA.BF.theta0, DA.BF.theta.group03.int01))[[1]]
DA.BF.theta.g03.i03 <- (bayes_factor(DA.BF.theta0, DA.BF.theta.group03.int03))[[1]]
DA.BF.theta.g03.i05 <- (bayes_factor(DA.BF.theta0, DA.BF.theta.group03.int05))[[1]]

# prior sd 0.05 for group beta and varying prior for interaction
DA.BF.theta.g05.i01 <- (bayes_factor(DA.BF.theta0, DA.BF.theta.group05.int01))[[1]]
DA.BF.theta.g05.i03 <- (bayes_factor(DA.BF.theta0, DA.BF.theta.group05.int03))[[1]]
DA.BF.theta.g05.i05 <- (bayes_factor(DA.BF.theta0, DA.BF.theta.group05.int05))[[1]]


BF.theta <- data.frame(group = c(0.1,0.1,0.1, 
                                 0.3, 0.3, 0.3, 
                                 0.5, 0.5, 0.5), interaction = rep(c(0.1,0.3,0.5),3), 
                        BF = c(DA.BF.theta.g01.i01, DA.BF.theta.g01.i03, DA.BF.theta.g01.i05,
                             DA.BF.theta.g03.i01, DA.BF.theta.g03.i03, DA.BF.theta.g03.i05,
                              DA.BF.theta.g05.i01, DA.BF.theta.g05.i03, DA.BF.theta.g05.i05)) 

# scatterplot3d(BF.theta)
# scatterplot3d(x=BF.theta$group, y=BF.theta$BF, z=BF.theta$interaction,
#               grid=TRUE, box=FALSE, color="black", pch = 16, 
#               xlab = "Prior SD for the group" ~beta, ylab="BF",
#               zlab="Prior SD for the interaction" ~beta,
#               title(main="BF with different prior SD for the group and the 
#           interaction adjustments in the parameter" ~theta))
#                      

##### 2.3. Pb --------------------------
# No group adjustment for Pb
DA.BF.Pb0 <- stan(data=DA.standata,file="./BayesFactor/models/DA/DA.BF.Pb0.stan",
                     iter=40000, warmup =1000, chains=3)
DA.BF.Pb0.bridge <- bridge_sampler(DA.BF.Pb0) 
saveRDS(DA.BF.Pb0.bridge, "./BayesFactor/marginal_lik/DA/DA.BF.Pb0.rds")

# group adjustment for Pb with prior 0,0.1 for the beta 
DA.BF.Pb.prior01 <- stan(data=DA.standata,file="./BayesFactor/models/DA/DA.BF.Pb.group01.stan",
                           iter=40000, warmup =1000, chains=3)
DA.BF.Pb.prior01.bridge <- bridge_sampler(DA.BF.Pb.prior01) 
saveRDS(DA.BF.Pb.prior01.bridge, "./BayesFactor/marginal_lik/DA/DA.BF.Pb.group01.rds") 

# group adjustment for Pb with prior 0,0.3 for the beta 
DA.BF.Pb.prior03 <- stan(data=DA.standata,file="./BayesFactor/models/DA/DA.BF.Pb.group03.stan",
                  iter=40000, warmup =1000, chains=3)
DA.BF.Pb.prior03.bridge <- bridge_sampler(DA.BF.Pb.prior03) 
saveRDS(DA.BF.Pb.prior03.bridge, "./BayesFactor/marginal_lik/DA/DA.BF.Pb.group03.rds") 

# group adjustment for Pb with prior 0,0.5 for the beta 
DA.BF.Pb.prior05 <- stan(data=DA.standata,file="./BayesFactor/models/DA/DA.BF.Pb.group05.stan",
                           iter=40000, warmup =1000, chains=3)
DA.BF.Pb.prior05.bridge <- bridge_sampler(DA.BF.Pb.prior05) 
saveRDS(DA.BF.Pb.prior05.bridge, "./BayesFactor/marginal_lik/DA/DA.BF.Pb.group05.rds") 


# load the fit and perform BF
DA.BF.Pb0 <- readRDS("./Bayesfactor/marginal_lik/DA/DA.BF.Pb0.rds")
DA.BF.Pb01 <- readRDS("./BayesFactor/marginal_lik/DA/DA.BF.Pb.group01.rds")
DA.BF.Pb03 <- readRDS("./BayesFactor/marginal_lik/DA/DA.BF.Pb.group03.rds")
DA.BF.Pb05 <- readRDS("./BayesFactor/marginal_lik/DA/DA.BF.Pb.group05.rds")

DA.BF.Pb.01 <- (bayes_factor(DA.BF.Pb0, DA.BF.Pb01))[[1]]
DA.BF.Pb.03 <- (bayes_factor(DA.BF.Pb0, DA.BF.Pb03))[[1]]
DA.BF.Pb.05 <- (bayes_factor(DA.BF.Pb0, DA.BF.Pb05))[[1]]

BF.Pb <- data.frame(priors= priors, param=rep("Pb",length(priors)),
                    BF = c(DA.BF.Pb.01, DA.BF.Pb.03, DA.BF.Pb.05))

##### 2.4. Delta -----------------------
# No group adjustment for delta
DA.BF.delta0 <- stan(data=DA.standata,file="./BayesFactor/models/DA/DA.BF.delta0.stan",
                  iter=40000, warmup =1000, chains=3)
DA.BF.delta0.bridge <- bridge_sampler(DA.BF.delta0) 
saveRDS(DA.BF.delta0.bridge, "./BayesFactor/marginal_lik/DA/DA.BF.delta0.rds")

# group adjustment for delta with prior 0,0.1 for the beta 
DA.BF.delta01 <- stan(data=DA.standata,file="./BayesFactor/models/DA/DA.BF.delta01.stan",
                           iter=40000, warmup =1000, chains=3)
DA.BF.delta.prior01.bridge <- bridge_sampler(DA.BF.delta.prior01) 
saveRDS(DA.BF.delta.prior01.bridge, "./BayesFactor/marginal_lik/DA/DA.BF.delta01.rds") 

# group adjustment for delta with prior 0,0.3 for the beta 
DA.BF.delta03 <- stan(data=DA.standata,file="./BayesFactor/models/DA/DA.BF.delta03.stan",
                           iter=40000, warmup =1000, chains=3)
DA.BF.delta03.bridge <- bridge_sampler(DA.BF.delta03) 
saveRDS(DA.BF.delta03.bridge, "./BayesFactor/marginal_lik/DA/DA.BF.delta03.rds") 

# group adjustment for delta with prior 0,0.5 for the beta 
DA.BF.delta05 <- stan(data=DA.standata,file="./BayesFactor/models/DA/DA.BF.delta05.stan",
                           iter=40000, warmup =1000, chains=3)
DA.BF.delta05.bridge <- bridge_sampler(DA.BF.delta05) 
saveRDS(DA.BF.delta05.bridge, "./BayesFactor/marginal_lik/DA/DA.BF.delta05.rds") 

# load the fit and perform BF
DA.BF.delta0 <- readRDS("./Bayesfactor/marginal_lik/DA/DA.BF.delta0.rds")
DA.BF.delta.prior01 <- readRDS("./BayesFactor/marginal_lik/DA/DA.BF.delta.group01.rds")
DA.BF.delta.prior03 <- readRDS("./BayesFactor/marginal_lik/DA/DA.BF.delta.group03.rds")
DA.BF.delta.prior05 <- readRDS("./BayesFactor/marginal_lik/DA/DA.BF.delta.group05.rds")

DA.BF.delta.01 <- (bayes_factor(DA.BF.delta0, DA.BF.delta.prior01))[[1]]
DA.BF.delta.03 <- (bayes_factor(DA.BF.delta0, DA.BF.delta.prior03))[[1]]
DA.BF.delta.05 <- (bayes_factor(DA.BF.delta0, DA.BF.delta.prior05))[[1]]

BF.delta <- data.frame(priors= priors, param=rep("delta",length(priors)),
                       BF = c(DA.BF.delta.01, DA.BF.delta.03, DA.BF.delta.05))

##### 2.5. Sigma -----------------------
# No group adjustment for sigma
DA.BF.sigma0 <- stan(data=DA.standata,file="./BayesFactor/models/DA/DA.BF.sigma0.stan",
                     iter=40000, warmup =1000, chains=3)
DA.BF.sigma0.bridge <- bridge_sampler(DA.BF.sigma0) 
saveRDS(DA.BF.sigma0.bridge, "./BayesFactor/marginal_lik/DA/DA.BF.sigma0.rds")

# group adjustment for sigma with prior 0,0.1 for the beta 
DA.BF.sigma.prior01 <- stan(data=DA.standata,file="./BayesFactor/models/DA/DA.BF.sigma01.stan",
                     iter=40000, warmup =1000, chains=3)
DA.BF.sigma.prior01.bridge <- bridge_sampler(DA.BF.sigma.prior01) 
saveRDS(DA.BF.sigma.prior01.bridge, "./BayesFactor/marginal_lik/DA/DA.BF.sigma.prior01.rds")

# group adjustment for sigma with prior 0,0.3 for the beta 
DA.BF.sigma.prior03 <- stan(data=DA.standata,file="./BayesFactor/models/DA/DA.BF.sigma03.stan",
                              iter=40000, warmup =1000, chains=3)
DA.BF.sigma.prior03.bridge <- bridge_sampler(DA.BF.sigma.prior03) 
saveRDS(DA.BF.sigma.prior03.bridge, "./BayesFactor/marginal_lik/DA/DA.BF.sigma.prior03.rds")

# group adjustment for sigma with prior 0,0.5 for the beta 
DA.BF.sigma.prior05 <- stan(data=DA.standata,file="./BayesFactor/models/DA/DA.BF.sigma05.stan",
                              iter=40000, warmup =1000, chains=3)
DA.BF.sigma.prior05.bridge <- bridge_sampler(DA.BF.sigma.prior05) 
saveRDS(DA.BF.sigma.prior05.bridge, "./BayesFactor/marginal_lik/DA/DA.BF.sigma.prior05.rds")

# load the fits and compute BF
DA.BF.sigma0 <- readRDS("./BayesFactor/marginal_lik/DA/DA.BF.sigma0.rds")
DA.BF.sigma.prior01 <- readRDS("./BayesFactor/marginal_lik/DA/DA.BF.sigma.group01.rds")
DA.BF.sigma.prior03 <- readRDS("./BayesFactor/marginal_lik/DA/DA.BF.sigma.group03.rds")
DA.BF.sigma.prior05 <- readRDS("./BayesFactor/marginal_lik/DA/DA.BF.sigma.group05.rds")

DA.BF.sigma.01 <- (bayes_factor(DA.BF.sigma0, DA.BF.sigma.prior01))[[1]]
DA.BF.sigma.03 <- (bayes_factor(DA.BF.sigma0, DA.BF.sigma.prior03))[[1]]
DA.BF.sigma.05 <- (bayes_factor(DA.BF.sigma0, DA.BF.sigma.prior05))[[1]]

BF.sigma <- data.frame(priors= priors, param=rep("sigma",length(priors)),
                       BF = c(DA.BF.sigma.01, DA.BF.sigma.03, DA.BF.sigma.05))

# plot all BF together
BF.delta.Pb <- rbind(BF.delta, BF.Pb)
BF.delta.Pb.plot <- BF.delta.Pb%>% 
  ggplot( aes(x=priors, y=BF,color=param, shape=param)) +
  geom_line() + geom_point()  +
  ylim(-1,12)+ theme_minimal() + 
  geom_hline(aes(yintercept=1), linetype="dotted") +
  xlab("Prior SD") + scale_color_grey() + 
  labs(color="Parameter", shape="Parameter")+
  ggtitle(" ")

BF.DA.mu.sigma <- rbind(BF.mu, BF.sigma)
BF.DA.mu.sigma.plot <- BF.DA.mu.sigma %>% 
  ggplot( aes(x=priors, y=BF,color=param, shape=param)) +
  geom_line() + geom_point() + theme_minimal() +
  xlab("Prior SD") + scale_color_grey() + 
  labs(color="Parameter", shape="Parameter") +
  ggtitle(" ")

grid.arrange(BF.DA.mu.sigma.plot, BF.delta.Pb.plot)

