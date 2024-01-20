rm(list = ls())
library(glmmTMB)
library(dplyr)
library(foreach)
library(glm2)
library(lme4)

num_clus <- 10
num_timepoints <- 10

true_beta = 1
true_G = 0.1

chol_start = log(sqrt(true_G)) #fix working G for glmmTMB

# Generate the random effects
true_alpha <- rnorm(num_clus, sd = sqrt(true_G)) %>% matrix(nrow=num_clus)

# Generate the observations
simdat <- data.frame(ID = factor(rep(1:num_clus, each = num_timepoints)), 
                     time = factor(rep(1:num_timepoints, num_clus)))

simdat$eta <- c(true_beta + true_alpha[simdat$ID,])
simdat$y <- rgamma(n = nrow(simdat), shape = 1, scale = exp(simdat$eta))

##--------------
## GLMM fits via Laplace. 
##--------------

fit_resp <- glmmTMB(formula = y ~ 1 + (1|ID), family = Gamma(link='log'),data=simdat, 
                    map = list(theta = factor(rep(NA,length(chol_start)))), start = list(theta = chol_start))
fit_lme <- glmer(formula = y ~ 1 + (1|ID), family = Gamma(link='log'),data=simdat)
ranef(fit_resp) %>% as.data.frame
ranef(fit_lme) %>% as.data.frame









