rm(list = ls())
library(nlme)
library(mgcv)
library(glmmTMB)
#library(e1071)
library(dplyr)
library(xtable)
library(mvtnorm)
library(MASS)
library(Matrix)
library(parallel)
library(foreach)
library(doParallel)
#library(tidyverse)
library(GGally)
library(lme4)
library(glm2)
#library(doMC)
library(xtable)
#registerDoMC(cores=28) # this should equal ncpus

#setwd("D:/Study new/PhD/R code")

#set.seed(2)
cl <- makeCluster(detectCores()-1)
# #cl <- detectCores()-2
registerDoParallel(cl)



sim_number = 1000
num_obs = 500
tau = 1 #controls strength of spatial correlation; high -> high correlation
sigma_sq_g = 0.5 #also controls strength of spatial correlation
sigma_sq_u = 1 #covariate measurement error
beta0 = 1
beta1 = 0.2

adj.fits = rep(0,sim_number)
naive.fits = rep(0,sim_number)
knownx.fits = rep(0,sim_number)
no_spatial.fits = rep(0,sim_number)

adj.fits.se = rep(0,sim_number)
no_spatial.fits.se = rep(0,sim_number)




results = foreach(k = 1:sim_number,.packages = c("mvtnorm","mgcv","Matrix")) %dopar% {
  
  #generate locations on a line, and construct correlation/covariance matrix
  
  S = runif(num_obs)
  R = matrix(nrow=num_obs,ncol=num_obs)
  for (i in 1:nrow(R)) {
    for (j in 1:ncol(R)) {
      R[i,j] = exp(- abs(S[i]-S[j])/tau )
    }
  }
  G = sigma_sq_g*R
  
  #generate X based on S
  
  X = 1/(1+S) + 3*exp(-50*(S-0.3)^2) + 2*exp(-25*(S-0.7)^2) 
  
  
  # fit W semi-parametrically; unconditional on W
  W = X + rnorm(n =  num_obs, mean = 0, sd = sqrt(sigma_sq_u))
  W.fit = gam(W ~  s(S,bs="cr",k=5))
  pred.X = W.fit$fitted.values
  
  
  #generate from poisson log link glmm
  
  alpha = c(rmvnorm(1, sigma = G)) #put this inside loop for unconditional, outside for conditional
  eta = beta0 + X*beta1 + alpha
  Y = rpois(n = num_obs, lambda = exp(eta))
  
  
  
  #get estimate of beta1 using proposed method
  
  Y.fit = gam(Y ~  pred.X + s(S,bs="cr",k=5), family = poisson)
  adj.beta1 = Y.fit$coefficients[2]
  Y.fit_summ = summary(Y.fit)
  Y.fit.se = Y.fit_summ$se[2]
  
  
  #estimate of beta1 if X was known
  
  knownx.fit = gam(Y ~  X + s(S,bs="cr",k=5), family = poisson)
  knownx.beta1 = knownx.fit$coefficients[2]
  
  
  #naive estimate of beta1, no adjustment for covariate error (but does adjust for spatial structure)
  
  naive.fit = gam(Y ~ W + s(S,bs="cr",k=5), family = poisson)
  naive.beta1 = naive.fit$coefficients[2]
  
  
  #adjust for covariate error, but not spatial structure
  
  # no_spatial.fit = glm(Y ~ pred.X, family = poisson)
  # no_spatial.beta1 = no_spatial.fit$coefficients[2]
  # no_spatial.fit_summ = summary(no_spatial.fit)
  # no_spatial.fit.se = no_spatial.fit_summ$coefficients[2,2]
  
  
  no_spatial.fit = gam(Y ~ pred.X, family = poisson)
  no_spatial.beta1 = no_spatial.fit$coefficients[2]
  no_spatial.fit_summ = summary(no_spatial.fit)
  no_spatial.fit.se = no_spatial.fit_summ$se[2]
  
  
  out = list(adj.fits = adj.beta1, naive.fits = naive.beta1, knownx.fits = knownx.beta1,
             no_spatial.fits = no_spatial.beta1, adj.fits.se = Y.fit.se, no_spatial.fits.se = no_spatial.fit.se)
}


for (k in 1:sim_number) {
  adj.fits[k] = results[[k]]$adj.fits
  naive.fits[k] = results[[k]]$naive.fits
  knownx.fits[k] = results[[k]]$knownx.fits
  no_spatial.fits[k] = results[[k]]$no_spatial.fits
  adj.fits.se[k] = results[[k]]$adj.fits.se
  no_spatial.fits.se[k] = results[[k]]$no_spatial.fits.se
}

mean(adj.fits)
mean(naive.fits)
mean(knownx.fits)
mean(no_spatial.fits)


sqrt(var(adj.fits))
mean(adj.fits.se)
sqrt(var(no_spatial.fits))
mean(no_spatial.fits.se)

stopCluster(cl)

