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

#set.seed(23)
cl <- makeCluster(detectCores()-1)
# #cl <- detectCores()-2
registerDoParallel(cl)



sim_number = 1000
num_obs = 500
tau = 1
sigma_sq_g = 1
sigma_sq_e = 0.4
sigma_sq_u = 0.5
beta0 = 10
beta1 = 2

adj.fits = rep(0,sim_number)
ols.fits = rep(0,sim_number)
gls.fits = rep(0,sim_number)
knownx.fits = rep(0,sim_number)
adjols.fits = rep(0,sim_number)

adj.fits.se = rep(0,sim_number)
adjols.fits.se = rep(0,sim_number)



results = foreach(k = 1:sim_number,.packages = c("mvtnorm","mgcv","Matrix")) %dopar% {
  
  
  #generate spatial locations on a line
  
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
  
  
  #unconditional on W; fit model to W
  W = X + rnorm(n =  num_obs, mean = 0, sd = sqrt(sigma_sq_u))
  W.fit = gam(W ~  s(S))
  pred.X = W.fit$fitted.values

  alpha = c(rmvnorm(1, sigma = G)) #put this inside loop for unconditional, outside for conditional
  epsilon = rnorm(n = num_obs,mean = 0,sd = sqrt(sigma_sq_e))
  Y = beta0 + X*beta1 + alpha + epsilon
  
  #proposed method

  Y.fit = gam(Y ~  pred.X + s(S,bs="cr",k=3))
  adj.beta1 = Y.fit$coefficients[2]
  Y.fit_summ = summary(Y.fit)
  Y.fit.se = Y.fit_summ$se[2]
  
  #naive ols
  
  ols.fit = lm(Y ~ W)
  ols.beta1 = ols.fit$coefficients[2]
  
  
  #naive gls
  
  omega = G + sigma_sq_e*diag(ncol(G))
  Wstar = cbind(1,W)
  gls.beta = solve( crossprod(Wstar,solve(omega))%*%Wstar ) %*% crossprod(Wstar, solve(omega)) %*% Y
  gls.beta1 = gls.beta[2]
  
  
  #known x (experimenting)
  
  knownx.fit = gam(Y ~  X + s(S,bs="cr",k=3))
  knownx.beta1 = knownx.fit$coefficients[2]
  
  
  #ols with pred.X
  
  adjols.fit = lm(Y ~ pred.X)
  adjols.beta1 = adjols.fit$coefficients[2]
  adjols.fit_summ = summary(adjols.fit)
  adjols.fit.se = adjols.fit_summ$coefficients[2,2]
  
  out = list(adj.fits = adj.beta1, ols.fits = ols.beta1, gls.fits = gls.beta1 , knownx.fits = knownx.beta1,
             adjols.fits = adjols.beta1, adj.fits.se = Y.fit.se, adjols.fits.se = adjols.fit.se )
}


for (k in 1:sim_number) {
  adj.fits[k] = results[[k]]$adj.fits
  ols.fits[k] = results[[k]]$ols.fits
  gls.fits[k] = results[[k]]$gls.fits
  knownx.fits[k] = results[[k]]$knownx.fits
  adjols.fits[k] = results[[k]]$adjols.fits
  adj.fits.se[k] = results[[k]]$adj.fits.se
  adjols.fits.se[k] = results[[k]]$adjols.fits.se
}

mean(adj.fits)
mean(ols.fits)
mean(gls.fits)
mean(knownx.fits)
median(knownx.fits)
mean(adjols.fits)


sqrt(var(adj.fits))
mean(adj.fits.se)
sqrt(var(adjols.fits))
mean(adjols.fits.se)

stopCluster(cl)
