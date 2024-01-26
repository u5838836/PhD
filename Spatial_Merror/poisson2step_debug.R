rm(list=ls())
library(TMB)
library(mvtnorm)
library(autoFRK)
library(mgcv)
library(spaMM)
#library(tidyverse)
library(RandomFields)
library(sp)
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(fields)
library(optimx)

setwd("D:/Study new/PhD/R code/Spatial_Merror")
#gdbsource("poisson2D_2step.R",interactive=TRUE)
TMB::compile("poisson_2step.cpp","-O1 -g",DLLFLAGS="")


dyn.load(dynlib("poisson_2step"))

  #Dataset generation
  
  set.seed(61)
  
  n <- 1000 # Bumped up sample size, as one of the things we want to show is that scalability of FRK/GAMs as an approach
  true_coefs <- c(1,1)
  xy <- data.frame(x = runif(n), y = runif(n))
  
  ## Generate covariates
  num_X = 1
  sigma_e = sqrt(0.5)
  X <- RFsimulate(model = RMmatern(nu = 1, var = 1, scale = 0.3), x = xy$x, y = xy$y, n = 1)@data[,1] # Using a stationary spatial field
  W <- X + rnorm(n, sd = sigma_e )
  
  
  ## Generate response
  spranef <- RFsimulate(model = RMmatern(nu = 0.5, var = 1, scale = 0.3), x = xy$x, y = xy$y, n = 1)@data[,1] # Gaussian field with exponential covariance function 
  
  
  eta <- cbind(1,X) %*% true_coefs + spranef
  simy <- rpois(n, lambda = exp(eta))
  
  
  
  #Generate autoFRK basis functions
  sp_basisfunctions1 <- mrts(xy, k = 10, x = xy) %>% 
    as.matrix
  sp_basisfunctions1 <- sp_basisfunctions1[,-1] #no intercept
  
  sp_basisfunctions2 <- mrts(xy, k = 20, x = xy) %>% 
    as.matrix
  
  #Fit the model for W via autoFRK()
  fit_Wmodel <- autoFRK(Data = matrix(W,ncol=1), loc = xy, G = sp_basisfunctions2)
  X_pred <- predict(fit_Wmodel)$pred
  
  
  #Get some 'reasonable' starting values
  
  mod1 = glm(simy ~ W + sp_basisfunctions1 , family = poisson)
  mod2 = glm(simy ~ X_pred , family = poisson)
  
  cholSigma1_vec = diag(log(abs(mod1$coefficients[-(1:2)])))[lower.tri(diag(nrow =  ncol(sp_basisfunctions1)), diag = TRUE)]
  
  
  # cholSigma1_vec = diag(nrow =  ncol(sp_basisfunctions1))[lower.tri(diag(nrow =  ncol(sp_basisfunctions1)), diag = TRUE)]
  # cholSigma2_vec = diag(nrow =  ncol(sp_basisfunctions2))[lower.tri(diag(nrow =  ncol(sp_basisfunctions2)), diag = TRUE)]
  
  
  inputdata=list(sp_basis1 = sp_basisfunctions1 , y = as.matrix(simy), predX = X_pred) 
  parameterlist=list(beta = mod1$coefficients[1:2] , cholSigma1 = cholSigma1_vec , theta1 = mod1$coefficients[-(1:2)] )
  parameterlist=list(beta = fit$par[1:2] , cholSigma1 = fit$par[-(1:2)] , theta1 = rep(0,ncol(sp_basisfunctions1)) )
  

  
  obj <- MakeADFun(data=inputdata,DLL = "poisson_2step",
                   parameters = parameterlist, random = c("theta1"), silent=T)
  
  
  
  skip_to_next <- FALSE
  
  
  # tryCatch(fit <- optimr(par = obj$par, fn = obj$fn, gr = obj$gr, method = "Rcgmin", hessian = TRUE,
  #                       control=list(maxit=500, trace=1, reltol = '1e-8')) , error = function(e) { skip_to_next <<- TRUE} )

  tryCatch(fit <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr,
                         control=list(trace = 1, iter.max = 1000, eval.max = 1000, rel.tol = '1e-10')), error = function(e) { skip_to_next <<- TRUE})

  if(skip_to_next) { next } 
  
  fit$convergence

  fit$par
  