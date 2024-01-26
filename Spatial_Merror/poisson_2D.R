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

setwd("D:/Study new/PhD/R code/Spatial_Merror")
#gdbsource("poisson_2D.R",interactive=TRUE) #make sure this is commented out, or you will get an infinite loop when debugging
TMB::compile("poisson_sglmm.cpp","-O1 -g",DLLFLAGS="")


dyn.load(dynlib("poisson_sglmm"))


#Dataset generation


set.seed(10)
n <- 1000 # Bumped up sample size, as one of the things we want to show is that scalability of FRK/GAMs as an approach
true_coefs <- c(1,1)
xy <- data.frame(x = runif(n), y = runif(n))

## Generate covariates
num_X = 1
sigma_e = sqrt(0.1)
X <- RFsimulate(model = RMmatern(nu = 1, var = 1, scale = 0.3), x = xy$x, y = xy$y, n = 1)@data[,1] # Using a stationary spatial field
W <- X + rnorm(n, sd = sigma_e )


## Generate response
spranef <- RFsimulate(model = RMmatern(nu = 0.5, var = 1, scale = 0.3), x = xy$x, y = xy$y, n = 1)@data[,1] # Gaussian field with exponential covariance function 


eta <- cbind(1,X) %*% true_coefs + spranef
simy <- rpois(n, lambda = exp(eta))



#Generate autoFRK basis functions
sp_basisfunctions1 <- mrts(xy, k = 5, x = xy) %>% 
  as.matrix
sp_basisfunctions1 <- sp_basisfunctions1[,-1] #no intercept

sp_basisfunctions2 <- mrts(xy, k = 20, x = xy) %>% 
  as.matrix

#Get some 'reasonable' starting values

mod = lm(W ~ -1 + sp_basisfunctions2)
mod1 = glm(simy ~ W + sp_basisfunctions1 , family = poisson)
modsum = summary(mod)

cholSigma1_vec = diag(log(abs(mod1$coefficients[-(1:2)])))[lower.tri(diag(nrow =  ncol(sp_basisfunctions1)), diag = TRUE)]
cholSigma2_vec = diag(log(abs(mod$coefficients)))[lower.tri(diag(nrow =  ncol(sp_basisfunctions2)), diag = TRUE)]

# cholSigma1_vec = diag(nrow =  ncol(sp_basisfunctions1))[lower.tri(diag(nrow =  ncol(sp_basisfunctions1)), diag = TRUE)]
# cholSigma2_vec = diag(nrow =  ncol(sp_basisfunctions2))[lower.tri(diag(nrow =  ncol(sp_basisfunctions2)), diag = TRUE)]


inputdata=list(sp_basis1 = sp_basisfunctions1, sp_basis2 = sp_basisfunctions2 , y = as.matrix(simy), w = as.matrix(W)) 
parameterlist=list(beta = true_coefs , cholSigma1 = cholSigma1_vec , cholSigma2 = cholSigma2_vec , 
                   theta1 = mod1$coefficients[-(1:2)] , X_thetas = mod$coefficients , logsigma_e = log(modsum$sigma) )

#mod1$coefficients[1:2] instead of true_coefs could be used, but

obj <- MakeADFun(data=inputdata,DLL = "poisson_sglmm",
                 parameters = parameterlist, random = c("theta1", "X_thetas"), silent=T)


#set a lower limit for the optimization because the diagonals of the Cholesky decomposition cannot be negative. 


fit <- optim(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS", hessian = TRUE, control=list(maxit=1000)) #obj$fn is already the (negative) marginal log-likelihood - not the joint.

#these will be the (psuedo) MLEs 

fit$par
fit$convergence

#make sure you do the following steps AFTER the optim() step.

rep1 <- sdreport(obj)
rep1 #MLEs
summary(rep1) #MLEs and also "estimated random effects"
as.list(rep1,"Est",FALSE) #easier way to look at the results
as.list(rep1,"Est",TRUE) #reports Sigma here, because of the ADREPORT(Sigma) line in the .cpp file.

#for some more options when reporting results, see the help files.









