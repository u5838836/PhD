##--------------------
## Example of autoFRK use for spatial measurement error setting
##--------------------
rm(list = ls())
library(autoFRK)
library(mgcv)
library(spaMM)
#library(mvtnorm)
library(tidyverse)
library(RandomFields)
library(sp)

set.seed(102021)
n <- 1000 # Bumped up sample size, as one of the things we want to show is that scalability of FRK/GAMs as an approach
true_coefs <- c(1,1)
xy <- data.frame(x = runif(n), y = runif(n))

## Generate covariates
#X <- 1/(1+xy) + 3*exp(-50*(xy[,1]-0.3)^2) + 2*exp(-25*(xy[,2]-0.7)^2) # From Huque et al. (2016)
#X <- apply(X, 1, prod)
#W <- X + rnorm(n, sd = 1)
X <- RFsimulate(model = RMmatern(nu = 1, var = 1, scale = 0.3), x = xy$x, y = xy$y, n = 1)@data[,1] # Using a stationary spatial field
W <- X + rnorm(n, sd = sqrt(0.1))
ggplot(data.frame(xy, z = X), aes(x = x, y = y, color = z)) +
   geom_point() +
   scale_color_viridis_c() +
   theme_bw()
var(X)
qplot(X, W) +
   geom_abline(intercept = 0, slope = 1) +
   theme_bw()


## Generate response
spranef <- RFsimulate(model = RMmatern(nu = 0.5, var = 1, scale = 0.3), x = xy$x, y = xy$y, n = 1)@data[,1] # Gaussian field with exponential covariance function 
ggplot(data.frame(xy, z = spranef), aes(x = x, y = y, color = z)) +
   geom_point() +
   scale_color_viridis_c() +
   theme_bw()

eta <- cbind(1,X) %*% true_coefs + spranef
#simy <- rnorm(n, mean = eta, sd = 0.1)
#simy <- rbinom(n, size = 5, p = plogis(eta))
simy <- rpois(n, lambda = exp(eta))

dat <- data.frame(xy, X = X, W = W, response = simy, size = 10)
head(dat)

ggplot(dat, aes(x = x, y = y, color = response)) +
   geom_point() +
   scale_color_viridis_c() +
   theme_bw()


##-------------------------------
## Models which assume you know the true covariate
#fit_X_ns <- glm(cbind(response, size-response) ~ X, data = dat, family =  "binomial")
fit_X_ns <- glm(response ~ X, data = dat, family =  "poisson")
summary(fit_X_ns)

#fit_X_sp <- gam(cbind(response, size-response) ~ X + s(x, y), data = dat, family = "binomial", method = "ML") #, bs = "gp"
fit_X_sp <- gam(response ~ X + s(x, y), data = dat, family = "poisson", method = "ML") #, bs = "gp"
summary(fit_X_sp)


## Models which assume you only know the measured covariate
#fit_ns <- glm(cbind(response, size-response) ~ W, data = dat, family = "binomial")
fit_ns <- glm(response ~ W, data = dat, family = "poisson")
summary(fit_ns)

#fit_sp <- fitme(cbind(response, size-response) ~ W + Matern(1|x+y), data = dat, family = "binomial") # This can take a while, but it partly demonstrates why spatial (G)LMMs are tough to fit computationally!
#summary(fit_sp)

#fit_mgcv <- gam(cbind(response, size-response) ~ W + s(x,y), data = dat, family = "binomial", method = "ML") #, bs = "gp"
fit_mgcv <- gam(response ~ W + s(x,y), data = dat, family = "poisson", method = "ML") #, bs = "gp"
summary(fit_mgcv)
##-------------------------------


## Proposed method (not really!)
# Example of spatial basis functions first to visualize them  
sp_basisfunctions <- mrts(xy, k = 25, x = xy) %>% 
   as.matrix
colnames(sp_basisfunctions) <- paste0("sp", 1:25)
sp_dat <- data.frame(xy, sp_basisfunctions) %>% 
   pivot_longer(sp1:sp25, names_to = "basisfunction")
sp_dat$basisfunction <- fct_inorder(sp_dat$basisfunction)
ggplot(sp_dat, aes(x = x, y = y, color = value)) +
   geom_point() +
   scale_color_viridis_c() +
   facet_wrap(. ~ basisfunction, nrow = 5) +
   theme_bw()


## Step 1: Fit measurement error model using the resolution adaptive basis functions
sp_basisfunctions <- mrts(xy, k = 50, x = xy)
fit_Wmodel <- autoFRK(Data = matrix(W,ncol=1), loc = xy, G = sp_basisfunctions)
dat$X_pred <- predict(fit_Wmodel)$pred
pairs(dat %>% dplyr::select(X,W,X_pred))
dat %>% 
   dplyr::select(X,W,X_pred) %>% 
   t %>% 
   dist


## Step 2: Using corrected X, fit the spatial model. Here I did not bother to code up fitting resolution adaptive basis functions, but one probably should!
## Also not sure iteration between steps should be done or something fancier; depends on what the marginal log-likelihood maths returns!
#fit_ns_correct <- glm(cbind(response, size-response) ~ X_pred, data = dat, family = "binomial")
fit_ns_correct <- glm(response ~ X_pred, data = dat, family = "poisson")
summary(fit_ns_correct)

#fit_sp_correct <- fitme(cbind(response, size-response) ~ X_pred + Matern(1|x+y), data = dat, family = "binomial") # This can take a while, but it partly demonstrates why spatial (G)LMMs are tough to fit computationally!
#summary(fit_sp_correct)

#fit_mgcv_correct <- gam(cbind(response, size-response) ~ X_pred + s(x,y), data = dat, family = "poisson", method = "ML) #, bs = "gp"
fit_mgcv_correct <- gam(response ~ X_pred + s(x,y), data = dat, family = "poisson", method = "ML") #, bs = "gp"
summary(fit_mgcv_correct)

