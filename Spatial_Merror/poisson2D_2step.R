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
library(FRK)

setwd("D:/Study new/PhD/R code/Spatial_Merror")
#gdbsource("poisson2D_2step.R",interactive=TRUE)
TMB::compile("poisson_2step.cpp","-O1 -g",DLLFLAGS="")
dyn.load(dynlib("poisson_2step"))

true_coefs = c(1,1)
n = 1000

onedataset_sim = function(k0, basisfun_truemod = TRUE , stationary_x = FALSE){
  
  set.seed(134*k0)
  
  #Dataset generation
  
  xy <- data.frame(x = runif(n), y = runif(n))

  
  #Generate autoFRK basis functions
  
  x_basisfunctions <- mrts(xy, k = 30, x = xy) %>% 
    as.matrix
  
  sp_basisfunctions <- mrts(xy, k = 5, x = xy) %>% 
    as.matrix
  sp_basisfunctions <- sp_basisfunctions[,-1] #no intercept
  

  
  
  ## Generate covariates and spatial random effects
  num_X = 1
  sigma_e = sqrt(0.5)
  
  if (basisfun_truemod) {
    
    x_basisfun <- mrts(xy, k = 30, x = xy) %>% 
      as.matrix
    
    alpha_x <- rmvnorm(n = 1, mu = rep(0,ncol(x_basisfun)), Sigma = 0.05*diag(nrow = ncol(x_basisfun))) %>%
      as.vector
    
    X <- x_basisfun%*%alpha_x
    
    rho_basisfun <- mrts(xy, k = 5, x = xy) %>%
      as.matrix
    rho_basisfun_true <- rho_basisfun_true[,-1]
    
    alpha_rho <- rmvnorm(n = 1, mu = rep(0,ncol(rho_basisfun)), Sigma = 0.05*diag(nrow = ncol(rho_basisfun))) %>%
      as.vector
    
    spranef <- rho_basisfun%*%alpha_rho
    
  } else {
    
    if (stationary_x) {
      X <- RFsimulate(model = RMmatern(nu = 1, var = 1, scale = 0.3), x = xy$x, y = xy$y, n = 1)@data[,1]
    } else {
      
      w_1 <- rnorm(1,sd=5)
      w_2 <- rnorm(1,sd=3)

      f_1 <- apply(xy, 1, function(x){cos(pi*norm(x-c(0,1),type = "F")) })
      f_2 <- apply(xy, 1, function(x){cos(2*pi*norm(x-c(3/4,1/4),type = "F")) })

      X <- (f_1*w_1 + f_2*f_2)
      
      # X <- RFsimulate(model = RMmatern(nu = 1, var = 1, scale = 0.3), x = xy$x, y = xy$y, n = 1)@data[,1] +
      #   0.5*xy$x + 0.5*xy$y
      
      # X <- RFsimulate(model = RMmatern(nu = 1, var = 1, scale = 0.3), x = xy$x, y = xy$y, n = 1)@data[,1] +
      #   seq(0,1,length.out = 1000)
    }
    
    spranef <- RFsimulate(model = RMmatern(nu = 0.5, var = 1, scale = 0.5), x = xy$x, y = xy$y, n = 1)@data[,1] 
    
    #spranef <- rep(0,n)
    
    # rho_basisfun <- mrts(xy, k = 5, x = xy) %>%
    #   as.matrix
    # 
    # alpha_rho <- rmvnorm(n = 1, mu = rep(0,ncol(rho_basisfun)), Sigma = diag(nrow = ncol(rho_basisfun))) %>%
    #   as.vector
    # 
    # spranef <- rho_basisfun%*%alpha_rho
    
}
    

  W <- X + rnorm(n, sd = sigma_e )
  
  
  ## Generate response
    
  
  eta <- cbind(1,X) %*% true_coefs + spranef
  simy <- rpois(n, lambda = exp(eta))
  
  
  
  #Fit the model for W via autoFRK()
  fit_Wmodel <- autoFRK(Data = matrix(W,ncol=1), loc = xy, G = x_basisfunctions)
  X_pred <- predict(fit_Wmodel)$pred
  
  
  #Get some 'reasonable' starting values
  
  naive_mod = glm(simy ~ W + sp_basisfunctions , family = poisson)
  no.spatial_mod = glm(simy ~ X_pred , family = poisson)
  
  cholSigma_vec = diag(log(abs(naive_mod$coefficients[-(1:2)])))[lower.tri(diag(nrow =  ncol(sp_basisfunctions)), diag = TRUE)]
  
  
  # cholSigma1_vec = diag(nrow =  ncol(sp_basisfunctions))[lower.tri(diag(nrow =  ncol(sp_basisfunctions)), diag = TRUE)]
  # cholSigma2_vec = diag(nrow =  ncol(x_basisfunctions))[lower.tri(diag(nrow =  ncol(x_basisfunctions)), diag = TRUE)]
  
  
  inputdata=list(sp_basis = sp_basisfunctions , y = as.matrix(simy), predX = X_pred) 
  parameterlist=list(beta = naive_mod$coefficients[1:2] , cholSigma = cholSigma_vec , alpha_rho = naive_mod$coefficients[-(1:2)] )
  
  #naive_mod$coefficients[1:2] instead of true_coefs could be used, but
  
  obj <- MakeADFun(data=inputdata,DLL = "poisson_2step",
                   parameters = parameterlist, random = c("alpha_rho"), silent=T)
  
  
  # fit <- optimr(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS", hessian = TRUE,
  #                       control=list(maxit=500, trace=0))
  
  fit <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr,
                         control=list(trace = 0, iter.max = 1000, eval.max = 1000))


  #Fit using FRK
  
  spdat = xy 
  coordinates(spdat) = ~x+y
  
  # BAUs <- auto_BAUs(manifold = plane(), data = spdat,
  #                   nonconvex_hull = FALSE, cellsize = c(0.03, 0.03), type="grid")
  
  BAUs <- BAUs_from_points(spdat, offset = 1e-10)
  BAUs$fs <- 1 # scalar fine-scale covariance matrix
  FRK_xbasis <- auto_basis(manifold = plane(), data = spdat, nres = 2, type = "bisquare")
  
  spdat = cbind(xy,W)
  coordinates(spdat) = ~x+y
  
  FRK_xfit <- SRE(f = W ~ 1, data = list(spdat), basis = FRK_xbasis, BAUs = BAUs, include_fs = FALSE, K_type = 'precision')
  FRK_xfit <- SRE.fit(FRK_xfit, method = 'TMB')
  pred <- predict(FRK_xfit,newdata=spdat)
  FRKpredX <- pred$newdata$p_mu
  
  spdat = cbind(xy,simy)
  coordinates(spdat) = ~x+y
  BAUs@data <- BAUs@data %>% mutate(FRKpredX = FRKpredX)
  
  FRK_yfit <- SRE(f = simy ~ FRKpredX, data = list(spdat), basis = FRK_xbasis, BAUs = BAUs,
                  include_fs = FALSE, K_type = 'precision', response = "poisson", link = "log")
  FRK_yfit <- SRE.fit(FRK_yfit, method = 'TMB')
  
  out <- list(fit_naive = naive_mod, fit_correct = fit, no.spatial_fit = no.spatial_mod, FRK_fit = FRK_yfit)
  return(out)
}



num_sims = 50
allfits = list()
for (k0 in 1:num_sims) {
  
  skip_to_next <- FALSE
  
  tryCatch(allfits[[k0]] <- onedataset_sim(k0 = k0), error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next } 
}

null_indices = sapply(allfits, function(x) is.null(x))
allfits <- allfits[!null_indices]


##--------------------------
## Process and explore results
##-------------------------
all_naive <- sapply(allfits, function(x) x$fit_naive$coefficients[1:2]) %>%
  t %>%
  as.data.frame

all_correct <- sapply(allfits, function(x) x$fit_correct$par[1:2]) %>%
  t %>%
  as.data.frame

all_no.spatial <- sapply(allfits, function(x) x$no.spatial_fit$coefficients) %>%
  t %>%
  as.data.frame

all_FRK <- sapply(allfits, function(x) x$FRK_fit@alphahat[,1]) %>%
  t %>%
  as.data.frame

colnames(all_naive) <- colnames(all_correct) <- 
  colnames(all_no.spatial) <- colnames(all_FRK) <- c("intercept", "slope")

all_coefs <- rbind(all_naive, all_correct, all_no.spatial, all_FRK) %>%
  mutate(type = rep(c("naive","correct","no_spatial","FRK"), each = num_sims)) %>%
  pivot_longer(intercept:slope)


all_coefs %>% #filter(name=="slope",value>-10) %>%
  ggplot(aes(x = type, y = value)) +
  geom_boxplot(aes(fill = type)) +
  #facet_grid(. ~ name) +
  theme_bw() +
  geom_hline(aes(yintercept = true_coefs),
             data.frame(name = c("intercept","slope"), true_coefs = c(1,1)),linetype="dotted") +
  labs(title = "Boxplots", x = 'Procedure' , y = 'Estimate') +
  scale_fill_brewer(palette = 1) +
  theme(legend.position="right") 
  


all_coefs %>% #filter(name=="slope",value>-10) %>% 
  ggplot(aes(value)) +
  geom_histogram(aes(y=..density..), alpha=0.5, binwidth = 0.25) +
  geom_density(col='red',alpha=.2, fill="#FF6666", adjust = 0.1) +
  facet_grid(type ~ .) +
  #facet_grid(type ~ name) +
  theme_bw() +
  labs(title = "Histograms", x = 'Estimated Value' , y = 'Density') +
  coord_cartesian(xlim = c(-5,5),ylim = c(0,1))



all_coefs %>% #filter(value>-10) %>%
  group_by(name, type) %>%
  summarise(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE))











#make sure you do the following steps AFTER the optim() step.

# rep1 <- sdreport(obj)
# rep1 #MLEs
# summary(rep1) #MLEs and also "estimated random effects"
# as.list(rep1,"Est",FALSE) #easier way to look at the results
# as.list(rep1,"Est",TRUE) #reports Sigma here, because of the ADREPORT(Sigma) line in the .cpp file.
# 
# #for some more options when reporting results, see the help files.
