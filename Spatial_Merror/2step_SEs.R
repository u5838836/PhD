rm(list=ls())
library(TMB)
library(mvtnorm)
library(autoFRK)
library(mgcv)
library(spaMM)
#library(tidyverse)
library(RandomFields)
library(sp)
library(ggplot2)
library(tidyr)
library(forcats)
library(fields)
library(optimx)
library(FRK)
library(Matrix)
library(MASS)
library(dplyr)


setwd("D:/Study new/PhD/R code/Spatial_Merror")
#gdbsource("poisson2D_2step.R",interactive=TRUE)
TMB::compile("poisson_2step.cpp","-O1 -g",DLLFLAGS="")
dyn.load(dynlib("poisson_2step"))

true_coefs = c(3,1)
n = 1024

#xy <- data.frame(x = runif(n), y = runif(n))
xy <- expand.grid(x = ((-sqrt(n)/2):(sqrt(n)/2))[-(sqrt(n)/2+1)], y = ((-sqrt(n)/2):(sqrt(n)/2))[-(sqrt(n)/2+1)])

#Generate autoFRK basis functions for fitting

x_basisfunctions <- mrts(xy, k = 30, x = xy) %>% 
  as.matrix

sp_basisfunctions <- mrts(xy, k = 5, x = xy) %>% 
  as.matrix
sp_basisfunctions <- sp_basisfunctions[,-1] #no intercept

#autoFRK basis functions for the true model

x_basisfun_true <- mrts(xy, k = 30, x = xy) %>% 
  as.matrix

rho_basisfun_true <- mrts(xy, k = 5, x = xy) %>%
  as.matrix
rho_basisfun_true <- rho_basisfun_true[,-1]



onedataset_sim = function(k0, basisfun_truemod = TRUE , stationary_x = TRUE){
  
  #set.seed(73*k0)
  
  #Dataset generation
  
  ## Generate covariates and spatial random effects
  
  
  if (basisfun_truemod) {
    #set.seed(23)
    
    alpha_x <- rmvnorm(n = 1, mu = rep(0,ncol(x_basisfun_true)), Sigma = diag(1/(1:ncol(x_basisfun_true)))) %>%
      as.vector
    
    X <- x_basisfun_true%*%alpha_x
    
    alpha_rho <- rmvnorm(n = 1, mu = rep(0,ncol(rho_basisfun_true)), Sigma = diag(1/(1:ncol(rho_basisfun_true)))) %>%
      as.vector
    
    spranef <- rho_basisfun_true%*%alpha_rho
    
    #set.seed(NULL)
    
  } else {
    
    if (stationary_x) {
      X <- RFsimulate(model = RMmatern(nu = 1, var = 0.5, scale = 5), x = xy$x, y = xy$y, n = 1)@data[,1]
      
      # set.seed(2343)
      # 
      # alpha_x <- rmvnorm(n = 1, mu = rep(0,ncol(x_basisfun_true)), Sigma = diag(1/(1:ncol(x_basisfun_true)))) %>%
      #   as.vector
      # 
      # X <- x_basisfun_true%*%alpha_x
      # 
      # set.seed(NULL)
    } else {
      
      ######-------A few different possible ways to generate (non-stationary) X-------######
      
      # w_1 <- rnorm(1,sd=1)
      # w_2 <- rnorm(1,sd=0.5)
      # 
      # f_1 <- apply(xy, 1, function(x){cos(pi*norm(x-c(0,1),type = "F")) })
      # f_2 <- apply(xy, 1, function(x){cos(2*pi*norm(x-c(3/4,1/4),type = "F")) })
      # 
      # X <- (f_1*w_1 + f_2*w_2)
      
      
      
      # X <- RFsimulate(model = RMmatern(nu = 1, var = 0.5, scale = 5), x = xy$x, y = xy$y, n = 1)@data[,1] +
      #   0.01*xy$x + 0.01*xy$y
      
      # X <- RFsimulate(model = RMmatern(nu = 1, var = 1, scale = 0.3), x = xy$x, y = xy$y, n = 1)@data[,1] +
      #   seq(0,1,length.out = 1024)
      
      
      
      #X <- 2*exp(-(xy$x^2+xy$y^2)/sqrt(n))
      
      u = (xy$x + 0.5*sqrt(n))/sqrt(n)
      v = (xy$y + 0.5*sqrt(n))/sqrt(n)
      X = 0.5*(1.35 + exp(u)*13*(u - 0.6)^2*exp(-v)*sin(7*v))
      #X = u*sin(4*pi*v) + 1
      # X = 2*(pi**0.3*0.4)*(1.2*exp(-(u-0.2)^2/0.3^2-(v-0.3)^2/0.4^2)+
      #                0.8*exp(-(u-0.7)^2/0.3^2-(v-0.8)^2/0.4^2))
      
      
      # u = 10*xy$x/sqrt(n)
      # v = 10*xy$y/sqrt(n)
      # X = 5*exp( -2.5*((u-2)^2 + (v-2)^2 ) ) + exp(-(u^2 + v^2)/3)
      
    }
    
    spranef <- RFsimulate(model = RMmatern(nu = 0.5, var = 0.5, scale = 1), x = xy$x, y = xy$y, n = 1)@data[,1] 
    
    #spranef <- rep(0,n)
    
    # set.seed(2345)
    # 
    # alpha_rho <- rmvnorm(n = 1, mu = rep(0,ncol(rho_basisfun_true)), Sigma = diag(1/(1:ncol(rho_basisfun_true)))) %>%
    #   as.vector
    # 
    # spranef <- rho_basisfun_true%*%alpha_rho
    # 
    # set.seed(NULL)
    
    #spranef <- 2*exp(-(xy$x^2+xy$y^2)/sqrt(n))
    
    # u = xy$x
    # v = xy$y
    
    # u = (xy$x + 0.5*sqrt(n))/sqrt(n)
    # v = (xy$y + 0.5*sqrt(n))/sqrt(n)
    # spranef = (exp(u)*(u - 0.6)^2*exp(-v)*sin(7*v))
    
    
    # spranef = u*sin(4*pi*v)
    # spranef = 2*(pi**0.3*0.4)*(1.2*exp(-(u-0.2)^2/0.3^2-(v-0.3)^2/0.4^2)+
    #                0.8*exp(-(u-0.7)^2/0.3^2-(v-0.8)^2/0.4^2))
    
  }
  
  sigma_e = sqrt(0.5)
  W <- X + rnorm(n, sd = sigma_e )
  
  #Fit the model for W via autoFRK()
  fit_Wmodel <- autoFRK(Data = matrix(W,ncol=1), loc = xy, G = x_basisfunctions)
  X_pred <- predict(fit_Wmodel)$pred
  
  
  ## Generate response
  
  eta <- cbind(1,X) %*% true_coefs + spranef
  simy <- rpois(n, lambda = exp(eta))
  
  
  #Get some 'reasonable' starting values
  
  naive_mod = glm(simy ~ W + sp_basisfunctions , family = poisson)
  no.spatial_mod = glm(simy ~ X_pred , family = poisson)
  
  cholSigma_vec = diag(log(abs(naive_mod$coefficients[-(1:2)])))[lower.tri(diag(nrow =  ncol(sp_basisfunctions)), diag = TRUE)]
  
  
  inputdata=list(sp_basis = sp_basisfunctions , y = as.matrix(simy), predX = X_pred) 
  parameterlist=list(beta = naive_mod$coefficients[1:2] , cholSigma = cholSigma_vec , alpha_rho = naive_mod$coefficients[-(1:2)] )
  
  
  obj <- MakeADFun(data=inputdata,DLL = "poisson_2step",
                   parameters = parameterlist, random = c("alpha_rho"), silent=T)
  
  
  # fit <- optimr(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS", hessian = TRUE,
  #                       control=list(maxit=500, trace=0))
  
  fit <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr,
                control=list(trace = 0, iter.max = 1000, eval.max = 1000))
  
  
  ###-------Get SEs--------###
  
  #Examine outputted SEs from TMB; *should* be correct when model is correct and true X is used
  
  report = sdreport(obj,getJointPrecision = TRUE, getReportCovariance = FALSE, skip.delta.method = TRUE)
  report_summary=summary(report)
  params = report_summary %>% 
    as.data.frame %>% 
    filter(rownames(report_summary)%in%c("alpha_rho","beta")) %>%
    select("Estimate") %>% .$Estimate
  
  sigmarho_hat = as.list(report,"Est",TRUE)$Sigma
  
  #Non-parametric estimate of variance of score
  
  sigmax_hat <- fit_Wmodel$M
  alphax_hat <- fit_Wmodel$w
  phi_hat <- fit_Wmodel$s
  #crossprod(x_basisfunctions,W - X_pred)/phi_hat - ginv(sigmax_hat)%*%alphax_hat
  varscore_x = 0
  varscore_y = 0
  
  eta = params[1] + params[2]*X_pred + sp_basisfunctions%*%params[-(1:2)] 
  iterative_weights = poisson()$var(exp(eta)) %>% as.vector %>% diag.spam
  
  for (j in 1:n) {
    varscore_x = varscore_x + tcrossprod(x_basisfunctions[j,]*(W[j] - X_pred[j])/phi_hat ) # - ginv(sigmax_hat)%*%alphax_hat/n
    
    score_y0 = simy[j] - exp(eta[j]) #-
      # 0.5*sum(diag(solve(crossprod(sp_basisfunctions,iterative_weights)%*%sp_basisfunctions +
      #             ginv(sigmarho_hat))%*%
      # crossprod(sp_basisfunctions,iterative_weights)%*%sp_basisfunctions))/n
    
    score_y1 = X_pred[j]*(simy[j] - exp(eta[j])) #- 
      # 0.5*sum(diag(solve(crossprod(sp_basisfunctions,diag.spam(as.vector(X_pred))%*%iterative_weights)%*%
      #                      sp_basisfunctions + ginv(sigmarho_hat))%*%
      #                crossprod(sp_basisfunctions,iterative_weights)%*%sp_basisfunctions))/n
    
    score_y_rho = sp_basisfunctions[j,]*(simy[j] - exp(eta[j])) # - ginv(sigmarho_hat)%*%params[-(1:2)]/n
    score_y = c(score_y0,score_y1,score_y_rho)
    varscore_y = varscore_y + tcrossprod(score_y)
  }
  
  ###adjustment terms for the log determinant and penalty term
  
  # varscore_y[1:2,1:2] = varscore_y[1:2,1:2] +
  #   tcrossprod(c( -0.5*sum(diag(solve(crossprod(sp_basisfunctions,iterative_weights)%*%sp_basisfunctions +
  #                                     ginv(sigmarho_hat))%*%
  #                               crossprod(sp_basisfunctions,iterative_weights)%*%sp_basisfunctions)) ,
  #               - 0.5*sum(diag(solve(crossprod(sp_basisfunctions,diag.spam(as.vector(X_pred))%*%iterative_weights)%*%
  #                                      sp_basisfunctions + ginv(sigmarho_hat))%*%
  #                                crossprod(sp_basisfunctions,iterative_weights)%*%sp_basisfunctions))))


  varscore_y[-(1:2),-(1:2)] = varscore_y[-(1:2),-(1:2)] + ginv(sigmarho_hat,10*(.Machine$double.eps)^(1/3))
  varscore_x = varscore_x + ginv(sigmax_hat)
  
  #observed negative Hessians (joint)
  
  Xstar = cbind(1,X_pred,sp_basisfunctions)
  Y_hess = crossprod(Xstar,iterative_weights)%*%Xstar + 
    bdiag.spam(matrix(0,nrow=2,ncol=2),ginv(sigmarho_hat,10*(.Machine$double.eps)^(1/3)))
  
  X_hess = crossprod(x_basisfunctions)/phi_hat + ginv(sigmax_hat) #X part/block
  XY_int_hess = fit$par[2]*apply(iterative_weights%*%x_basisfunctions,2,sum) 
  XY_slope_hess = fit$par[2]*apply(as.vector(X_pred)*iterative_weights%*%x_basisfunctions,2,sum) #-
    #as.vector(crossprod(x_basisfunctions, simy - exp(eta)))
  
  XY_rho_hess = NULL
  for (k in 1:ncol(sp_basisfunctions)) {
    XY_rho_hess = rbind(XY_rho_hess,
                        fit$par[2]*apply(sp_basisfunctions[,k]*iterative_weights%*%x_basisfunctions,2,sum)  )
  }
  
  #alpha_index = nrow(report$jointPrecision) - ncol(sp_basisfunctions) +1
  #jointPrecision = bdiag(report$jointPrecision[c(1:2,alpha_index:nrow(report$jointPrecision)),c(1:2,alpha_index:nrow(report$jointPrecision))], X_hess)
  
  jointPrecision = bdiag(Y_hess, X_hess)
  jointPrecision[1:(2+ncol(sp_basisfunctions)),(3+ncol(sp_basisfunctions)):ncol(jointPrecision)] =
    rbind(XY_int_hess,XY_slope_hess,XY_rho_hess)
  
  out <- list(fit_naive = naive_mod, fit_correct = fit, no.spatial_fit = no.spatial_mod, report = report,
              varscore = bdiag.spam(varscore_y,varscore_x), jointPrecision = jointPrecision)
  return(out)
}



num_sims = 100
# allfits = list()
# for (k0 in 1:num_sims) {
#   
#   skip_to_next <- FALSE
#   
#   tryCatch(allfits[[k0]] <- onedataset_sim(k0 = k0), error = function(e) { skip_to_next <<- TRUE})
#   
#   if(skip_to_next) { next } 
#   
#   message(paste0("Finished Dataset ",k0))
# }


tic <- proc.time()
allfits <- foreach(i=1:num_sims,.packages = c("MASS","mvtnorm",
                                              "doParallel","glm2","GGally","Matrix","TMB","dplyr", "autoFRK",
                                              "mgcv","spaMM","RandomFields","sp","fields")) %dopar% {
                                                
                                                setwd("D:/Study new/PhD/R code/Spatial_Merror")        
                                                dyn.load(dynlib("poisson_2step"))
                                                skip_to_next <- FALSE
                                                tryCatch(result <- onedataset_sim(k0 = i), error = function(e) { skip_to_next <<- TRUE})
                                                #result <- onedataset_sim(k0 = i)
                                                if(skip_to_next) { return(NULL) } else { return(result) }
                                              }
toc <- proc.time()


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


colnames(all_naive) <- colnames(all_correct) <- 
  colnames(all_no.spatial) <- c("intercept", "slope")

all_coefs <- rbind(all_naive, all_correct, all_no.spatial) %>%
  mutate(type = rep(c("naive","correct","no_spatial"), each = num_sims - sum(null_indices))) %>%
  pivot_longer(intercept:slope)


all_coefs %>% filter(name == "slope") %>%
  ggplot(aes(x = type, y = value)) +
  geom_boxplot(aes(fill = type)) +
  #facet_grid(. ~ name) +
  theme_bw() +
  geom_hline(aes(yintercept = true_coefs),
             data.frame(name = c("intercept","slope"), true_coefs = c(1,1)),linetype="dotted") +
  labs(title = "Boxplots", x = 'Procedure' , y = 'Estimate') +
  scale_fill_brewer(palette = 1) +
  theme(legend.position="right") 



all_coefs %>% filter(name == "slope") %>%
  ggplot(aes(value)) +
  geom_histogram(aes(y=..density..), alpha=0.5, binwidth = 0.25) +
  geom_density(col='red',alpha=.2, fill="#FF6666", adjust = 3) +
  facet_grid(type ~ .) +
  #facet_grid(type ~ name) +
  theme_bw() +
  labs(title = "Histograms", x = 'Estimated Value' , y = 'Density') +
  coord_cartesian(xlim = c(-1,3),ylim = c(0,1))



all_coefs %>% 
  group_by(name, type) %>%
  summarise(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE))


sapply(allfits, function(x) summary(x$report)[2,2]) %>% mean(na.rm = TRUE)

sapply(allfits, function(x){
  invnegHess = as.matrix(solve(x$jointPrecision))
  tcrossprod(invnegHess%*%x$varscore,invnegHess)[2,2]
    }) %>% mean(na.rm = TRUE) %>% sqrt


sapply(allfits, function(x){
  invnegHess = as.matrix(solve(x$jointPrecision))
  tcrossprod(invnegHess%*%x$varscore,invnegHess)[2,2]
}) %>% as.data.frame %>%
  ggplot(aes(x=.)) +
  geom_histogram(aes(y=..density..), alpha=0.5, bins = 20, fill = 'blue') +
  labs(title = "Histogram", x = 'Estimated Value' , y = 'Density') +
  theme_bw()


stopCluster(cl)

#make sure you do the following steps AFTER the optim() step.

# rep1 <- sdreport(obj)
# rep1 #MLEs
# summary(rep1) #MLEs and also "estimated random effects"
# as.list(rep1,"Est",FALSE) #easier way to look at the results
# as.list(rep1,"Est",TRUE) #reports Sigma here, because of the ADREPORT(Sigma) line in the .cpp file.
# 
# #for some more options when reporting results, see the help files.


