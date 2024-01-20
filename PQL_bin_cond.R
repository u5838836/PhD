rm(list = ls())
library(glmmTMB)
library(dplyr)
library(xtable)
library(mvtnorm)
library(MASS)
library(Matrix)
library(parallel)
library(foreach)
library(doParallel)
library(GGally)
library(lme4)
library(glm2)
library(doMC)
registerDoMC(cores=28)
set.seed(259)

#The number of clusters and cluster sizes we wish to explore for the overall sim
clus_test = c(25,50,100,200)
timepoint_test = c(25,50,100,200)

#create storage for quantities of interest

fixedcoverage1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
fixedcoverage2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
randomcoverage1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
randomcoverage2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))

fixedPQLshapiro1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
fixedPQLshapiro2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
randomPQLshapiro1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
randomPQLshapiro2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))


######------Function for simulating and fitting one dataset-------########
onedatasetsim <- function(i) {     
  
  ##--------------
  ## Simulate data
  ##--------------
  
  # Generate the observations
  simdat <- data.frame(ID = factor(rep(1:num_clus, each = num_resp*num_timepoints)), 
                       time = factor(rep(rep(1:num_timepoints, num_resp), num_clus)),
                       resp = factor(rep(rep(1:num_resp, each = num_timepoints), num_clus)))
  
  # Generate covariate 
  set.seed(253)
  bigX <- lapply(1:num_clus, function(j) kronecker(diag(num_resp), cbind(1,rmvnorm(num_timepoints, mean = rep(0,num_X),sigma = 1*diag(num_X)))) ) %>%
    do.call(rbind, .)
  bigZ <- bigX[,1:(num_Z),drop=FALSE]
  set.seed(NULL)
  
  #(True) linear predictor for each observation and simulated observations
  
  simdat$eta <- c(bigX %*% c(t(true_beta)) + rowSums(bigZ * true_alpha[simdat$ID,]))
  simdat$y <- rbinom(n = nrow(simdat), size=1, prob = exp(simdat$eta)/(1+exp(simdat$eta)))
  
  ##--------------
  ## Do some univariate GLMM fits first via Laplace.
  ##--------------
  bigX <- as.matrix(bigX)
  bigZ <- as.matrix(bigZ)
  smallX <- bigX[,1:(num_X+1)] 
  smallX <- smallX[rowSums(smallX) != 0,] 
  smallZ <- bigZ[,1:num_Z]
  smallZ <- smallZ[rowSums(smallZ) != 0,] 
  
  fit_resps=list()
  for (j in 1:num_resp) {
    fit_resps[[j]] <- glmmTMB(formula = y ~  smallX - 1 + (smallZ-1|ID), family = "binomial",data = subset(simdat, resp == j), 
                              map = list(theta = factor(rep(NA,length(chol_start)))), start = list(theta = chol_start))
  }

  ##--------------
  ## Univariate fits via PQL ######
  ##--------------     
  univariate_glmmPQL <- function(starting_fit, resp, max_iter = 10000) {        
    
    #set starting values for parameters
    cw_beta <- (fixef(starting_fit))$cond
    cw_G <- working_G
    cw_alpha_mat <- as.matrix(ranef(starting_fit)$cond$ID)
    
    ## Optimization parameters - error tolerance
    err <- 100
    counter <- 0
    num_Xplusint <- num_X+1
    while(err > 1e-3 & counter < max_iter) {  
      
      cw_Ginv <- solve(cw_G)
      ## Maximise wrt alpha -- on a per cluster basis
      update_alpha <- function(j) {
        
        selindices <- which(simdat$ID == j & simdat$resp==resp)
        subX <- bigX[selindices,(num_Xplusint*resp - num_Xplusint + 1):(num_Xplusint*resp),drop = FALSE]
        subZ <- bigZ[selindices,(num_Z*resp - num_Z + 1):(num_Z*resp),drop = FALSE]
        suby <- simdat$y[selindices]
        
        innerlogL_alpha <- function(alpha) {
          cw_eta <- as.vector(subX %*% cw_beta + subZ %*% alpha)           
          likcontrib <- -sum(dbinom(suby, size=1, prob = exp(cw_eta)/(1+exp(cw_eta)), log = TRUE)) + 0.5*crossprod(alpha, cw_Ginv) %*% alpha
          return(as.vector(likcontrib))
        }
        
        innergrad_alpha <- function(alpha) {
          cw_eta <- as.vector(subX %*% cw_beta + subZ %*% alpha)           
          out <- -crossprod(subZ, suby - exp(cw_eta)/(1+exp(cw_eta))) + cw_Ginv %*% alpha     
          return(as.vector(out))
        }
        
        do_update <- optim(cw_alpha_mat[j,], fn = innerlogL_alpha, gr = innergrad_alpha, control = list(trace = 0, maxit = max_iter), method = "BFGS")
        return(do_update$par)
      }     
      new_alpha_mat <- foreach(j = 1:num_clus, .combine = "rbind") %do% update_alpha(j = j)
      
      ## Maximise wrt beta
      subX <- as.matrix(bigX[which(simdat$resp==resp),(num_Xplusint*resp - num_Xplusint + 1):(num_Xplusint*resp),drop = FALSE])
      subZ <- as.matrix(bigZ[which(simdat$resp==resp),(num_Z*resp - num_Z + 1):(num_Z*resp),drop = FALSE])
      suby <- simdat$y[which(simdat$resp==resp)]
      make_offset <- rowSums(subZ * new_alpha_mat[simdat$ID[simdat$resp == resp],])
      update_beta <- glm2(suby ~ subX - 1, offset = make_offset, family = "binomial")
      new_beta <- update_beta$coefficients
      rm(update_beta, suby)
      
      ## G is known
      new_G = working_G
      
      #error for this iteration
      err <- sum((new_beta - cw_beta)^2) + 0.5*sum((new_G - cw_G)^2) + sum((new_alpha_mat - cw_alpha_mat)^2)
      
      #update current estimates
      cw_beta <- new_beta
      cw_G <- as.matrix(new_G)
      cw_alpha_mat <- new_alpha_mat
      counter <- counter + 1
    }
    
    out <- list(beta = new_beta, G = new_G, alpha = new_alpha_mat)
    return(out)
  }
  
  #fit the glmm via PQL
  unifits=list()
  for (j in 1:num_resp) {
    unifits[[j]] <- univariate_glmmPQL(starting_fit = fit_resps[[j]], resp = j)
  }
  
  out <- list(unifits = unifits,true_alphas = true_alpha,lme4fit = fit_resps[[1]]) 
  return(out)
}

#######-------Running the sims------#######

for (r in 1:length(clus_test)) {
  for (s in 1:length(timepoint_test)) {
    num_resp <- 1
    num_clus <- clus_test[r]
    num_timepoints <- timepoint_test[s]
    
    num_X <- 1 # Number of fixed effects, excludes intercept
    num_Z <- 2
    num_sims = 1000
    
    #Set true values
    true_beta = c(-0.1,0.1)
    true_G = 1*diag(2)
    working_G = 1*diag(2)
    
    #Cholesky decomposition of working G to use in glmmTMB()
    working_G_chol=t(chol(cov2cor(working_G)))
    working_G_chol=diag(diag(1/working_G_chol))%*%working_G_chol
    chol_start = log(sqrt(diag(working_G)))
    chol_start = c(chol_start,working_G_chol[lower.tri(working_G_chol)])
    
    # Generate the 'random' effects outside foreach so they are conditioned on. Reparametrise to satisfy the STZ constraint.
    true_alpha <- rmvnorm(num_clus, sigma = true_G)
    mean_true_alpha = apply(true_alpha,2,mean)
    true_beta = true_beta + mean_true_alpha
    true_alpha = sweep(true_alpha,2,mean_true_alpha)
    
    
    tic <- proc.time()
    results <- foreach(i=1:num_sims,.packages = c("MASS","mvtnorm","doParallel","glm2","GGally","Matrix","glmmTMB","dplyr")) %dopar% onedatasetsim(i = i)
    toc <- proc.time()
    
    ###extract the 1000 prediction gaps and estimates for the first cluster
    alphadiffs = rep(0,num_sims)
    alphas = rep(0,num_sims)
    alphadiffs1 = rep(0,num_sims)
    alphas1 = rep(0,num_sims)
    
    
    for (i in 1:num_sims) {
      univpreds = results[[i]]$unifits[[1]]$alpha
      alphadiffs[i] = (univpreds - results[[i]]$true_alphas)[1,1]
      alphas[i] = univpreds[1,1]
      alphadiffs1[i] = (univpreds - results[[i]]$true_alphas)[1,2]
      alphas1[i] = univpreds[1,2]
    }
    
    
    ###extract fixed effects estimates
    
    betadiffs=rep(0,num_sims)
    betadiffs1=rep(0,num_sims)
    
    for (i in 1:num_sims) {
      univbetas = results[[i]]$unifits[[1]]$beta
      betadiffs[i] = (univbetas - true_beta)[1]
      betadiffs1[i] = (univbetas - true_beta)[2]
    }
    
    #intercepts - shapiro-wilk p-values
    fixedPQLshapiro1[r,s] = shapiro.test(betadiffs)$p.value
    randomPQLshapiro1[r,s] = shapiro.test(alphadiffs)$p.value
    
    #slopes - shapiro-wilk p-values
    fixedPQLshapiro2[r,s] = shapiro.test(betadiffs1)$p.value
    randomPQLshapiro2[r,s] = shapiro.test(alphadiffs1)$p.value
    
    # Generate same covariates
    set.seed(253)
    bigX <- lapply(1:num_clus, function(i) kronecker(diag(num_resp), cbind(1,rmvnorm(num_timepoints, mean = rep(0,num_X),sigma = 1*diag(num_X)))) ) %>%
      do.call(rbind, .)
    bigZ <- bigX[,1:(num_Z),drop=FALSE]
    set.seed(NULL)
    
    # Generate true eta and use this to calculate the appropriate variances to use in building coverage intervals for fixed and random effects
    simdat <- data.frame(ID = factor(rep(1:num_clus, each = num_resp*num_timepoints)), 
                         time = factor(rep(rep(1:num_timepoints, num_resp), num_clus)),
                         resp = factor(rep(rep(1:num_resp, each = num_timepoints), num_clus)))
    
    eta <- c(bigX %*% c(t(true_beta)) + rowSums(bigZ * true_alpha[simdat$ID,]))
    
    #variance for fixed effects
    betavar=0
    for (i in 0:(num_clus-1)) {
      
      p = exp(eta[(i*num_timepoints+1):((i+1)*num_timepoints)])/(1+exp(eta[(i*num_timepoints+1):((i+1)*num_timepoints)]))
      
      betavar = betavar + solve(crossprod(bigX[(i*num_timepoints+1):((i+1)*num_timepoints),],
                                          diag(p*(1-p)))%*% 
                                  bigX[(i*num_timepoints+1):((i+1)*num_timepoints),]/num_timepoints)
    }
    betavar = betavar/num_clus
    
    #variance for random effects
    p = exp(eta[1:num_timepoints])/(1+exp(eta[1:num_timepoints]))
    alphavar = solve(crossprod(bigZ[1:num_timepoints,],
                               diag(p*(1-p)))%*%bigZ[1:num_timepoints,]/num_timepoints)
    
    
    #####Nominal 95% Coverage intervals and coverage probability#####
    
    fixedcount=0
    fixedcount1=0
    randcount=0
    randcount1=0
    
    for (i in 1:num_sims) {
      fixed_CI = c(results[[i]]$unifits[[1]]$beta[1] - 1.96*sqrt(betavar[1,1]/(num_timepoints*num_clus)), 
                   results[[i]]$unifits[[1]]$beta[1] + 1.96*sqrt(betavar[1,1]/(num_timepoints*num_clus)))
      if (true_beta[1]<fixed_CI[2]&true_beta[1]>fixed_CI[1]) {
        fixedcount=fixedcount+1
      }
      
      fixed_CI1 = c(results[[i]]$unifits[[1]]$beta[2] - 1.96*sqrt(betavar[2,2]/(num_timepoints*num_clus)), 
                    results[[i]]$unifits[[1]]$beta[2] + 1.96*sqrt(betavar[2,2]/(num_timepoints*num_clus)))
      if (true_beta[2]<fixed_CI1[2]&true_beta[2]>fixed_CI1[1]) {
        fixedcount1=fixedcount1+1
      }
      
      
      rand_CI = c(results[[i]]$unifits[[1]]$alpha[1,1] - 1.96*sqrt(alphavar[1,1]/num_timepoints), 
                  results[[i]]$unifits[[1]]$alpha[1,1] + 1.96*sqrt(alphavar[1,1]/num_timepoints))
      if (true_alpha[1,1]<rand_CI[2]&true_alpha[1,1]>rand_CI[1]) {
        randcount=randcount+1
      }
      
      rand_CI1 = c(results[[i]]$unifits[[1]]$alpha[1,2] - 1.96*sqrt(alphavar[2,2]/num_timepoints), 
                   results[[i]]$unifits[[1]]$alpha[1,2] + 1.96*sqrt(alphavar[2,2]/num_timepoints))
      if (true_alpha[1,2]<rand_CI1[2]&true_alpha[1,2]>rand_CI1[1]) {
        randcount1=randcount1+1
      }
    }
    
    fixedcoverage1[r,s] = fixedcount/num_sims
    fixedcoverage2[r,s] = fixedcount1/num_sims
    randomcoverage1[r,s] = randcount/num_sims
    randomcoverage2[r,s] = randcount1/num_sims
  }
}

#plots

fixedcoverage1 %>% c %>% cbind(rep(c(25,50,100,200),4)) %>% cbind(rep(c(25,50,100,200),each=4)) %>%
  cbind(rep("fixed_int")) %>% as.data.frame -> fixedcoverage1
colnames(fixedcoverage1) <- c("coverage","num_clus","clus_size","quantity")

fixedcoverage2 %>% c %>% cbind(rep(c(25,50,100,200),4)) %>% cbind(rep(c(25,50,100,200),each=4)) %>%
  cbind(rep("fixed_slope")) %>% as.data.frame -> fixedcoverage2
colnames(fixedcoverage2) <- c("coverage","num_clus","clus_size","quantity")

randomcoverage1 %>% c %>% cbind(rep(c(25,50,100,200),4)) %>% cbind(rep(c(25,50,100,200),each=4)) %>%
  cbind(rep("rand_int")) %>% as.data.frame -> randomcoverage1
colnames(randomcoverage1) <- c("coverage","num_clus","clus_size","quantity")

randomcoverage2 %>% c %>% cbind(rep(c(25,50,100,200),4)) %>% cbind(rep(c(25,50,100,200),each=4)) %>%
  cbind(rep("rand_slope")) %>% as.data.frame -> randomcoverage2
colnames(randomcoverage2) <- c("coverage","num_clus","clus_size","quantity")

coverages = rbind(fixedcoverage1,fixedcoverage2,randomcoverage1,randomcoverage2)
coverages[,-4] =  apply(coverages[,-4],c(1,2),as.numeric)

quantity.labs <- c("Fixed Intercept", "Fixed Slope", "Random Intercept", "Random Slope")
names(quantity.labs) <- c("fixed_int", "fixed_slope", "rand_int", "rand_slope")

p1 = coverages %>% ggplot(aes(x = clus_size, y = coverage)) +
  geom_line(aes(col = as.factor(num_clus))) +
  facet_wrap(~ quantity, labeller = labeller(quantity = quantity.labs)) +
  theme_bw() +
  labs(x = 'Cluster size' , y = 'Empirical Coverage', col = 'Number of Clusters') + 
  geom_hline(yintercept = 0.95, lty = 2) +
  coord_cartesian(ylim = c(0,1)) +
  scale_colour_viridis_d()


fixedPQLshapiro1 %>% c %>% cbind(rep(c(25,50,100,200),4)) %>% cbind(rep(c(25,50,100,200),each=4)) %>%
  cbind(rep("fixed_int")) %>% as.data.frame -> fixedPQLshapiro1
colnames(fixedPQLshapiro1) <- c("PQLshapiro","num_clus","clus_size","quantity")

fixedPQLshapiro2 %>% c %>% cbind(rep(c(25,50,100,200),4)) %>% cbind(rep(c(25,50,100,200),each=4)) %>%
  cbind(rep("fixed_slope")) %>% as.data.frame -> fixedPQLshapiro2
colnames(fixedPQLshapiro2) <- c("PQLshapiro","num_clus","clus_size","quantity")

randomPQLshapiro1 %>% c %>% cbind(rep(c(25,50,100,200),4)) %>% cbind(rep(c(25,50,100,200),each=4)) %>%
  cbind(rep("rand_int")) %>% as.data.frame -> randomPQLshapiro1
colnames(randomPQLshapiro1) <- c("PQLshapiro","num_clus","clus_size","quantity")

randomPQLshapiro2 %>% c %>% cbind(rep(c(25,50,100,200),4)) %>% cbind(rep(c(25,50,100,200),each=4)) %>%
  cbind(rep("rand_slope")) %>% as.data.frame -> randomPQLshapiro2
colnames(randomPQLshapiro2) <- c("PQLshapiro","num_clus","clus_size","quantity")

PQLshapiros = rbind(fixedPQLshapiro1,fixedPQLshapiro2,randomPQLshapiro1,randomPQLshapiro2)
PQLshapiros[,-4] =  apply(PQLshapiros[,-4],c(1,2),as.numeric)

quantity.labs <- c("Fixed Intercept", "Fixed Slope", "Random Intercept", "Random Slope")
names(quantity.labs) <- c("fixed_int", "fixed_slope", "rand_int", "rand_slope")

p2 = PQLshapiros %>% ggplot(aes(x = clus_size, y = PQLshapiro)) +
  geom_line(aes(col = as.factor(num_clus))) +
  facet_wrap(~ quantity, labeller = labeller(quantity = quantity.labs)) +
  theme_bw() +
  labs(x = 'Cluster size' , y = 'P-value', col = 'Number of Clusters') + 
  geom_hline(yintercept = 0.05, lty = 2) +
  coord_cartesian(ylim = c(0,1)) +
  scale_colour_viridis_d()

p1
ggsave("cond_bin_coverage.pdf", width = 7, height = 7)
p2
ggsave("cond_bin_shapiro.pdf", width = 7, height = 7)

save(coverages,PQLshapiros, file = "cond_bin.Rdata")







