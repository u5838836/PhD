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
set.seed(261)

#The number of clusters and cluster sizes we wish to explore for the overall sim
clus_test = c(25,50,100,200)
timepoint_test = c(25,50,100,200)

#create storage for quantities of interest

fixedcoverage1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
fixedcoverage2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
randomcoverage1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
randomcoverage2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
lcombcoverage1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
lcombcoverage2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))

fixedPQLshapiro1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
fixedPQLshapiro2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
randomPQLshapiro1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
randomPQLshapiro2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
diffsPQLshapiro1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
diffsPQLshapiro2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))

corrected_randomcoverage1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
corrected_randomcoverage2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))

frobenius_norm = matrix(nrow=length(clus_test),ncol=length(timepoint_test))

########--------Define the sim function for one dataset--------########

onedatasetsim <- function(i) {     
  ##--------------
  ## Simulate data
  ##--------------
  # Generate the random effects
  true_alpha <- rmvnorm(num_clus, sigma = true_G)
  
  # Generate the observations
  simdat <- data.frame(ID = factor(rep(1:num_clus, each = num_resp*num_timepoints)), 
                       time = factor(rep(rep(1:num_timepoints, num_resp), num_clus)),
                       resp = factor(rep(rep(1:num_resp, each = num_timepoints), num_clus)))
  
  # Generate covariate -- assuming all responses have the same set of covariates, which is usually the case
  set.seed(2530)
  bigX <- lapply(1:num_clus, function(i) kronecker(diag(num_resp), cbind(1,rmvnorm(num_timepoints, mean = rep(0,num_X),sigma = 1*diag(num_X)))) ) %>%
    do.call(rbind, .)
  bigZ <- bigX[,1:(num_Z),drop=FALSE]
  set.seed(NULL)
  
  #(True) linear predictor for each observation and simulated observations
  simdat$eta <- c(bigX %*% c(t(true_beta)) + rowSums(bigZ * true_alpha[simdat$ID,]))
  simdat$y <- rbinom(n = nrow(simdat), size = 1, prob = exp(simdat$eta)/(1+exp(simdat$eta)))
  
  ##--------------
  ## Do some univariate GLMM fits first via Laplace
  ##--------------
  bigX <- as.matrix(bigX)
  bigZ <- as.matrix(bigZ)
  smallX <- bigX[,1:(num_X+1)] 
  smallX <- smallX[rowSums(smallX) != 0,] 
  smallZ <- bigZ[,1:num_Z]
  smallZ <- smallZ[rowSums(smallZ) != 0,] 
  
  simdat = cbind(smallX, smallZ, simdat)
  colnames(simdat)[1:4] = c('x_int','x_slope','z_int','z_slope')
  
  fit_resps=list()
  for (j in 1:num_resp) {
    fit_resps[[j]] <- glmmTMB(formula = y ~  x_int + x_slope -1 + (z_int + z_slope - 1|ID) , family = "binomial", data = subset(simdat, resp == j), 
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
      
      ## Maximise wrt G using fixed point iterative update
      bigW <- as.vector(binomial()$var(binomial()$linkinv(subX %*% new_beta + make_offset)))       
      cw_G0 <- cw_G
      err2 <- 1 
      while(err2 > 0.01) {
        new_G <- sapply(1:num_clus, function(j) {
          selindices <- (num_timepoints*(j-1)):(num_timepoints*j)
          solve(crossprod(as.matrix(subZ[selindices,])*sqrt(bigW[selindices])) + solve(cw_G0)) + tcrossprod(new_alpha_mat[j,]) 
        })
        new_G <- matrix(rowMeans(new_G), nrow = num_Z)
        err2 <- norm(new_G - cw_G0, type = "F")
        cw_G0 <- new_G
      }
      
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
  
  unifits=list()
  for (j in 1:num_resp) {
    unifits[[j]] <- univariate_glmmPQL(starting_fit = fit_resps[[j]], resp = j)
  }
  
  out <- list(unifits = unifits,true_alphas = true_alpha,lme4fit = fit_resps[[1]]) 
  return(out)
}


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
    
    
    
    #### Get appropriate normal scale mixture quantiles via simulation - n.b. prediction gap depends on rate of n/m ####
    
    set.seed(2530)
    bigX <- lapply(1:num_clus, function(i) kronecker(diag(num_resp), cbind(1,rmvnorm(num_timepoints, mean = rep(0,num_X),sigma = 1*diag(num_X)))) ) %>%
      do.call(rbind, .)
    bigZ <- bigX[,1:(num_Z),drop=FALSE]
    
    smallX = bigX[1:num_timepoints,]
    smallZ = smallX
    
    Ss = matrix(nrow=10000,ncol=num_Z)
    SSS = matrix(nrow=10000,ncol=num_Z)
    
    if (num_clus<=num_timepoints) {
      
      for (i in 1:10000) {
        bs=rmvnorm(num_clus,sigma = true_G)
        b1=bs[1, , drop = FALSE]
        
        eta <- c(smallX %*% c(t(true_beta)) + smallZ %*% t(b1))
        p = exp(eta)/(1+exp(eta))
        
        alphavar = solve(crossprod(smallZ,
                                   bdiag(diag(p*(1-p))))%*%smallZ/num_timepoints)
        
        S = rmvnorm(1, mean = apply(bs, 2, sum)/sqrt(num_clus) , sigma = as.matrix(alphavar) )
        Ss[i,] = S
        SSS[i,] = rmvnorm(1 , sigma = as.matrix(alphavar) )
      }
      
    } else {
      if (num_clus>num_timepoints) {
        for (i in 1:10000) {
          b1=rmvnorm(1,sigma = true_G)
          
          eta <- c(smallX %*% c(t(true_beta)) + smallZ %*% t(b1))
          p = exp(eta)/(1+exp(eta))
          
          alphavar = solve(crossprod(smallZ,
                                     bdiag(diag(p*(1-p))))%*%smallZ/num_timepoints)
          
          S = rmvnorm(1,sigma = as.matrix(alphavar) )
          Ss[i,] = S
        }
      }
    }
    
    if (num_clus==num_timepoints) {
      
      mix_quantiles1 = quantile(Ss[,1],probs = c(0.025,0.975))/sqrt(num_timepoints)
      mix_quantiles2 = quantile(Ss[,2],probs = c(0.025,0.975))/sqrt(num_timepoints)
      lcomb_quantiles1 = quantile(SSS[,1],probs = c(0.025,0.975))/sqrt(num_timepoints)
      lcomb_quantiles2 = quantile(SSS[,2],probs = c(0.025,0.975))/sqrt(num_timepoints)
      
      
    } else {
      if (num_clus < num_timepoints) {
        
        
        mix_quantiles1 = c(-1.96*sqrt(true_G[1,1]/(num_clus)), 1.96*sqrt(true_G[1,1]/(num_clus)))
        mix_quantiles2 = c(-1.96*sqrt(true_G[2,2]/(num_clus)), 1.96*sqrt(true_G[2,2]/(num_clus)))
        
        lcomb_quantiles1 = quantile(SSS[,1],probs = c(0.025,0.975))/sqrt(num_timepoints)
        lcomb_quantiles2 = quantile(SSS[,2],probs = c(0.025,0.975))/sqrt(num_timepoints)
        
      } else {
        
        mix_quantiles1 = quantile(Ss[,1],probs = c(0.025,0.975))/sqrt(num_timepoints)
        mix_quantiles2 = quantile(Ss[,2],probs = c(0.025,0.975))/sqrt(num_timepoints)
        lcomb_quantiles1 = mix_quantiles1
        lcomb_quantiles2 = mix_quantiles2
        
      }
    }
    
    ####### higher order corrected intervals for prediction gap - does not depend on relative rates of n and m #########
    
    SSSS = matrix(nrow=10000,ncol=num_Z)
    
    for (i in 1:10000) {
      bs=rmvnorm(num_clus,sigma = true_G)
      b1=bs[1, , drop = FALSE]
      
      eta <- c(smallX %*% c(t(true_beta)) + smallZ %*% t(b1))
      p = exp(eta)/(1+exp(eta))
      
      alphavar = solve(crossprod(smallZ,
                                 bdiag(diag(p*(1-p))))%*%smallZ/num_timepoints)
      
      S = rmvnorm(1, mean = apply(bs, 2, sum)/num_clus , sigma = as.matrix(alphavar)/num_timepoints )
      SSSS[i,] = S
    }
    
    corrected_quantiles1 = quantile(SSSS[,1],probs = c(0.025,0.975))
    corrected_quantiles2 = quantile(SSSS[,2],probs = c(0.025,0.975))
    
    
    tic <- proc.time()
    results <- foreach(i=1:num_sims,.packages = c("MASS","mvtnorm","doParallel","glm2","GGally","Matrix","glmmTMB","dplyr")) %dopar% {
      skip_to_next <- FALSE
      tryCatch(result <- onedatasetsim(i = i), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next) { return(NULL) } else { return(result) }
    }
    toc <- proc.time()
    
    null_indices = sapply(results, function(x) is.null(x))
    results <- results[!null_indices]
    num_sims = num_sims - sum(null_indices)
    
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
    
    #shapiro-wilk p.values
    fixedPQLshapiro1[r,s] = shapiro.test(betadiffs)$p.value
    fixedPQLshapiro2[r,s] = shapiro.test(betadiffs1)$p.value
    randomPQLshapiro1[r,s] = shapiro.test(alphas)$p.value
    randomPQLshapiro2[r,s] = shapiro.test(alphas1)$p.value
    diffsPQLshapiro1[r,s] = shapiro.test(alphadiffs)$p.value
    diffsPQLshapiro2[r,s] = shapiro.test(alphadiffs1)$p.value
    
    ########## Build coverage intervals for fixed and random effects using normal and scale mixture normal quantiles ###########
    
    normcount=0
    normcount1=0
    mixcount=0
    mixcount1=0
    lcombcount = 0
    lcombcount1 = 0 
    
    corrected_count = 0
    corrected_count1 = 0
    
    for (i in 1:num_sims) {
      norm_CI = c(results[[i]]$unifits[[1]]$beta[1] - 1.96*sqrt(true_G[1,1]/num_clus), 
                  results[[i]]$unifits[[1]]$beta[1] + 1.96*sqrt(true_G[1,1]/num_clus))
      if (true_beta[1]<norm_CI[2]&true_beta[1]>norm_CI[1]) {
        normcount=normcount+1
      }
      
      norm_CI1 = c(results[[i]]$unifits[[1]]$beta[2] - 1.96*sqrt(true_G[2,2]/num_clus), 
                   results[[i]]$unifits[[1]]$beta[2] + 1.96*sqrt(true_G[2,2]/num_clus))
      if (true_beta[2]<norm_CI1[2]&true_beta[2]>norm_CI1[1]) {
        normcount1=normcount1+1
      }
      
      mix_PI1 = c(results[[i]]$unifits[[1]]$alpha[1,1] - mix_quantiles1[2] ,results[[i]]$unifits[[1]]$alpha[1,1] - mix_quantiles1[1])
      
      if (results[[i]]$true_alphas[1,1]<mix_PI1[2]&results[[i]]$true_alphas[1,1]>mix_PI1[1]) {
        mixcount=mixcount+1
      }
      
      mix_PI2 = c(results[[i]]$unifits[[1]]$alpha[1,2] - mix_quantiles2[2] ,results[[i]]$unifits[[1]]$alpha[1,2] - mix_quantiles2[1])
      
      if (results[[i]]$true_alphas[1,2]<mix_PI2[2]&results[[i]]$true_alphas[1,2]>mix_PI2[1]) {
        mixcount1=mixcount1+1
      }
      
      lcomb_PI1 = c(results[[i]]$unifits[[1]]$beta[1] + results[[i]]$unifits[[1]]$alpha[1,1] - lcomb_quantiles1[2] ,
                    results[[i]]$unifits[[1]]$beta[1] + results[[i]]$unifits[[1]]$alpha[1,1] - lcomb_quantiles1[1])
      
      if ((true_beta[1] + results[[i]]$true_alphas[1,1])<lcomb_PI1[2]& (true_beta[1] + results[[i]]$true_alphas[1,1]) >lcomb_PI1[1]) {
        lcombcount=lcombcount+1
      }
      
      lcomb_PI2 = c(results[[i]]$unifits[[1]]$beta[2] + results[[i]]$unifits[[1]]$alpha[1,2] - lcomb_quantiles2[2] ,
                    results[[i]]$unifits[[1]]$beta[2] + results[[i]]$unifits[[1]]$alpha[1,2] - lcomb_quantiles2[1])
      
      if ((true_beta[2] + results[[i]]$true_alphas[1,2])<lcomb_PI2[2]& (true_beta[2] + results[[i]]$true_alphas[1,2])>lcomb_PI2[1]) {
        lcombcount1=lcombcount1+1
      }
      
      ############----------Corrected for higher order, prediction gap-----------##################
      
      corrected_PI1 = c(results[[i]]$unifits[[1]]$alpha[1,1] - corrected_quantiles1[2] ,results[[i]]$unifits[[1]]$alpha[1,1] - corrected_quantiles1[1])
      
      if (results[[i]]$true_alphas[1,1]<corrected_PI1[2]&results[[i]]$true_alphas[1,1]>corrected_PI1[1]) {
        corrected_count=corrected_count+1
      }
      
      corrected_PI2 = c(results[[i]]$unifits[[1]]$alpha[1,2] - corrected_quantiles2[2] ,results[[i]]$unifits[[1]]$alpha[1,2] - corrected_quantiles2[1])
      
      if (results[[i]]$true_alphas[1,2]<corrected_PI2[2]&results[[i]]$true_alphas[1,2]>corrected_PI2[1]) {
        corrected_count1=corrected_count1+1
      }
    }
    
    
    fixedcoverage1[r,s] = normcount/num_sims
    fixedcoverage2[r,s] = normcount1/num_sims
    randomcoverage1[r,s] = mixcount/num_sims
    randomcoverage2[r,s] = mixcount1/num_sims
    lcombcoverage1[r,s] = lcombcount/num_sims
    lcombcoverage2[r,s] = lcombcount1/num_sims
    
    corrected_randomcoverage1[r,s] = corrected_count/num_sims
    corrected_randomcoverage2[r,s] = corrected_count1/num_sims
    
    
    ############---------- Frobenius norm of est_G - true_G ------------##############
    
    frobenius = 0
    for (i in 1:num_sims) {
      est_G = crossprod(results[[i]]$unifits[[1]]$alpha)/num_clus
      frobenius = frobenius + sqrt(sum((est_G - true_G)^2))
    }
    frobenius = frobenius/num_sims
    frobenius_norm[r,s] = frobenius
    
    
  }
}


#plots

fixedcoverage1 %>% c %>% cbind(rep(c(25,50,100,200),4)) %>% cbind(rep(c(25,50,100,200),each=4)) %>%
  cbind(rep("fixed_int")) %>% as.data.frame -> fixedcoverage1
colnames(fixedcoverage1) <- c("coverage","num_clus","clus_size","quantity")

fixedcoverage2 %>% c %>% cbind(rep(c(25,50,100,200),4)) %>% cbind(rep(c(25,50,100,200),each=4)) %>%
  cbind(rep("fixed_slope")) %>% as.data.frame -> fixedcoverage2
colnames(fixedcoverage2) <- c("coverage","num_clus","clus_size","quantity")

corrected_randomcoverage1 %>% c %>% cbind(rep(c(25,50,100,200),4)) %>% cbind(rep(c(25,50,100,200),each=4)) %>%
  cbind(rep("rand_int")) %>% as.data.frame -> corrected_randomcoverage1
colnames(corrected_randomcoverage1) <- c("coverage","num_clus","clus_size","quantity")

corrected_randomcoverage2 %>% c %>% cbind(rep(c(25,50,100,200),4)) %>% cbind(rep(c(25,50,100,200),each=4)) %>%
  cbind(rep("rand_slope")) %>% as.data.frame -> corrected_randomcoverage2
colnames(corrected_randomcoverage2) <- c("coverage","num_clus","clus_size","quantity")

lcombcoverage1 %>% c %>% cbind(rep(c(25,50,100,200),4)) %>% cbind(rep(c(25,50,100,200),each=4)) %>%
  cbind(rep("lcomb_int")) %>% as.data.frame -> lcombcoverage1
colnames(lcombcoverage1) <- c("coverage","num_clus","clus_size","quantity")

lcombcoverage2 %>% c %>% cbind(rep(c(25,50,100,200),4)) %>% cbind(rep(c(25,50,100,200),each=4)) %>%
  cbind(rep("lcomb_slope")) %>% as.data.frame -> lcombcoverage2
colnames(lcombcoverage2) <- c("coverage","num_clus","clus_size","quantity")

coverages = rbind(fixedcoverage1,fixedcoverage2,corrected_randomcoverage1,corrected_randomcoverage2,lcombcoverage1,lcombcoverage2)
coverages[,-4] =  apply(coverages[,-4],c(1,2),as.numeric)
coverages$quantity = factor(coverages$quantity,levels = c("fixed_int","fixed_slope","rand_int","rand_slope",'lcomb_int','lcomb_slope'))

quantity.labs <- c("Fixed Intercept", "Fixed Slope", "Random Intercept", "Random Slope", "Fixed + Random Intercept", "Fixed + Random Slope")
names(quantity.labs) <- c("fixed_int", "fixed_slope", "rand_int", "rand_slope", "lcomb_int", "lcomb_slope")

p1 = coverages %>% ggplot(aes(x = clus_size, y = coverage)) +
  geom_line(aes(col = as.factor(num_clus))) +
  facet_wrap(~ quantity, labeller = labeller(quantity = quantity.labs)) +
  theme_bw() +
  labs(x = 'Cluster size' , y = 'Empirical Coverage', col = 'Number of Clusters') + 
  geom_hline(yintercept = 0.95, lty = 2) +
  coord_cartesian(ylim = c(0.8,1)) +
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

diffsPQLshapiro1 %>% c %>% cbind(rep(c(25,50,100,200),4)) %>% cbind(rep(c(25,50,100,200),each=4)) %>%
  cbind(rep("diff_int")) %>% as.data.frame -> diffsPQLshapiro1
colnames(diffsPQLshapiro1) <- c("PQLshapiro","num_clus","clus_size","quantity")

diffsPQLshapiro2 %>% c %>% cbind(rep(c(25,50,100,200),4)) %>% cbind(rep(c(25,50,100,200),each=4)) %>%
  cbind(rep("diff_slope")) %>% as.data.frame -> diffsPQLshapiro2
colnames(diffsPQLshapiro2) <- c("PQLshapiro","num_clus","clus_size","quantity")

PQLshapiros = rbind(fixedPQLshapiro1,fixedPQLshapiro2,randomPQLshapiro1,randomPQLshapiro2,diffsPQLshapiro1,diffsPQLshapiro2)
PQLshapiros[,-4] =  apply(PQLshapiros[,-4],c(1,2),as.numeric)
PQLshapiros$quantity = factor(PQLshapiros$quantity,levels = c("fixed_int","fixed_slope","rand_int","rand_slope",'diff_int','diff_slope'))

quantity.labs <- c("Fixed Intercept", "Fixed Slope", "Random Intercept", "Random Slope", "Predgap Intercept", "Predgap Slope")
names(quantity.labs) <- c("fixed_int", "fixed_slope", "rand_int", "rand_slope", "diff_int", "diff_slope")

p2 = PQLshapiros %>% ggplot(aes(x = clus_size, y = PQLshapiro)) +
  geom_line(aes(col = as.factor(num_clus))) +
  facet_wrap(~ quantity, labeller = labeller(quantity = quantity.labs)) +
  theme_bw() +
  labs(x = 'Cluster size' , y = 'P-value', col = 'Number of Clusters') + 
  geom_hline(yintercept = 0.05, lty = 2) +
  coord_cartesian(ylim = c(0,1)) +
  scale_colour_viridis_d()

p1
ggsave("uncond_bin_coverage_estG.pdf", width = 7, height = 7)
p2
ggsave("uncond_bin_shapiro_estG.pdf", width = 7, height = 7)

frobenius_norm %>% xtable %>% print

save(coverages,PQLshapiros, file = "uncond_bin_estG.Rdata")

