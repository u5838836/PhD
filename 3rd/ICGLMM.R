rm(list = ls())
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
#registerDoMC(cores=28) # this should equal ncpus

setwd("D:/Study new/PhD/R code/3rd")

set.seed(261)
cl <- makeCluster(detectCores()-1)
#cl <- detectCores()-2
registerDoParallel(cl)


clus_test = c(25,50,100,200)
timepoint_test = c(25,50,100,200)

# fixedcoverage1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
# fixedcoverage2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
# randomcoverage1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
# randomcoverage2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
lcombcoverage1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
lcombcoverage2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
# TMB_fixedcoverage1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
# TMB_fixedcoverage2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
TMB_randomcoverage1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
TMB_randomcoverage2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
TMB_lcombcoverage1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
TMB_lcombcoverage2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
# lme_fixedcoverage1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
# lme_fixedcoverage2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
lme_randomcoverage1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
lme_randomcoverage2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
lme_lcombcoverage1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
lme_lcombcoverage2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))


# fixedPQLshapiro1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
# fixedPQLshapiro2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
# randomPQLshapiro1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
# randomPQLshapiro2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
# diffsPQLshapiro1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
# diffsPQLshapiro2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
# 
# fixedTMBshapiro1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
# fixedTMBshapiro2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
# randomTMBshapiro1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
# randomTMBshapiro2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
# diffsTMBshapiro1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
# diffsTMBshapiro2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
# 
# fixedlmeshapiro1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
# fixedlmeshapiro2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
# randomlmeshapiro1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
# randomlmeshapiro2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
# diffslmeshapiro1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
# diffslmeshapiro2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))


our_intlength1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
our_intlength2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
TMB_intlength1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
TMB_intlength2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
lme_intlength1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
lme_intlength2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))

corrected_randomcoverage1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
corrected_randomcoverage2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))

TMB_intlength_var1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
TMB_intlength_var2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
lme_intlength_var1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
lme_intlength_var2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))

fixedsdest_TMB1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
fixedsdest_TMB2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))

predgap_varest_TMB1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
predgap_varest_TMB2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
predgap_varest_lme1 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))
predgap_varest_lme2 = matrix(nrow=length(clus_test),ncol=length(timepoint_test))


########--------Define the sim function for one dataset--------########

onedatasetsim <- function(i) {     
  
  #message("Onto simulated dataset", i)
  ###Generate the observations
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
  # for (j in 1:(num_resp-1)) {
  #   bigZ <- cbind(bigZ,bigX[,(j*(num_X+1)+1):(j*(num_X+1)+num_Z)])
  # }
  set.seed(NULL)
  
  simdat$eta <- c(bigX %*% c(t(true_beta)) + rowSums(bigZ * true_alpha[simdat$ID,]))
  simdat$y <- rpois(n = nrow(simdat), lambda = exp(simdat$eta))
  #rm(bigX, bigZ)

  ##--------------
  ## Do some univariate GLMM fits first via Laplace. These are needed for all the methods anyway
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
    fit_resps[[j]] <- glmmTMB(formula = y ~  x_int + x_slope -1 + (z_int + z_slope - 1|ID), family = "poisson",data = subset(simdat, resp == j), 
                              map = list(theta = factor(rep(NA,length(chol_start)))), start = list(theta = chol_start))
  }
  
  fit_glmer <- glmer(formula = y ~  x_int + x_slope -1 + (z_int + z_slope - 1|ID), family = "poisson",data = simdat)
  
  out <- list(true_alphas = true_alpha,TMBfit = fit_resps[[1]],lmefit = fit_glmer) 
  return(out)
}



########-----------Run the sims--------###########

for (r in 1:length(clus_test)) {
  for (s in 1:length(timepoint_test)) {
    num_resp <- 1
    num_clus <- clus_test[r]
    num_timepoints <- timepoint_test[s]
    num_X <- 1 # Excludes intercepts
    num_Z <- 2
    num_sims = 3

    true_beta = c(-0.1,0.1)
    true_G = 1*diag(2)
    
    true_G_chol=t(chol(cov2cor(true_G)))
    true_G_chol=diag(diag(1/true_G_chol))%*%true_G_chol
    chol_start = log(sqrt(diag(true_G)))
    chol_start = c(chol_start,true_G_chol[lower.tri(true_G_chol)])
    
    #### Get mixture quantiles ####
    
    set.seed(2530)
    bigX <- lapply(1:num_clus, function(i) kronecker(diag(num_resp), cbind(1,rmvnorm(num_timepoints, mean = rep(0,num_X),sigma = 1*diag(num_X)))) ) %>%
      do.call(rbind, .)
    bigZ <- bigX[,1:(num_Z),drop=FALSE]
    
    smallX = bigX[1:num_timepoints,]
    smallZ = smallX
    
    #######higher order corrected intervals#########
    
    SSS = matrix(nrow=10000,ncol=num_Z)
    SSSS = matrix(nrow=10000,ncol=num_Z)
    
    for (i in 1:10000) {
      bs=rmvnorm(num_clus,sigma = true_G)
      b1=bs[1, , drop = FALSE]
      
      eta <- c(smallX %*% c(t(true_beta)) + smallZ %*% t(b1))
      
      alphavar = solve(crossprod(smallZ,
                                 bdiag(diag(exp(eta))))%*%smallZ/num_timepoints)
      
      S = rmvnorm(1, sigma = as.matrix(alphavar)/num_timepoints)
      SS = S + rmvnorm(1, sigma = true_G/num_clus)
      SSS[i,] = S
      SSSS[i,] = SS
    }
    
    corrected_quantiles1 = quantile(SSSS[,1],probs = c(0.025,0.975))
    corrected_quantiles2 = quantile(SSSS[,2],probs = c(0.025,0.975))
    lcomb_quantiles1 = quantile(SSS[,1],probs = c(0.025,0.975))
    lcomb_quantiles2 = quantile(SSS[,2],probs = c(0.025,0.975))
    
    
    tic <- proc.time()
    results <- foreach(i=1:num_sims,.packages = c("MASS","mvtnorm","doParallel","glm2","GGally","Matrix","glmmTMB","dplyr","lme4")) %dopar% {
      skip_to_next <- FALSE
      tryCatch(result <- onedatasetsim(i = i), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next) { return(NULL) } else { return(result) }
    }
    toc <- proc.time()
    
    null_indices = sapply(results, function(x) is.null(x))
    results <- results[!null_indices]
    num_sims = num_sims - sum(null_indices)
    
    
    #save(results,true_beta,true_G,file=paste0(num_clus,num_timepoints,"pois_uncond.Rdata"))
    
    
    ###random effects MSEs
    # univalphaMSE=0
    # alphadiffs = rep(0,num_sims)
    # alphas = rep(0,num_sims)
    # alphadiffs1 = rep(0,num_sims)
    # alphas1 = rep(0,num_sims)
    
    TMBalphaMSE=0
    TMBdiffs = rep(0,num_sims)
    TMBpreds = rep(0,num_sims)
    TMBdiffs1 = rep(0,num_sims)
    TMBpreds1 = rep(0,num_sims)
    
    lmealphaMSE=0
    lmediffs = rep(0,num_sims)
    lmepreds = rep(0,num_sims)
    lmediffs1 = rep(0,num_sims)
    lmepreds1 = rep(0,num_sims)
    
    for (i in 1:num_sims) {
      # univpreds = ranef(results[[i]]$TMBfit)$cond$ID
      # univalphaMSE =  univalphaMSE + apply((univpreds - results[[i]]$true_alphas)^2,2,mean)
      # alphadiffs[i] = (univpreds - results[[i]]$true_alphas)[1,1]
      # alphas[i] = univpreds[1,1]
      # alphadiffs1[i] = (univpreds - results[[i]]$true_alphas)[1,2]
      # alphas1[i] = univpreds[1,2]
      
      TMBalphaMSE =  TMBalphaMSE + apply((as.matrix(ranef(results[[i]]$TMBfit)$cond$ID) - results[[i]]$true_alphas)^2,2,mean)
      TMBdiffs[i] = (as.matrix(ranef(results[[i]]$TMBfit)$cond$ID)  - results[[i]]$true_alphas)[1,1]
      TMBpreds[i] = (as.matrix(ranef(results[[i]]$TMBfit)$cond$ID))[1,1]
      TMBdiffs1[i] = (as.matrix(ranef(results[[i]]$TMBfit)$cond$ID)  - results[[i]]$true_alphas)[1,2]
      TMBpreds1[i] = (as.matrix(ranef(results[[i]]$TMBfit)$cond$ID))[1,2]
      
      lmealphaMSE =  lmealphaMSE + apply((as.matrix(ranef(results[[i]]$lmefit)$ID) - results[[i]]$true_alphas)^2,2,mean)
      lmediffs[i] = (as.matrix(ranef(results[[i]]$lmefit)$ID)  - results[[i]]$true_alphas)[1,1]
      lmepreds[i] = (as.matrix(ranef(results[[i]]$lmefit)$ID))[1,1]
      lmediffs1[i] = (as.matrix(ranef(results[[i]]$lmefit)$ID)  - results[[i]]$true_alphas)[1,2]
      lmepreds1[i] = (as.matrix(ranef(results[[i]]$lmefit)$ID))[1,2]
    }
    
    #univalphaMSE=univalphaMSE/num_sims
    TMBalphaMSE=TMBalphaMSE/num_sims
    lmealphaMSE=lmealphaMSE/num_sims
    
    ###fixed effects MSEs
    
    # univbetaMSE=0
    # betadiffs=rep(0,num_sims)
    # betadiffs1=rep(0,num_sims)
    
    TMBbetaMSE=0
    TMBbetadiffs=rep(0,num_sims)
    TMBbetadiffs1=rep(0,num_sims)
    
    lmebetaMSE=0
    lmebetadiffs=rep(0,num_sims)
    lmebetadiffs1=rep(0,num_sims)
    
    for (i in 1:num_sims) {
      # univbetas = fixef(results[[i]]$TMBfit)$cond
      # univbetaMSE =  univbetaMSE + (univbetas - true_beta)^2
      # betadiffs[i] = (univbetas - true_beta)[1]
      # betadiffs1[i] = (univbetas - true_beta)[2]
      
      TMBbetas = (fixef(results[[i]]$TMBfit))$cond
      TMBbetaMSE =  TMBbetaMSE + (TMBbetas - true_beta)^2
      TMBbetadiffs[i] = (TMBbetas - true_beta)[1]
      TMBbetadiffs1[i] = (TMBbetas - true_beta)[2]
      
      lmebetas = (fixef(results[[i]]$lmefit))
      lmebetaMSE =  lmebetaMSE + (lmebetas - true_beta)^2
      lmebetadiffs[i] = (lmebetas - true_beta)[1]
      lmebetadiffs1[i] = (lmebetas - true_beta)[2]
    }
    
    #univbetaMSE=univbetaMSE/num_sims
    TMBbetaMSE=TMBbetaMSE/num_sims
    lmebetaMSE=lmebetaMSE/num_sims
    
    ###Some plots
    
    # pdf("PQL_knownG_dists.pdf",family="Times",height=7,width=10)
    # par(mfrow=c(2,3))
    # 
    # hist(betadiffs,probability = TRUE) # fixed intercept
    # lines(density(betadiffs),col='red')
    # hist(alphadiffs,probability = TRUE) # random intercept estimate - random intercept
    # lines(density(alphadiffs),col='red')
    # hist(alphas,probability = TRUE,main='Histogram of alphahats') # random intercept estimate
    # lines(density(alphas),col='red')
    # 
    # hist(betadiffs1,probability = TRUE) # fixed slope
    # lines(density(betadiffs1),col='red')
    # hist(alphadiffs1,probability = TRUE) # random slope estimate - random slope
    # lines(density(alphadiffs1),col='red')
    # hist(alphas1,probability = TRUE,main='Histogram of alphahats1') # random slope estimate
    # lines(density(alphas1),col='red')
    # 
    # dev.off()
    # 
    # 
    # pdf("TMB_known_G.pdf",family="Times",height=7,width=10)
    # par(mfrow=c(2,3))
    # 
    # hist(TMBbetadiffs,probability = TRUE) # fixed intercept
    # lines(density(TMBbetadiffs),col='red')
    # hist(TMBdiffs,probability = TRUE) # random intercept estimate - random intercept
    # lines(density(TMBdiffs),col='red')
    # hist(TMBpreds,probability = TRUE) # random intercept estimate
    # lines(density(TMBpreds),col='red')
    # 
    # hist(TMBbetadiffs1,probability = TRUE) # fixed intercept
    # lines(density(TMBbetadiffs1),col='red')
    # hist(TMBdiffs1,probability = TRUE) # random slope estimate - random slope
    # lines(density(TMBdiffs1),col='red')
    # hist(TMBpreds1,probability = TRUE) # random slope estimate
    # lines(density(TMBpreds1),col='red')
    # 
    # dev.off()
    
    # fixedPQLshapiro1[r,s] = shapiro.test(betadiffs)$p.value
    # fixedPQLshapiro2[r,s] = shapiro.test(betadiffs1)$p.value
    # randomPQLshapiro1[r,s] = shapiro.test(alphas)$p.value
    # randomPQLshapiro2[r,s] = shapiro.test(alphas1)$p.value
    # diffsPQLshapiro1[r,s] = shapiro.test(alphadiffs)$p.value
    # diffsPQLshapiro2[r,s] = shapiro.test(alphadiffs1)$p.value
    # 
    # fixedTMBshapiro1[r,s] = shapiro.test(TMBbetadiffs)$p.value
    # fixedTMBshapiro2[r,s] = shapiro.test(TMBbetadiffs1)$p.value
    # randomTMBshapiro1[r,s] = shapiro.test(TMBpreds)$p.value
    # randomTMBshapiro2[r,s] = shapiro.test(TMBpreds1)$p.value
    # diffsTMBshapiro1[r,s] = shapiro.test(TMBdiffs)$p.value
    # diffsTMBshapiro2[r,s] = shapiro.test(TMBdiffs1)$p.value
    # 
    # fixedlmeshapiro1[r,s] = shapiro.test(lmebetadiffs)$p.value
    # fixedlmeshapiro2[r,s] = shapiro.test(lmebetadiffs1)$p.value
    # randomlmeshapiro1[r,s] = shapiro.test(lmepreds)$p.value
    # randomlmeshapiro2[r,s] = shapiro.test(lmepreds1)$p.value
    # diffslmeshapiro1[r,s] = shapiro.test(lmediffs)$p.value
    # diffslmeshapiro2[r,s] = shapiro.test(lmediffs1)$p.value
    
    
    # kurtosis(alphadiffs)
    # kurtosis(alphadiffs1)
    
    
    ###TMB predictor test results
    
    
    # kurtosis(TMBdiffs)
    # kurtosis(TMBdiffs1)
    
    # num_clus*var(betadiffs)
    # num_clus*var(betadiffs1)
    # num_clus*var(TMBbetadiffs)
    # num_clus*var(TMBbetadiffs1)
    # 
    # var(alphas)
    # var(alphas1)
    # var(TMBpreds)
    # var(TMBpreds1)
    
    #####Confidence intervals#####
    
    # normcount=0
    # normcount1=0
    # mixcount=0
    # mixcount1=0
    lcombcount = 0
    lcombcount1 = 0 
    
    # TMBcount=0
    # TMBcount1=0
    TMB_mixcount=0
    TMB_mixcount1=0
    TMB_lcombcount = 0
    TMB_lcombcount1 = 0 
    
    # lmecount=0
    # lmecount1=0
    lme_mixcount=0
    lme_mixcount1=0
    lme_lcombcount = 0
    lme_lcombcount1 = 0 
    
    corrected_count = 0
    corrected_count1 = 0
    
    ourintlength1 = rep(0,num_sims)
    ourintlength2 = rep(0,num_sims)
    TMBintlength1 = rep(0,num_sims)
    TMBintlength2 = rep(0,num_sims)
    lmeintlength1 = rep(0,num_sims)
    lmeintlength2 = rep(0,num_sims)
    
    
    for (i in 1:num_sims) {
      # norm_CI = c(fixef(results[[i]]$TMBfit)$cond[1] - 1.96*sqrt(true_G[1,1]/num_clus), 
      #             fixef(results[[i]]$TMBfit)$cond[1] + 1.96*sqrt(true_G[1,1]/num_clus))
      # if (true_beta[1]<norm_CI[2]&true_beta[1]>norm_CI[1]) {
      #   normcount=normcount+1
      # }
      # 
      # norm_CI1 = c(fixef(results[[i]]$TMBfit)$cond[2] - 1.96*sqrt(true_G[2,2]/num_clus), 
      #              fixef(results[[i]]$TMBfit)$cond[2] + 1.96*sqrt(true_G[2,2]/num_clus))
      # if (true_beta[2]<norm_CI1[2]&true_beta[2]>norm_CI1[1]) {
      #   normcount1=normcount1+1
      # }
      # 
      # 
      # mix_PI1 = c(ranef(results[[i]]$TMBfit)$cond$ID[1,1] - mix_quantiles1[2] ,ranef(results[[i]]$TMBfit)$cond$ID[1,1] - mix_quantiles1[1])
      # 
      # if (results[[i]]$true_alphas[1,1]<mix_PI1[2]&results[[i]]$true_alphas[1,1]>mix_PI1[1]) {
      #   mixcount=mixcount+1
      # }
      # 
      # mix_PI2 = c(ranef(results[[i]]$TMBfit)$cond$ID[1,2] - mix_quantiles2[2] ,ranef(results[[i]]$TMBfit)$cond$ID[1,2] - mix_quantiles2[1])
      # 
      # if (results[[i]]$true_alphas[1,2]<mix_PI2[2]&results[[i]]$true_alphas[1,2]>mix_PI2[1]) {
      #   mixcount1=mixcount1+1
      # }
      
      
      lcomb_PI1 = c(fixef(results[[i]]$TMBfit)$cond[1] + ranef(results[[i]]$TMBfit)$cond$ID[1,1] - lcomb_quantiles1[2] ,
                    fixef(results[[i]]$TMBfit)$cond[1] + ranef(results[[i]]$TMBfit)$cond$ID[1,1] - lcomb_quantiles1[1])
      
      if ((true_beta[1] + results[[i]]$true_alphas[1,1])<lcomb_PI1[2]& (true_beta[1] + results[[i]]$true_alphas[1,1]) >lcomb_PI1[1]) {
        lcombcount=lcombcount+1
      }
      
      lcomb_PI2 = c(fixef(results[[i]]$TMBfit)$cond[2] + ranef(results[[i]]$TMBfit)$cond$ID[1,2] - lcomb_quantiles2[2] ,
                    fixef(results[[i]]$TMBfit)$cond[2] + ranef(results[[i]]$TMBfit)$cond$ID[1,2] - lcomb_quantiles2[1])
      
      if ((true_beta[2] + results[[i]]$true_alphas[1,2])<lcomb_PI2[2]& (true_beta[2] + results[[i]]$true_alphas[1,2])>lcomb_PI2[1]) {
        lcombcount1=lcombcount1+1
      }
      
      ############----------Corrected for higher order-----------##################
      
      corrected_PI1 = c(ranef(results[[i]]$TMBfit)$cond$ID[1,1] - corrected_quantiles1[2] ,ranef(results[[i]]$TMBfit)$cond$ID[1,1] - corrected_quantiles1[1])
      
      if (results[[i]]$true_alphas[1,1]<corrected_PI1[2]&results[[i]]$true_alphas[1,1]>corrected_PI1[1]) {
        corrected_count=corrected_count+1
      }
      
      corrected_PI2 = c(ranef(results[[i]]$TMBfit)$cond$ID[1,2] - corrected_quantiles2[2] ,ranef(results[[i]]$TMBfit)$cond$ID[1,2] - corrected_quantiles2[1])
      
      if (results[[i]]$true_alphas[1,2]<corrected_PI2[2]&results[[i]]$true_alphas[1,2]>corrected_PI2[1]) {
        corrected_count1=corrected_count1+1
      }
    
      
      ############-------------Intervals with naive glmmTMB-----------------##########
      
      
      
      # TMBCI = confint(results[[i]]$TMBfit,estimate=FALSE)[1,]
      # TMBCI1 = confint(results[[i]]$TMBfit,estimate=FALSE)[2,]
      
      # TMBCI = c(fixef(results[[i]]$TMBfit)$cond[1]- 1.96*sqrt(diag(vcov(results[[i]]$TMBfit)$cond))[1],
      #           fixef(results[[i]]$TMBfit)$cond[1]+ 1.96*sqrt(diag(vcov(results[[i]]$TMBfit)$cond))[1])
      # TMBCI1 = c(fixef(results[[i]]$TMBfit)$cond[2]- 1.96*sqrt(diag(vcov(results[[i]]$TMBfit)$cond))[2],
      #            fixef(results[[i]]$TMBfit)$cond[2]+ 1.96*sqrt(diag(vcov(results[[i]]$TMBfit)$cond))[2])
      # 
      # if (true_beta[1]<TMBCI[2]&true_beta[1]>TMBCI[1]) {
      #   TMBcount=TMBcount+1
      # }
      # if (true_beta[2]<TMBCI1[2]&true_beta[2]>TMBCI1[1]) {
      #   TMBcount1=TMBcount1+1
      # }
      
      
      TMBalpha_int_sd = results[[i]]$TMBfit %>% ranef %>% as.data.frame %>% .$condsd %>% .[1]
      TMBalpha_slope_sd = results[[i]]$TMBfit %>% ranef %>% as.data.frame %>% .$condsd %>% .[num_clus+1]
      
      
      TMB_mix_PI1 = c(ranef(results[[i]]$TMBfit)$cond$ID[1,1] - 1.96*TMBalpha_int_sd , ranef(results[[i]]$TMBfit)$cond$ID[1,1] + 1.96*TMBalpha_int_sd)
      
      if (results[[i]]$true_alphas[1,1]<TMB_mix_PI1[2]&results[[i]]$true_alphas[1,1]>TMB_mix_PI1[1]) {
        TMB_mixcount=TMB_mixcount+1
      }
      
      TMB_mix_PI2 = c(ranef(results[[i]]$TMBfit)$cond$ID[1,2] - 1.96*TMBalpha_slope_sd , ranef(results[[i]]$TMBfit)$cond$ID[1,2] + 1.96*TMBalpha_slope_sd)
      
      if (results[[i]]$true_alphas[1,2]<TMB_mix_PI2[2]&results[[i]]$true_alphas[1,2]>TMB_mix_PI2[1]) {
        TMB_mixcount1=TMB_mixcount1+1
      }
      
      
      lcomb1_est = predict(results[[i]]$TMBfit,data.frame(x_int = 1 , x_slope = 0, z_int = 1 , z_slope = 0,  ID = 1) , se.fit = TRUE)
      lcomb2_est = predict(results[[i]]$TMBfit,data.frame(x_int = 0 , x_slope = 1, z_int = 0 , z_slope = 1,  ID = 1) , se.fit = TRUE)
      
      
      TMB_lcomb_PI1 = c(lcomb1_est$fit - 1.96*lcomb1_est$se.fit , lcomb1_est$fit + 1.96*lcomb1_est$se.fit)
      
      if ((true_beta[1] + results[[i]]$true_alphas[1,1])<TMB_lcomb_PI1[2]& (true_beta[1] + results[[i]]$true_alphas[1,1]) >TMB_lcomb_PI1[1]) {
        TMB_lcombcount=TMB_lcombcount+1
      }
      
      TMB_lcomb_PI2 = c(lcomb2_est$fit - 1.96*lcomb2_est$se.fit , lcomb2_est$fit + 1.96*lcomb2_est$se.fit)
      
      
      if ((true_beta[2] + results[[i]]$true_alphas[1,2])<TMB_lcomb_PI2[2]& (true_beta[2] + results[[i]]$true_alphas[1,2])>TMB_lcomb_PI2[1]) {
        TMB_lcombcount1=TMB_lcombcount1+1
      }
      
      
      
      ############-------------Intervals with naive lme4-----------------##########
      
      # lmeCI = c(fixef(results[[i]]$lmefit)[1]- 1.96*sqrt(diag(vcov(results[[i]]$lmefit)))[1],
      #           fixef(results[[i]]$lmefit)[1]+ 1.96*sqrt(diag(vcov(results[[i]]$lmefit)))[1])
      # lmeCI1 = c(fixef(results[[i]]$lmefit)[2]- 1.96*sqrt(diag(vcov(results[[i]]$lmefit)))[2],
      #            fixef(results[[i]]$lmefit)[2]+ 1.96*sqrt(diag(vcov(results[[i]]$lmefit)))[2])
      # 
      # # lmeCI = confint(results[[i]]$lmefit)[1,]
      # # lmeCI1 = confint(results[[i]]$lmefit)[2,]
      # 
      # if (true_beta[1]<lmeCI[2]&true_beta[1]>lmeCI[1]) {
      #   lmecount=lmecount+1
      # }
      # if (true_beta[2]<lmeCI1[2]&true_beta[2]>lmeCI1[1]) {
      #   lmecount1=lmecount1+1
      # }
      
      lmealpha_int_sd = results[[i]]$lmefit %>% ranef %>% as.data.frame %>% .$condsd %>% .[1]
      lmealpha_slope_sd = results[[i]]$lmefit %>% ranef %>% as.data.frame %>% .$condsd %>% .[num_clus+1]
      
      
      lme_mix_PI1 = c(ranef(results[[i]]$lmefit)$ID[1,1] - 1.96*lmealpha_int_sd , ranef(results[[i]]$lmefit)$ID[1,1] + 1.96*lmealpha_int_sd)
      
      if (results[[i]]$true_alphas[1,1]<lme_mix_PI1[2]&results[[i]]$true_alphas[1,1]>lme_mix_PI1[1]) {
        lme_mixcount=lme_mixcount+1
      }
      
      lme_mix_PI2 = c(ranef(results[[i]]$lmefit)$ID[1,2] - 1.96*lmealpha_slope_sd , ranef(results[[i]]$lmefit)$ID[1,2] + 1.96*lmealpha_slope_sd)
      
      if (results[[i]]$true_alphas[1,2]<lme_mix_PI2[2]&results[[i]]$true_alphas[1,2]>lme_mix_PI2[1]) {
        lme_mixcount1=lme_mixcount1+1
      }
      
      lcomb1_est = predict(results[[i]]$lmefit,data.frame(x_int = 1 , x_slope = 0, z_int = 1 , z_slope = 0,  ID = 1))
      lcomb2_est = predict(results[[i]]$lmefit,data.frame(x_int = 0 , x_slope = 1, z_int = 0 , z_slope = 1,  ID = 1))
      
      lme_lcomb_PI1 = c(lcomb1_est - 1.96*lmealpha_int_sd , lcomb1_est + 1.96*lmealpha_int_sd)
      
      if ((true_beta[1] + results[[i]]$true_alphas[1,1])<lme_lcomb_PI1[2]& (true_beta[1] + results[[i]]$true_alphas[1,1]) >lme_lcomb_PI1[1]) {
        lme_lcombcount=lme_lcombcount+1
      }
      
      lme_lcomb_PI2 = c(lcomb2_est - 1.96*lmealpha_slope_sd , lcomb2_est + 1.96*lmealpha_slope_sd)
      
      
      if ((true_beta[2] + results[[i]]$true_alphas[1,2])<lme_lcomb_PI2[2]& (true_beta[2] + results[[i]]$true_alphas[1,2])>lme_lcomb_PI2[1]) {
        lme_lcombcount1=lme_lcombcount1+1
      }
      
      ourintlength1[i] = corrected_PI1[2]-corrected_PI1[1]
      ourintlength2[i] = corrected_PI2[2]-corrected_PI2[1]
      TMBintlength1[i] = TMB_mix_PI1[2]-TMB_mix_PI1[1]
      TMBintlength2[i] = TMB_mix_PI2[2]-TMB_mix_PI2[1]
      lmeintlength1[i] = lme_mix_PI1[2]-lme_mix_PI1[1]
      lmeintlength2[i] = lme_mix_PI2[2]-lme_mix_PI2[1]
      
    }
    
    
    
    # fixedcoverage1[r,s] = normcount/num_sims
    # fixedcoverage2[r,s] = normcount1/num_sims
    # randomcoverage1[r,s] = mixcount/num_sims
    # randomcoverage2[r,s] = mixcount1/num_sims
    lcombcoverage1[r,s] = lcombcount/num_sims
    lcombcoverage2[r,s] = lcombcount1/num_sims
    
    
    
    # TMB_fixedcoverage1[r,s] = TMBcount/num_sims
    # TMB_fixedcoverage2[r,s] = TMBcount1/num_sims
    TMB_randomcoverage1[r,s] = TMB_mixcount/num_sims
    TMB_randomcoverage2[r,s] = TMB_mixcount1/num_sims
    TMB_lcombcoverage1[r,s] = TMB_lcombcount/num_sims
    TMB_lcombcoverage2[r,s] = TMB_lcombcount1/num_sims
    
    # lme_fixedcoverage1[r,s] = lmecount/num_sims
    # lme_fixedcoverage2[r,s] = lmecount1/num_sims
    lme_randomcoverage1[r,s] = lme_mixcount/num_sims
    lme_randomcoverage2[r,s] = lme_mixcount1/num_sims
    lme_lcombcoverage1[r,s] = lme_lcombcount/num_sims
    lme_lcombcoverage2[r,s] = lme_lcombcount1/num_sims
    
    
    corrected_randomcoverage1[r,s] = corrected_count/num_sims
    corrected_randomcoverage2[r,s] = corrected_count1/num_sims
    
    our_intlength1[r,s] = mean(ourintlength1)
    our_intlength2[r,s] = mean(ourintlength2)
    TMB_intlength1[r,s] = mean(TMBintlength1)
    TMB_intlength2[r,s] = mean(TMBintlength2)
    lme_intlength1[r,s] = mean(lmeintlength1)
    lme_intlength2[r,s] = mean(lmeintlength2)
    
    TMB_intlength_var1[r,s] = var(TMBintlength1)
    TMB_intlength_var2[r,s] = var(TMBintlength2)
    lme_intlength_var1[r,s] = var(lmeintlength1)
    lme_intlength_var2[r,s] = var(lmeintlength2)
    
    fixedsdest_TMB1[r,s] = sapply(results, function(x) sqrt(num_clus*diag(vcov(x$TMBfit)$cond))[1]) %>% mean
    fixedsdest_TMB2[r,s] = sapply(results, function(x) sqrt(num_clus*diag(vcov(x$TMBfit)$cond))[2]) %>% mean
    
    predgap_varest_TMB1[r,s] = sapply(results, function(x) x$TMBfit %>% ranef %>% as.data.frame %>% .$condsd %>% .[1]) %>% mean
    predgap_varest_TMB2[r,s] = sapply(results, function(x) x$TMBfit %>% ranef %>% as.data.frame %>% .$condsd %>% .[num_clus+1] ) %>% mean
    predgap_varest_lme1[r,s] = sapply(results, function(x) x$lmefit %>% ranef %>% as.data.frame %>% .$condsd %>% .[1]) %>% mean
    predgap_varest_lme2[r,s] = sapply(results, function(x) x$lmefit %>% ranef %>% as.data.frame %>% .$condsd %>% .[num_clus+1]) %>% mean
    
  }
}

# pois_coverage = rbind(fixedcoverage1,fixedcoverage2,randomcoverage1,randomcoverage2,lcombcoverage1,lcombcoverage2)
# save(pois_coverage, file="pois_coverage.Rdata")
# print(xtable(pois_coverage,digits=3), include.rownames=FALSE)

pois_coverage_TMB = rbind(TMB_randomcoverage1,TMB_randomcoverage2,TMB_lcombcoverage1,TMB_lcombcoverage2) #TMB_fixedcoverage1,TMB_fixedcoverage2,
#save(pois_coverage_TMB, file="pois_coverage_TMB.Rdata")
print(xtable(pois_coverage_TMB,digits=3), include.rownames=FALSE)

pois_coverage_lme = rbind(lme_randomcoverage1,lme_randomcoverage2,lme_lcombcoverage1,lme_lcombcoverage2) #lme_fixedcoverage1,lme_fixedcoverage2,
#save(pois_coverage_lme, file="pois_coverage_lme.Rdata")
print(xtable(pois_coverage_lme,digits=3), include.rownames=FALSE)

pois_corrected_coverage = rbind(corrected_randomcoverage1,corrected_randomcoverage2,lcombcoverage1,lcombcoverage2)
#save(pois_corrected_coverage, file="pois_coverage_corrected.Rdata")
print(xtable(pois_corrected_coverage,digits=3), include.rownames=FALSE)

# pois_shapiro = rbind(fixedPQLshapiro1,fixedPQLshapiro2,randomPQLshapiro1,randomPQLshapiro2,diffsPQLshapiro1,diffsPQLshapiro2,
#                      fixedTMBshapiro1,fixedTMBshapiro2,randomTMBshapiro1,randomTMBshapiro2,diffsTMBshapiro1,diffsTMBshapiro2,
#                      fixedlmeshapiro1,fixedlmeshapiro2,randomlmeshapiro1,randomlmeshapiro2,diffslmeshapiro1,diffslmeshapiro2)
# save(pois_shapiro, file="pois_shapiro.Rdata")
# print(xtable(pois_shapiro,digits=3), include.rownames=FALSE)


print(xtable(our_intlength1,digits=3), include.rownames=FALSE)
print(xtable(our_intlength2,digits=3), include.rownames=FALSE)
print(xtable(TMB_intlength1,digits=3), include.rownames=FALSE)
print(xtable(TMB_intlength2,digits=3), include.rownames=FALSE)
print(xtable(lme_intlength1,digits=3), include.rownames=FALSE)
print(xtable(lme_intlength2,digits=3), include.rownames=FALSE)
print(xtable(rbind(TMB_intlength_var1,TMB_intlength_var2),digits=3), include.rownames=FALSE)
print(xtable(rbind(lme_intlength_var1,lme_intlength_var2),digits=3), include.rownames=FALSE)

fixedsdest_TMB1 %>% xtable(digits=3) %>% print(include.rownames=FALSE)
fixedsdest_TMB2 %>% xtable(digits=3) %>% print(include.rownames=FALSE)

predgap_varest_TMB1 %>% xtable(digits=3) %>% print(include.rownames=FALSE)
predgap_varest_TMB2 %>% xtable(digits=3) %>% print(include.rownames=FALSE)
predgap_varest_lme1 %>% xtable(digits=3) %>% print(include.rownames=FALSE)
predgap_varest_lme2 %>% xtable(digits=3) %>% print(include.rownames=FALSE)

#normalised fixed effects standard error TMB all n,m. Estimated prediction gap variance TMB vs lme4 n,m. 

stopCluster(cl)


