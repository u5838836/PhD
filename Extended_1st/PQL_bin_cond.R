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
library(tidyr)
registerDoMC(cores=28)
set.seed(259)

#The number of clusters and cluster sizes we wish to explore for the overall sim
# clus_test = c(10,11,12,13,14) 
# timepoint_test = c(10,11,12,13,14) 
clus_test = c(25,50,100,200,400)
timepoint_test = c(25,50,100,200,400)

#create storage for quantities of interest

all_estimates = NULL

for (i in 1:5) {
  assign(paste0("fixedcoverage",i),matrix(nrow=length(clus_test),ncol=length(timepoint_test)))
  assign(paste0("randomcoverage",i),matrix(nrow=length(clus_test),ncol=length(timepoint_test)))
  assign(paste0("fixedPQLshapiro",i),matrix(nrow=length(clus_test),ncol=length(timepoint_test)))
  assign(paste0("randomPQLshapiro",i),matrix(nrow=length(clus_test),ncol=length(timepoint_test)))
}

diff_norm = matrix(nrow=length(clus_test),ncol=length(timepoint_test))


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
  bigX <- cbind(1,mvrnorm(n = num_clus*num_timepoints,mu = c(0,0),Sigma = matrix(c(1,0.5,0.5,1),nrow=2)),
                rnorm(num_clus*num_timepoints),rbinom(n = num_clus*num_timepoints,size = 1,prob = 0.5))
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
    while(err > 1e-5 & counter < max_iter) {  
      
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
      #new_G = crossprod(new_alpha_mat)/num_clus
      
      ## Maximise wrt G using fixed point iterative update
      # bigW <- as.vector(binomial()$var(binomial()$linkinv(subX %*% new_beta + make_offset)))       
      # cw_G0 <- cw_G
      # err2 <- 1 
      # while(err2 > 0.01) {
      #   new_G <- sapply(1:num_clus, function(j) {
      #     selindices <- (num_timepoints*(j-1)+1):(num_timepoints*j)
      #     solve(crossprod(as.matrix(subZ[selindices,])*sqrt(bigW[selindices])) + solve(cw_G0)) + tcrossprod(new_alpha_mat[j,]) 
      #   })
      #   new_G <- matrix(rowMeans(new_G), nrow = num_Z)
      #   err2 <- norm(new_G - cw_G0, type = "F")
      #   cw_G0 <- new_G
      # }
      
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
    
    num_X <- 4 # Number of fixed effects, excludes intercept
    num_Z <- 5
    num_sims = 1000
    
    #Set true values
    true_beta = c(-0.1,0.1,-0.1,0.1,0.1)
    true_G = 1*diag(5)
    working_G = 1*diag(5)
    
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
    
    for (i in 1:5) {
      assign(paste0("alphadiffs",i),rep(0,num_sims))
      assign(paste0("alphas",i),rep(0,num_sims))
      assign(paste0("betadiffs",i),rep(0,num_sims))
    }
    
    
    for (i in 1:num_sims) {
      univpreds = results[[i]]$unifits[[1]]$alpha
      univdiffs = univpreds - results[[i]]$true_alphas
      alphadiffs1[i] = univdiffs[1,1]
      alphadiffs2[i] = univdiffs[1,2]
      alphadiffs3[i] = univdiffs[1,3]
      alphadiffs4[i] = univdiffs[1,4]
      alphadiffs5[i] = univdiffs[1,5]
      alphas1[i] = univpreds[1,1]
      alphas2[i] = univpreds[1,2]
      alphas3[i] = univpreds[1,3]
      alphas4[i] = univpreds[1,4]
      alphas5[i] = univpreds[1,5]
    }
    
    
    ###extract fixed effects estimates
    
    for (i in 1:num_sims) {
      univbetas = results[[i]]$unifits[[1]]$beta
      univbetadiffs  = univbetas - true_beta
      betadiffs1[i] = univbetadiffs[1]
      betadiffs2[i] = univbetadiffs[2]
      betadiffs3[i] = univbetadiffs[3]
      betadiffs4[i] = univbetadiffs[4]
      betadiffs5[i] = univbetadiffs[5]
    }
    
    #fixed - shapiro-wilk p-values
    fixedPQLshapiro1[r,s] = shapiro.test(betadiffs1)$p.value
    fixedPQLshapiro2[r,s] = shapiro.test(betadiffs2)$p.value
    fixedPQLshapiro3[r,s] = shapiro.test(betadiffs3)$p.value
    fixedPQLshapiro4[r,s] = shapiro.test(betadiffs4)$p.value
    fixedPQLshapiro5[r,s] = shapiro.test(betadiffs5)$p.value

    #random - shapiro-wilk p-values
    randomPQLshapiro1[r,s] = shapiro.test(alphadiffs1)$p.value
    randomPQLshapiro2[r,s] = shapiro.test(alphadiffs2)$p.value
    randomPQLshapiro3[r,s] = shapiro.test(alphadiffs3)$p.value
    randomPQLshapiro4[r,s] = shapiro.test(alphadiffs4)$p.value
    randomPQLshapiro5[r,s] = shapiro.test(alphadiffs5)$p.value
    
    # Generate same covariates
    set.seed(253)
    bigX <- cbind(1,mvrnorm(n = num_clus*num_timepoints,mu = c(0,0),Sigma = matrix(c(1,0.5,0.5,1),nrow=2)),
                  rnorm(num_clus*num_timepoints),rbinom(n = num_clus*num_timepoints,size = 1,prob = 0.5))
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
    
    fixedcount1=0
    fixedcount2=0
    fixedcount3=0
    fixedcount4=0
    fixedcount5=0
    randcount1=0
    randcount2=0
    randcount3=0
    randcount4=0
    randcount5=0
    
    for (i in 1:num_sims) {
      
      
      for (j in 1:5) {
        
        fixed_CI = c(results[[i]]$unifits[[1]]$beta[j] - 1.96*sqrt(betavar[j,j]/(num_timepoints*num_clus)), 
                     results[[i]]$unifits[[1]]$beta[j] + 1.96*sqrt(betavar[j,j]/(num_timepoints*num_clus)))
        if (true_beta[j]<fixed_CI[2]&true_beta[j]>fixed_CI[1]) {
          assign(paste0("fixedcount",j),get(paste0("fixedcount",j))+1)
        }
        
        
        rand_CI = c(results[[i]]$unifits[[1]]$alpha[1,j] - 1.96*sqrt(alphavar[j,j]/num_timepoints), 
                    results[[i]]$unifits[[1]]$alpha[1,j] + 1.96*sqrt(alphavar[j,j]/num_timepoints))
        if (true_alpha[1,j]<rand_CI[2]&true_alpha[1,j]>rand_CI[1]) {
          assign(paste0("randcount",j),get(paste0("randcount",j))+1)
        }

      }
        

      }
      
    
    PQL_LA_diff = rep(0,num_sims)
    for (i in 1:num_sims) {
      PQL_LA_diff[i] = sqrt(sum((ranef(results[[i]]$lme4fit)$ID - results[[i]]$unifits[[1]]$alpha)^2) + sum((results[[i]]$lme4fit$fit$par - results[[i]]$unifits[[1]]$beta)^2))/num_clus
    }
    
    fixedcoverage1[r,s] = fixedcount1/num_sims
    fixedcoverage2[r,s] = fixedcount2/num_sims
    fixedcoverage3[r,s] = fixedcount3/num_sims
    fixedcoverage4[r,s] = fixedcount4/num_sims
    fixedcoverage5[r,s] = fixedcount5/num_sims
    randomcoverage1[r,s] = randcount1/num_sims
    randomcoverage2[r,s] = randcount2/num_sims
    randomcoverage3[r,s] = randcount3/num_sims
    randomcoverage4[r,s] = randcount4/num_sims
    randomcoverage5[r,s] = randcount5/num_sims

    diff_norm[r,s] = mean(PQL_LA_diff)
    cbind(betadiffs1,betadiffs2,betadiffs3,betadiffs4,betadiffs5,
          alphadiffs1,alphadiffs2,alphadiffs3,alphadiffs4,alphadiffs5) %>%
      as.data.frame %>% mutate(numclus = num_clus, clussize = num_timepoints) %>% rbind(all_estimates) -> all_estimates
  }
}

#plots


for (j in 1:5) {
  assign(paste0("fixedcoverage",j),get(paste0("fixedcoverage",j)) %>% c %>% cbind(rep(c(25,50,100,200,400),5)) %>% cbind(rep(c(25,50,100,200,400),each=5)) %>%
            cbind(rep(paste0("fixed",j))) %>% as.data.frame)

  assign(paste0("randomcoverage",j),get(paste0("randomcoverage",j)) %>% c %>% cbind(rep(c(25,50,100,200,400),5)) %>% cbind(rep(c(25,50,100,200,400),each=5)) %>%
               cbind(rep(paste0("random",j))) %>% as.data.frame )


  assign(paste0("fixedPQLshapiro",j),get(paste0("fixedPQLshapiro",j)) %>% c %>% cbind(rep(c(25,50,100,200,400),5)) %>% cbind(rep(c(25,50,100,200,400),each=5)) %>%
           cbind(rep(paste0("fixed",j))) %>% as.data.frame)


  assign(paste0("randomPQLshapiro",j),get(paste0("randomPQLshapiro",j)) %>% c %>% cbind(rep(c(25,50,100,200,400),5)) %>% cbind(rep(c(25,50,100,200,400),each=5)) %>%
           cbind(rep(paste0("random",j))) %>% as.data.frame)

}

colnames(fixedcoverage1) <- c("coverage","num_clus","clus_size","quantity")
colnames(randomcoverage1) <- c("coverage","num_clus","clus_size","quantity")
colnames(fixedPQLshapiro1) <- c("PQLshapiro","num_clus","clus_size","quantity")
colnames(randomPQLshapiro1) <- c("PQLshapiro","num_clus","clus_size","quantity")
colnames(fixedcoverage2) <- c("coverage","num_clus","clus_size","quantity")
colnames(randomcoverage2) <- c("coverage","num_clus","clus_size","quantity")
colnames(fixedPQLshapiro2) <- c("PQLshapiro","num_clus","clus_size","quantity")
colnames(randomPQLshapiro2) <- c("PQLshapiro","num_clus","clus_size","quantity")
colnames(fixedcoverage3) <- c("coverage","num_clus","clus_size","quantity")
colnames(randomcoverage3) <- c("coverage","num_clus","clus_size","quantity")
colnames(fixedPQLshapiro3) <- c("PQLshapiro","num_clus","clus_size","quantity")
colnames(randomPQLshapiro3) <- c("PQLshapiro","num_clus","clus_size","quantity")
colnames(fixedcoverage4) <- c("coverage","num_clus","clus_size","quantity")
colnames(randomcoverage4) <- c("coverage","num_clus","clus_size","quantity")
colnames(fixedPQLshapiro4) <- c("PQLshapiro","num_clus","clus_size","quantity")
colnames(randomPQLshapiro4) <- c("PQLshapiro","num_clus","clus_size","quantity")
colnames(fixedcoverage5) <- c("coverage","num_clus","clus_size","quantity")
colnames(randomcoverage5) <- c("coverage","num_clus","clus_size","quantity")
colnames(fixedPQLshapiro5) <- c("PQLshapiro","num_clus","clus_size","quantity")
colnames(randomPQLshapiro5) <- c("PQLshapiro","num_clus","clus_size","quantity")


coverages = rbind(fixedcoverage1,fixedcoverage2,fixedcoverage3,fixedcoverage4,fixedcoverage5,
                  randomcoverage1,randomcoverage2,randomcoverage3,randomcoverage4,randomcoverage5)
coverages[,-4] =  apply(coverages[,-4],c(1,2),as.numeric)

PQLshapiros = rbind(fixedPQLshapiro1,fixedPQLshapiro2,fixedPQLshapiro3,fixedPQLshapiro4,fixedPQLshapiro5,
                    randomPQLshapiro1,randomPQLshapiro2,randomPQLshapiro3,randomPQLshapiro4,randomPQLshapiro5)
PQLshapiros[,-4] =  apply(PQLshapiros[,-4],c(1,2),as.numeric)

quantity.labs <- c("Fixed1", "Fixed2","Fixed3","Fixed4","Fixed5",
                   "Rand1","Rand2","Rand3","Rand4","Rand5")
names(quantity.labs) <- c("fixed1", "fixed2","fixed3","fixed4","fixed5",
                          "random1","random2","random3","random4","random5")

p1 = coverages %>% ggplot(aes(x = clus_size, y = coverage)) +
  geom_line(aes(col = as.factor(num_clus))) +
  facet_wrap(~ quantity, labeller = labeller(quantity = quantity.labs), scales = 'free') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Bernoulli Responses",x = 'Cluster size' , y = 'Empirical Coverage', col = 'Number of Clusters') + 
  geom_hline(yintercept = 0.95, lty = 2) +
  #coord_cartesian(ylim = c(0.8,1)) +
  scale_colour_viridis_d() +
  theme(text = element_text(size=15), axis.text.x = element_text(angle = 90,vjust = 0.5))

p1
ggsave("cond_bin_coverage.pdf", width = 7, height = 7)


p2 = PQLshapiros %>% ggplot(aes(x = as.factor(clus_size), y = as.factor(num_clus))) +
  geom_raster(aes(fill = PQLshapiro)) +
  facet_wrap(~ quantity, labeller = labeller(quantity = quantity.labs)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Bernoulli Responses",x = 'Cluster size' , y = 'Number of Clusters', fill = 'P-value') + 
  scale_fill_distiller(palette = "Greens")+
  theme(text = element_text(size=15), axis.text.x = element_text(angle = 90,vjust = 0.5))

p2
ggsave("cond_bin_shapiro.pdf", width = 7, height = 7)


p3 = all_estimates %>% 
  pivot_longer(betadiffs1:alphadiffs5,names_to = 'quantity',values_to = 'estimate') %>%
  filter(quantity == "betadiffs3") %>%
  ggplot(aes(x = estimate)) +
  geom_histogram(aes(y = ..density..)) +
  geom_density(aes(y = ..density..), color = 'red') +
  facet_grid(as.factor(numclus) ~ as.factor(clussize)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Bernoulli Responses",x = '', y = '') + 
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = 0) +
  geom_vline(xintercept = 0,lty = 2, color = 'blue')

p3
ggsave("cond_bin_fixedhist.pdf", width = 7, height = 7)


p4 = all_estimates %>% 
  pivot_longer(betadiffs1:alphadiffs5,names_to = 'quantity',values_to = 'estimate') %>%
  filter(quantity == "alphadiffs3") %>%
  ggplot(aes(x = estimate)) +
  geom_histogram(aes(y = ..density..)) +
  geom_density(aes(y = ..density..), color = 'red') +
  facet_grid(as.factor(numclus) ~ as.factor(clussize)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Bernoulli Responses",x = '', y = '') + 
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = 0) +
  geom_vline(xintercept = 0,lty = 2, color = 'blue')

p4
ggsave("cond_bin_randhist.pdf", width = 7, height = 7)



diff_norm %>% xtable %>% print

save(diff_norm, file = "cond_bin_PQLvsLA.Rdata")

save(coverages,PQLshapiros, file = "cond_bin.Rdata")



