rm(list = ls())
library(glmmTMB)
library(dplyr)
library(foreach)
library(glm2)

num_clus <- 200
num_timepoints <- 100
num_sims = 3
fixed_covariate = rnorm(num_timepoints*num_clus) %>% matrix(nrow = num_clus,ncol = num_timepoints) %>% c

true_beta = 0.1
true_G = 0.1

chol_start = log(sqrt(true_G)) #fix working G for glmmTMB

PQL_sum_pred_error = rep(0,num_sims)
LA_sum_pred_error = rep(0,num_sims)
PQL_vec_pred_error = rep(0,num_sims)
LA_vec_pred_error = rep(0,num_sims)
PQL_beta_error = rep(0,num_sims)
LA_beta_error = rep(0,num_sims)


#Define PQL
univariate_glmmPQL <- function(starting_fit, max_iter = 10000) {
  
  #set starting values for parameters
  cw_beta <- (fixef(starting_fit))$cond
  cw_G <- true_G
  cw_alpha_mat <- as.matrix(ranef(starting_fit)$cond$ID)
  
  ## Optimization parameters - error tolerance
  err <- 100
  counter <- 0
  while(err > 1e-7 & counter < max_iter) {  
    
    cw_Ginv <- solve(cw_G)
    ## Maximise wrt alpha -- on a per cluster basis
    update_alpha <- function(j) {
      
      innerlogL_alpha <- function(alpha) {
        cw_eta <- fixed_covariate*cw_beta + alpha          
        likcontrib <- -sum(dpois(simdat$y[simdat$ID==j], lambda = exp(cw_eta), log = TRUE)) + 0.5*crossprod(alpha, cw_Ginv) %*% alpha
        return(as.vector(likcontrib))
      }
      
      innergrad_alpha <- function(alpha) {
        cw_eta <- fixed_covariate*cw_beta + alpha          
        out <- -sum(simdat$y[simdat$ID==j] - exp(cw_eta)) + cw_Ginv %*% alpha 
        return(as.vector(out))
      }
      
      do_update <- optim(cw_alpha_mat[j,], fn = innerlogL_alpha, gr = innergrad_alpha, control = list(trace = 0, maxit = max_iter), method = "BFGS")
      return(do_update$par)
    }     
    new_alpha_mat <- foreach(j = 1:num_clus, .combine = "rbind") %do% update_alpha(j = j)
    
    ## Maximise wrt beta
    make_offset <- rep(as.vector(new_alpha_mat),each=num_timepoints)
    update_beta <- glm2(simdat$y ~ -1 + fixed_covariate, offset = make_offset, family = "poisson")
    new_beta <- update_beta$coefficients
    rm(update_beta)
    
    ## G is known
    new_G = true_G
    
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



for (j in 1:num_sims) {
  # Generate the random effects
  true_alpha <- rnorm(num_clus, sd = sqrt(true_G)) %>% matrix(nrow=num_clus)
  
  # Generate the observations
  simdat <- data.frame(ID = factor(rep(1:num_clus, each = num_timepoints)), 
                       time = factor(rep(1:num_timepoints, num_clus)))
  
  simdat$eta <- c(fixed_covariate*true_beta + true_alpha[simdat$ID,])
  simdat$y <- rpois(n = nrow(simdat), lambda = exp(simdat$eta))
  
  ##--------------
  ## GLMM fits via Laplace. 
  ##--------------
  
  fit_resp <- glmmTMB(formula = y ~ -1 + fixed_covariate + (1|ID), family = "poisson",data=simdat, 
                      map = list(theta = factor(rep(NA,length(chol_start)))), start = list(theta = chol_start))
  
  sum(fit_resp$fit$parfull[-1])/true_G #this should be equal to the first derivative of the logdet term since X_i=Z_i
  
  dummy_Z = matrix(0,nrow = num_clus*num_timepoints,ncol = num_clus)
  for (i in 1:num_clus) {
    dummy_Z[((i-1)*num_timepoints+1):(i*num_timepoints),i] = 1
  }
  mu_est = exp(fixed_covariate*fit_resp$fit$parfull[1] + dummy_Z%*%fit_resp$fit$parfull[-1]) #estimated means for each obv
  
  sum(simdat$y - mu_est) - sum(fit_resp$fit$parfull[-1])/true_G #this is correctly very close to 0
  
  #get the first (total!) derivative of the logdet term
  tr_deriv = 0
  for (i in 1:num_clus) {
    tr_deriv = tr_deriv + num_timepoints*exp(fit_resp$fit$parfull[1] + 
                                               fit_resp$fit$parfull[1+i])/((num_timepoints*exp(fit_resp$fit$parfull[1] 
                                                                                               + fit_resp$fit$parfull[1+i]) + 1/true_G)^2)
  }
  sum(simdat$y - mu_est) - 0.5*tr_deriv/true_G
  
  #Fit PQL
  PQL_fit = univariate_glmmPQL(starting_fit = fit_resp)
  
  
  LA_sum_pred_error[j] = (sum(fit_resp$fit$parfull[-1]) - sum(true_alpha))^2 #Square prediction error of Laplace for sum of b_i
  PQL_sum_pred_error[j] = sum(PQL_fit$alpha - true_alpha)^2 #Square prediction error of PQL for sum of b_i
  LA_vec_pred_error[j] = sum((fit_resp$fit$parfull[-1] - true_alpha)^2)
  PQL_vec_pred_error[j] = sum((PQL_fit$alpha - true_alpha)^2)
  LA_beta_error[j] = (fit_resp$fit$parfull[1] - true_beta)^2
  PQL_beta_error[j] = (PQL_fit$beta - true_beta)^2
}

num_clus*true_G #Mean square prediction error of PQL for sum of b_i
mean(LA_sum_pred_error)
mean(PQL_sum_pred_error)
mean(LA_vec_pred_error)
mean(PQL_vec_pred_error)
mean(LA_beta_error)
mean(PQL_beta_error)

#Single fit

as.data.frame(ranef(fit_resp))
estb1 = as.data.frame(ranef(fit_resp))$condval[1]
estb1_sd = as.data.frame(ranef(fit_resp))$condsd[1]
int_est = fit_resp$fit$par
sqrt(1/(num_timepoints*(exp(int_est + estb1) + 1/(true_G^2*num_timepoints))) +
       (1/num_clus)*true_G*((1 - 1/(1 + num_timepoints*true_G*exp(int_est + estb1)))^2))
estb1_sd
sqrt(true_G/num_clus)

