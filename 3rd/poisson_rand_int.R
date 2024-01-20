rm(list = ls())
library(glmmTMB)
library(e1071) 
library(dplyr)
library(xtable)
library(mvtnorm)
library(MASS)
library(Matrix)
library(parallel)
library(foreach)
library(doParallel)
library(tidyverse)
library(GGally)
library(TMB)
library(glm2)
library(lme4)
setwd("D:/Study new/PhD/R code/3rd/Figures")

cl <- makeCluster(detectCores()-1)
#cl <- detectCores()-2
registerDoParallel(cl)

num_resp <- 1
num_clus <- 2
num_timepoints <- 200
num_X <- 0 # Excludes intercepts
num_Z <- 1
num_sims = 1000

true_beta = 0
true_G = 1

chol_start = log(sqrt(true_G))

onedatasetsim <- function(i) {     
  
  message("Onto simulated dataset", i)
  ###Generate the observations
  ##--------------
  ## Simulate data
  ##--------------
  
  # Generate the random effects
  true_alpha <- rnorm(num_clus, sd = sqrt(true_G)) %>% matrix(nrow=num_clus)
  
  # Generate the observations
  simdat <- data.frame(ID = factor(rep(1:num_clus, each = num_resp*num_timepoints)), 
                       time = factor(rep(rep(1:num_timepoints, num_resp), num_clus)),
                       resp = factor(rep(rep(1:num_resp, each = num_timepoints), num_clus)))
  
  # Generate covariate -- assuming all responses have the same set of covariates, which is usually the case
  bigZ <- rep(1,num_clus*num_timepoints)
  
  simdat$eta <- c(bigZ * true_alpha[simdat$ID,])
  simdat$y <- rpois(n = nrow(simdat), lambda = exp(simdat$eta))
  #rm(bigX, bigZ)
  
  
  ##--------------
  ## Do some univariate GLMM fits first via Laplace. These are needed for all the methods anyway
  ##--------------
  bigZ <- as.matrix(bigZ)
  smallZ <- bigZ[,1:num_Z,drop=FALSE]
  smallZ <- smallZ[rowSums(smallZ) != 0,] 
  
  fit_resps=list()
  for (j in 1:num_resp) {
    fit_resps[[j]] <- glmmTMB(formula = y ~ - 1 + (smallZ-1|ID), family = "poisson",data = subset(simdat, resp == j), 
                              map = list(theta = factor(rep(NA,length(chol_start)))), start = list(theta = chol_start))
  }
  
  fit_glmer <- glmer(y ~ -1 + (smallZ - 1|ID), family = poisson, data=simdat)
  
  out <- list(true_alphas = true_alpha,TMBfit = fit_resps[[1]], lme4fit = fit_glmer) 
  return(out)
}

package_string = c("MASS","mvtnorm","doParallel","glm2","GGally","Matrix","glmmTMB","dplyr",'lme4')

tic <- proc.time()
set.seed(261)
results <- foreach(i=1:num_sims,.packages = package_string) %dopar% onedatasetsim(i = i)
toc <- proc.time()

stopCluster(cl)

#save(results,true_G,file="poisson_rand_int.Rdata")


###random effects
# alphadiffs = rep(0,num_sims)
# alphas = rep(0,num_sims)

TMBdiffs = rep(0,num_sims)
TMBpreds = rep(0,num_sims)

lme4diffs = rep(0,num_sims)
lme4preds = rep(0,num_sims)

for (i in 1:num_sims) {
  # univpreds = results[[i]]$unifits[[1]]$alpha
  # alphadiffs[i] = (univpreds - results[[i]]$true_alphas)[1,1]
  # alphas[i] = univpreds[1,1]
  
  TMBdiffs[i] = (as.matrix(ranef(results[[i]]$TMBfit)$cond$ID)  - results[[i]]$true_alphas)[1,1]
  TMBpreds[i] = (as.matrix(ranef(results[[i]]$TMBfit)$cond$ID))[1,1]
  
  lme4diffs[i] = (as.matrix(ranef(results[[i]]$lme4fit)$ID)  - results[[i]]$true_alphas)[1,1]
  lme4preds[i] = (as.matrix(ranef(results[[i]]$lme4fit)$ID))[1,1]
}

# ###Some plots
# 
# pdf("pois_rand_int.pdf",family="Times",height=7,width=10)
# par(mfrow=c(2,2))
# 
# hist(alphadiffs,probability = TRUE,main='PQL',xlab = bquote(hat(alpha)[1] - dot(alpha)[1])) # random intercept estimate - random intercept
# lines(density(alphadiffs),col='red')
# hist(alphas,probability = TRUE,main='PQL',xlab = bquote(hat(alpha)[1])) # random intercept estimate
# lines(density(alphas),col='red')
# 
# hist(TMBdiffs,probability = TRUE,main='TMB',xlab = bquote(hat(alpha)[1] - dot(alpha)[1])) # random intercept estimate - random intercept
# lines(density(TMBdiffs),col='red')
# hist(TMBpreds,probability = TRUE,main='TMB',xlab = bquote(hat(alpha)[1])) # random intercept estimate
# lines(density(TMBpreds),col='red')
# 
# hist(lme4diffs,probability = TRUE,main='lme4',xlab = bquote(hat(alpha)[1] - dot(alpha)[1])) # random intercept estimate - random intercept
# lines(density(lme4diffs),col='red')
# hist(lme4preds,probability = TRUE,main='lme4',xlab = bquote(hat(alpha)[1])) # random intercept estimate
# lines(density(lme4preds),col='red')
# 
# dev.off()


# shapiro.test(alphadiffs)
# 
# shapiro.test(alphas)


###TMB predictor test results

shapiro.test(TMBdiffs)

shapiro.test(TMBpreds)

###lme4 predictor test results

shapiro.test(lme4diffs)

shapiro.test(lme4preds)


#num_timepoints*var(alphadiffs)
num_timepoints*var(TMBdiffs)
num_timepoints*var(lme4diffs)
exp(true_G/2)

#### Get mixture quantiles

Ss = rep(0,10000)

for (i in 1:10000) {
  b1=rnorm(1,sd = sqrt(true_G))
  S = rnorm(1,sd=sqrt(exp(-b1)))
  Ss[i] = S
}

mix_quantiles = quantile(Ss,probs = c(0.025,0.975))/sqrt(num_timepoints)

#norm_quantiles = 1.96*exp(true_G/8)/sqrt(num_timepoints)

#####Prediction intervals#####

normcount=0
mixcount=0
TMB_mixcount=0
lme4_mixcount=0
naive_count=0

mix_intlength = rep(0,num_sims)
TMB_intlength = rep(0,num_sims)
lme4_intlength = rep(0,num_sims)
norm_intlength = rep(0,num_sims)
naive_intlength = rep(0,num_sims)

for (i in 1:num_sims) {
  naive_PI = c(ranef(results[[i]]$TMBfit)$cond$ID[1,1] - 1.96*exp(true_G/4)/sqrt(num_timepoints),
              ranef(results[[i]]$TMBfit)$cond$ID[1,1] + 1.96*exp(true_G/4)/sqrt(num_timepoints))
  if (results[[i]]$true_alphas[1,1]<naive_PI[2]&results[[i]]$true_alphas[1,1]>naive_PI[1]) {
    naive_count=naive_count+1
  }
  
  naive_intlength[i] = naive_PI[2] - naive_PI[1]
  
  # norm_PI = c(ranef(results[[i]]$TMBfit)$cond$ID[1,1] - 1.96*sqrt(exp(-ranef(results[[i]]$TMBfit)$cond$ID[1,1])/num_timepoints), 
  #             ranef(results[[i]]$TMBfit)$cond$ID[1,1] + 1.96*sqrt(exp(-ranef(results[[i]]$TMBfit)$cond$ID[1,1])/num_timepoints))
  # if (results[[i]]$true_alphas[1,1]<norm_PI[2]&results[[i]]$true_alphas[1,1]>norm_PI[1]) {
  #   normcount=normcount+1
  # }
  
  
  norm_PI = c(ranef(results[[i]]$TMBfit)$cond$ID[1,1] - 1.96*sqrt(exp(-results[[i]]$true_alphas[1,1])/num_timepoints), 
              ranef(results[[i]]$TMBfit)$cond$ID[1,1] + 1.96*sqrt(exp(-results[[i]]$true_alphas[1,1])/num_timepoints))
  if (results[[i]]$true_alphas[1,1]<norm_PI[2]&results[[i]]$true_alphas[1,1]>norm_PI[1]) {
    normcount=normcount+1
  }
  
  norm_intlength[i] = norm_PI[2] - norm_PI[1]
  
  
  mix_PI = c(ranef(results[[i]]$TMBfit)$cond$ID[1,1] - mix_quantiles[2] ,ranef(results[[i]]$TMBfit)$cond$ID[1,1] - mix_quantiles[1])
  
  if (results[[i]]$true_alphas[1,1]<mix_PI[2]&results[[i]]$true_alphas[1,1]>mix_PI[1]) {
    mixcount=mixcount+1
  }
  
  mix_intlength[i] = mix_PI[2] - mix_PI[1]
  
  
  #######glmmTMB intervals#########
  
  TMBalpha_int_sd = results[[i]]$TMBfit %>% ranef %>% as.data.frame %>% .$condsd %>% .[1]
  
  TMB_mix_PI1 = c(ranef(results[[i]]$TMBfit)$cond$ID[1,1] - 1.96*TMBalpha_int_sd , ranef(results[[i]]$TMBfit)$cond$ID[1,1] + 1.96*TMBalpha_int_sd)
  
  if (results[[i]]$true_alphas[1,1]<TMB_mix_PI1[2]&results[[i]]$true_alphas[1,1]>TMB_mix_PI1[1]) {
    TMB_mixcount=TMB_mixcount+1
  }
  
  TMB_intlength[i] = TMB_mix_PI1[2] - TMB_mix_PI1[1]
  
  
  
  #######lme4 intervals#########
  
  lme4alpha_int_sd = results[[i]]$lme4fit %>% ranef %>% as.data.frame %>% .$condsd %>% .[1]
  
  lme4_mix_PI1 = c(ranef(results[[i]]$lme4fit)$ID[1,1] - 1.96*lme4alpha_int_sd , ranef(results[[i]]$lme4fit)$ID[1,1] + 1.96*lme4alpha_int_sd)
  
  if (results[[i]]$true_alphas[1,1]<lme4_mix_PI1[2]&results[[i]]$true_alphas[1,1]>lme4_mix_PI1[1]) {
    lme4_mixcount=lme4_mixcount+1
  }
  
  lme4_intlength[i] = lme4_mix_PI1[2] - lme4_mix_PI1[1]
  
  
}

normcount/num_sims
mixcount/num_sims
TMB_mixcount/num_sims
lme4_mixcount/num_sims
naive_count/num_sims

x <- sapply(results, function(x){
  sqrt(exp(-ranef(x$TMBfit)$cond$ID[1,1]))/sqrt(num_timepoints)
})

y <- sapply(results, function(x){
  x$TMBfit %>% ranef %>% as.data.frame %>% .$condsd %>% .[1]
})

cbind(x,y) %>% as.data.frame %>% 
  ggplot(aes(x = x, y = y)) +
  geom_point() +
  theme_classic() +
  labs(title = 'Cluster One', x = bquote(exp(-hat(alpha)[1]/2)/sqrt(m)) , y = 'Estimated SD')+
  theme(plot.title = element_text(size = 12, hjust = 0.5)) +
  geom_abline(intercept = 0, slope = 1, col = 'red')


mean(mix_intlength)
mean(TMB_intlength)
mean(lme4_intlength)
mean(norm_intlength)
mean(naive_intlength)
var(mix_intlength)*num_timepoints
var(TMB_intlength)*num_timepoints
var(lme4_intlength)*num_timepoints
var(norm_intlength)*num_timepoints
var(naive_intlength)*num_timepoints


#try same thing for second cluster

x2 <- sapply(results, function(x){
  sqrt(exp(-ranef(x$TMBfit)$cond$ID[2,1]))/sqrt(num_timepoints)
})

y2 <- sapply(results, function(x){
  x$TMBfit %>% ranef %>% as.data.frame %>% .$condsd %>% .[2]
})

cbind(x2,y2) %>% as.data.frame %>% 
  ggplot(aes(x = x2, y = y2)) +
  geom_point() +
  theme_classic() +
  labs(title = 'Cluster Two', x = bquote(exp(-hat(alpha)[2]/2)/sqrt(m)) , y = 'Estimated SD')+
  theme(plot.title = element_text(size = 12, hjust = 0.5)) +
  geom_abline(intercept = 0, slope = 1, col = 'red')


#compare first and second cluster

cbind(y,y2) %>% as_tibble %>% 
  ggplot(aes(x = y, y = y2)) +
  geom_point() +
  theme_classic() + 
  labs(title = 'Estimated SDs', x = 'Cluster One', y = 'Cluster Two') +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

facet_names <- c('V1' = "One",'V2' = "Two")

cbind(log(num_timepoints*y^2),log(num_timepoints*y2^2)) %>% as_tibble %>% pivot_longer(V1:V2,names_to = 'cluster',values_to = 'sd_est') %>%
  ggplot(aes(x=sd_est,fill=cluster)) +
  geom_histogram(aes(y = ..density..), alpha = 0.4, binwidth = 0.3) +
  geom_line(aes(y = ..density.., colour = 'black'), stat = 'density') + 
  facet_wrap(cluster~., labeller = as_labeller(facet_names)) +
  labs(title = 'Log Estimated Variances (Normalised)', x = 'Cluster') +
  theme_classic()+
  theme(plot.title = element_text(size = 12, hjust = 0.5),legend.position="none")+
  scale_fill_brewer(type = 'qual', palette = 3) 


######### Results to show in MS ##########


mixcount/num_sims
TMB_mixcount/num_sims
naive_count/num_sims
mean(mix_intlength)
mean(TMB_intlength)
var(TMB_intlength)*num_timepoints
shapiro.test(TMBdiffs)


