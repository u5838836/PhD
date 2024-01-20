rm(list = ls())
library(glmmTMB)
library(dplyr)
library(xtable)
library(mvtnorm)
#library(MASS)
library(Matrix)
library(parallel)
library(foreach)
library(doParallel)
#library(tidyverse)
library(GGally)
library(lme4)
library(glm2)
library(ALA)
library(faraway)
library(R2BayesX)
library(tidyr)


########################-------------Bolus Data-------------###############################
library(cold)
dat = bolus
setwd("D:/Study new/PhD/R code/3rd/bolus")
dat$id = as.factor(dat$id)

glmmtmb_mod = glmmTMB(y  ~ 1 + time + (1 + time |id), data = dat,family = poisson)
glmmtmb_ranef = glmmtmb_mod %>% ranef %>% as.data.frame
glmmtmb_ranef
glmmtmb_mod
lme_mod = glmer(y  ~ 1 + time + (1 + time |id), data = dat,family = poisson)
lme_ranef = lme_mod %>% ranef %>% as.data.frame
lme_ranef
lme_mod

varcovar = VarCorr(glmmtmb_mod)$cond[[1]]
varcovar_lme = VarCorr(lme_mod)$id
X_1 = glmmtmb_mod$frame %>% filter(id=="1") %>% dplyr::select(time)
X_1 = rep(1,nrow(X_1)) %>% cbind(X_1) %>% as.matrix

#########---------- normal scale mixture quantiles ------------##########
S = matrix(nrow=10000,ncol=ncol(X_1))
for (i in 1:10000) {
  b_1 = rmvnorm(n=1,sigma = varcovar_lme) %>% t %>% as.vector
  eta_1 = X_1%*%lme_mod@beta + X_1%*%b_1
  W_1 = diag(nrow=nrow(X_1))
  diag(W_1) = poisson()$var(poisson()$linkinv(eta_1))
  S[i,] = rmvnorm(1,sigma = solve(crossprod(X_1,W_1)%*%X_1))
}
int_interval = quantile(S[,1],probs = c(0.025,0.975))
treat_interval = quantile(S[,2],probs = c(0.025,0.975))
lme_sds_int = lme_ranef %>% filter(term == '(Intercept)') %>% .$condsd
lme_sds_treat = lme_ranef %>% filter(term == 'time') %>% .$condsd
glmmtmb_sds_int = glmmtmb_ranef %>% filter(term == '(Intercept)') %>% .$condsd
glmmtmb_sds_treat = glmmtmb_ranef %>% filter(term == 'time') %>% .$condsd

patient = 1:65
all_groups = tibble(patient,lme_sds_int,glmmtmb_sds_int,lme_sds_treat,glmmtmb_sds_treat) %>%
  pivot_longer(c(lme_sds_int,glmmtmb_sds_int,lme_sds_treat,glmmtmb_sds_treat),
               names_to = "package",values_to = "estimate") %>%
  mutate(type = rep(c('Intercept','Intercept','Slope','Slope'),65)) %>%
  mutate(Software = rep(c('lme4','glmmTMB','lme4','glmmTMB'),65)) %>%
  mutate(intlength = 2*1.96*estimate) %>%
  mutate(unconditional = rep(c(int_interval[2]-int_interval[1],int_interval[2]-int_interval[1],
                               treat_interval[2]-treat_interval[1],treat_interval[2]-treat_interval[1]),65)) %>%
  mutate(paired = rep(1:130,each=2))

pdf('intlength.pdf', width=7, height=3.5)
all_groups %>%
  ggplot(aes(x = patient, y = intlength)) +
  facet_wrap(~ type, ncol = 2, scales = "free") +
  geom_point(aes(shape=Software, color = Software), size = 2) +
  labs(x = "Patient Number" , y = "Interval Length", col = "Package", shape = "Package") +
  geom_hline(aes(yintercept=unconditional)) +
  #geom_vline(data=filter(all_groups, type=="Intercept"), aes(xintercept = 49),linetype = 2) +
  #geom_vline(data=filter(all_groups, type=="Intercept"), aes(xintercept = 40),linetype = 2) +
  #geom_vline(data=filter(all_groups, type=="Slope"), aes(xintercept = 49),linetype = 2) +
  #geom_vline(data=filter(all_groups, type=="Slope"), aes(xintercept = 58),linetype = 2) +
  geom_line(aes(group = paired),linetype = 3) +
  scale_color_manual(name = "Package", labels = c("glmmTMB", "lme4"), values = c(2,4)) +
  scale_shape_manual(name = "Package", labels = c("glmmTMB", "lme4"), values = c(19,17)) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 65, by = 5))
dev.off()

plot(glmmtmb_ranef$condval,lme_ranef$condval,
     ylab = "lme4 Predicted Random Effect", xlab = "glmmTMB Predicted Random Effect")
abline(0,1,col="red")

##############---------epilepsy data--------------###############
#no age

# setwd("D:/Study new/PhD/R code/3rd/Figures")
# 
# data(epilepsy)
# epilepsy$id = as.factor(epilepsy$id)
# epilepsy$treated = epilepsy$treat*epilepsy$expind
# epilepsy$treated_f = as.factor(epilepsy$treated)
# epilepsy$log_tadj = log(epilepsy$timeadj)
# epilepsy$seize_adj = epilepsy$seizures/epilepsy$timeadj
# epilepsy$period = rep(1:5)
# epilepsy$treat_f = as.factor(epilepsy$treat)
# #epilepsy = epilepsy %>% filter(treat==1)
# 
# #pdf('boxplot.pdf')
# epilepsy %>% ggplot(aes(x = treated_f,y = seize_adj)) +
#   geom_boxplot() +
#   labs(x = "Treated" , y = "Seizures Per Week") +
#   scale_x_discrete(labels = c("No","Yes")) +
#   theme_bw() 
# #dev.off()
# 
# #pdf('ts.pdf')
# epilepsy %>% group_by(treat_f,period) %>% summarise(mean = mean(seize_adj)) %>%
#   ggplot(aes(x=period,y=mean)) +
#   geom_line(aes(col = treat_f)) +
#   labs(x = "Period" , y = "Mean Seizures Per Week", col = "Group") +
#   scale_color_manual(labels = c("Control", "Treatment"), values = c("2", "4")) +
#   theme_bw()
# #dev.off()
# 
# epilepsy %>% group_by(id,period) %>% summarise(mean = mean(seize_adj)) %>%
#   ggplot(aes(x=period,y=mean)) +
#   geom_line(aes(col = id)) +
#   theme_bw()
# 
# glmmtmb_mod = glmmTMB(seizures  ~ 1 + treated + (1 + treated |id),
#                       data = epilepsy,family = poisson, offset = log_tadj)
# glmmtmb_ranef = glmmtmb_mod %>% ranef %>% as.data.frame
# glmmtmb_ranef
# glmmtmb_mod
# lme_mod = glmer(seizures  ~ 1 + treated + (1 + treated |id),
#                 data = epilepsy,family = poisson, offset = log_tadj)
# lme_ranef = lme_mod %>% ranef %>% as.data.frame
# lme_ranef
# lme_mod
# 
# varcovar = VarCorr(glmmtmb_mod)$cond[[1]]
# varcovar_lme = VarCorr(lme_mod)$id
# X_1 = glmmtmb_mod$frame %>% filter(id=="29") %>% select(treated)
# X_1 = rep(1,nrow(X_1)) %>% cbind(X_1) %>% as.matrix
# glmmtmb_mod$fit$parfull
# sigma(glmmtmb_mod)
# lme_mod@beta
# lme_mod
# summary(lme_mod)
# small_offset = glmmtmb_mod$frame %>% filter(id=="29") %>% select("(offset)")
# small_offset = small_offset$`(offset)`
# 
# #########---------- normal scale mixture quantiles ------------##########
# S = matrix(nrow=10000,ncol=ncol(X_1))
# for (i in 1:10000) {
#   b_1 = rmvnorm(n=1,sigma = varcovar_lme) %>% t %>% as.vector
#   eta_1 = X_1%*%lme_mod@beta + X_1%*%b_1 + small_offset
#   W_1 = diag(nrow=nrow(X_1))
#   diag(W_1) = poisson()$var(poisson()$linkinv(eta_1))
#   S[i,] = rmvnorm(1,sigma = solve(crossprod(X_1,W_1)%*%X_1))
# }
# int_interval = quantile(S[,1],probs = c(0.025,0.975))
# treat_interval = quantile(S[,2],probs = c(0.025,0.975))
# lme_sds_int = lme_ranef %>% filter(grp %in% 29:59 & term == '(Intercept)') %>% .$condsd 
# lme_sds_treat = lme_ranef %>% filter(grp %in% 29:59 & term == 'treated') %>% .$condsd 
# glmmtmb_sds_int = glmmtmb_ranef %>% filter(grp %in% 29:59 & term == '(Intercept)') %>% .$condsd 
# glmmtmb_sds_treat = glmmtmb_ranef %>% filter(grp %in% 29:59 & term == 'treated') %>% .$condsd 
# 
# 
# # pdf('intlength.pdf')
# # 
# # plot(lme_sds_int*1.96*2,ylab='Interval Length',xlab='Patient Number')
# # abline(h=int_interval[2]-int_interval[1])
# # dev.off()
# 
# patient = 29:59
# 
# treatment_group = tibble(patient,lme_sds_int,glmmtmb_sds_int,lme_sds_treat,glmmtmb_sds_treat) %>% 
#   pivot_longer(c(lme_sds_int,glmmtmb_sds_int,lme_sds_treat,glmmtmb_sds_treat),
#                names_to = "package",values_to = "estimate") %>%
#   mutate(type = rep(c('Intercept','Intercept','Slope','Slope'),31)) %>%
#   mutate(Software = rep(c('lme4','glmmTMB','lme4','glmmTMB'),31)) %>%
#   mutate(intlength = 2*1.96*estimate) %>%
#   mutate(unconditional = rep(c(int_interval[2]-int_interval[1],int_interval[2]-int_interval[1],
#                                treat_interval[2]-treat_interval[1],treat_interval[2]-treat_interval[1]),31)) %>%
#   mutate(paired = rep(1:62,each=2))
# 
# #pdf('intlength.pdf', width=7, height=7)
# tibble(patient,lme_sds_int,glmmtmb_sds_int,lme_sds_treat,glmmtmb_sds_treat) %>% 
#   pivot_longer(c(lme_sds_int,glmmtmb_sds_int,lme_sds_treat,glmmtmb_sds_treat),
#                names_to = "package",values_to = "estimate") %>%
#   mutate(type = rep(c('Intercept','Intercept','Slope','Slope'),31)) %>%
#   mutate(Software = rep(c('lme4','glmmTMB','lme4','glmmTMB'),31)) %>%
#   mutate(intlength = 2*1.96*estimate) %>%
#   mutate(unconditional = rep(c(int_interval[2]-int_interval[1],int_interval[2]-int_interval[1],
#                                treat_interval[2]-treat_interval[1],treat_interval[2]-treat_interval[1]),31)) %>%
#   mutate(paired = rep(1:62,each=2)) %>%
#   ggplot(aes(x = patient, y = intlength)) +
#   facet_wrap(~ type, ncol = 1, scales = "free") +
#   geom_point(aes(shape=Software), size = 2)+
#   labs(x = "Patient Number" , y = "Interval Length", col = "Package") +
#   geom_hline(aes(yintercept=unconditional)) +
#   geom_vline(data=filter(treatment_group, type=="Intercept"), aes(xintercept = 49),linetype = 2) + 
#   geom_vline(data=filter(treatment_group, type=="Intercept"), aes(xintercept = 40),linetype = 2) + 
#   geom_vline(data=filter(treatment_group, type=="Slope"), aes(xintercept = 49),linetype = 2) + 
#   geom_vline(data=filter(treatment_group, type=="Slope"), aes(xintercept = 58),linetype = 2) + 
#   geom_line(aes(group = paired),linetype = 3) +
#   scale_color_manual(labels = c("glmmTMB", "lme4"), values = c(2,4)) +
#   theme_classic() +
#   scale_x_continuous(breaks = seq(29, 59, by = 1))
# #dev.off()
# 
# 
# plot(lme_sds_treat*1.96*2,ylab='Interval Length',xlab='Patient')
# abline(h=treat_interval[2]-treat_interval[1])
# mean((lme_sds_int*1.96*2)^2)
# (int_interval[2]-int_interval[1])^2
# mean((lme_sds_treat*1.96*2)^2)
# (treat_interval[2]-treat_interval[1])^2
# 
# 
# #Control group - no unconditional results
# 
# lme_sds_int0 = lme_ranef %>% filter(grp %in% 1:28 & term == '(Intercept)') %>% .$condsd 
# lme_sds_treat0 = lme_ranef %>% filter(grp %in% 1:28 & term == 'treated') %>% .$condsd 
# glmmtmb_sds_int0 = glmmtmb_ranef %>% filter(grp %in% 1:28 & term == '(Intercept)') %>% .$condsd 
# glmmtmb_sds_treat0 = glmmtmb_ranef %>% filter(grp %in% 1:28 & term == 'treated') %>% .$condsd 
# 
# 
# patient0 = 1:28
# 
# #pdf('intlength0.pdf', width=7, height=7)
# tibble(patient0,lme_sds_int0,glmmtmb_sds_int0,lme_sds_treat0,glmmtmb_sds_treat0) %>% 
#   pivot_longer(c(lme_sds_int0,glmmtmb_sds_int0,lme_sds_treat0,glmmtmb_sds_treat0),
#                names_to = "package",values_to = "estimate") %>%
#   mutate(type = rep(c('Intercept','Intercept','Slope','Slope'),28)) %>%
#   mutate(Software = rep(c('lme4','glmmTMB','lme4','glmmTMB'),28)) %>%
#   mutate(intlength = 2*1.96*estimate) %>%
#   mutate(paired = rep(1:56,each=2)) %>%
#   ggplot(aes(x = patient0, y = intlength)) +
#   facet_wrap(~ type, ncol = 1, scales = "free") +
#   geom_point(aes(shape=Software), size = 2) +
#   geom_line(aes(group = paired), linetype = 3) +
#   labs(x = "Patient Number" , y = "Interval Length", col = "Package") +
#   scale_color_manual(labels = c("glmmTMB", "lme4"), values = c(2,4)) +
#   theme_classic() +
#   scale_x_continuous(breaks = seq(1, 28, by = 1))
# #dev.off()
# 
# plot(glmmtmb_ranef$condval,lme_ranef$condval, 
#      ylab = "lme4 Predicted Random Effect", xlab = "glmmTMB Predicted Random Effect")
# abline(0,1,col="red")




# Grouse ticks data
# data(grouseticks)
# grouseticks %>% group_by(LOCATION) %>% summarise(len = length(LOCATION)) %>% View
# grouseticks %>% group_by(BROOD) %>% summarise(len = length(BROOD)) %>% View
# glmmtmb_mod = glmmTMB(TICKS ~ 1 + (1|LOCATION),grouseticks,"poisson")
# lme4_mod = glmmTMB(TICKS ~ 1 + (1|LOCATION),grouseticks,"poisson")
# glmmtmb_mod1 = glmmTMB(TICKS ~ 1 + (1|BROOD),grouseticks,"poisson")
# lme4_mod1 = glmmTMB(TICKS ~ 1 + (1|BROOD),grouseticks,"poisson")
# 
# glmmtmb_mod %>% ranef %>% as.data.frame
# glmmtmb_mod1 %>% ranef %>% as.data.frame

########## AIDS DATA ############

#Check if for any two clusters i and j, X_i = X_j
# rm(list=ls())
# data(cd4)
# cd4$CD4 = ( exp(cd4$logCD4) - 1) %>% sapply(as.integer)
# cd4
# X_i = list()
# ind = 1
# for (i in levels(cd4$id)) {
#   X_i[[ind]] = cd4 %>% filter(id==i) %>% select(treatment, gender, week) #age,
#   ind = ind + 1
# }
# matequal <- function(x, y){dim(x) == dim(y) && all(x == y) }
# 
# pairwise_comp = matrix(nrow = length(X_i), ncol = length(X_i))
# for (i in 1:length(X_i)) {
#   for (j in i:length(X_i)) {
#     pairwise_comp[i,j] = matequal(X_i[[i]], X_i[[j]])
#   }
# }
# all_comp = pairwise_comp[upper.tri(pairwise_comp)]
# sum(all_comp)
# glmmtmb_mod = glmmTMB(CD4 ~ 1+ treatment + gender + week + (1 + treatment + gender + week |id),cd4,"poisson")
# glmmtmb_mod %>% ranef %>% as.data.frame

####### different possible datasets from faraway, ALA packagem lme4 ########

#many suffer from 'lack of (good) research question', some suffer 'model might not be good fit'
#ideally want the predictions (of the random effects) to be meaningful also

# data(lead) #gamma glmm possible here. ALA
# data(muscatine) #possible binomial, use only age as covariate. But base age + time is probably the more optimal model. ALA
# data(cake) #possible...gamma/poisson. But research question a little silly? lme4
# data(Penicillin) #poisson/gamma might be bad fit. lme4
# data(sleepstudy) #possible gamma. lme4
# 
# 
# #not many observations
# data(Dyestuff) #poisson/gamma might be bad fit, not many obs.
# data(Pastes) #gamma might be bad fit, not many obs
# data(cbpp) #possible if size covariate ignored, but I think size is the point of the data
# 
# 
# 
# 
# 
# glmmtmb_mod = glmmTMB(lead ~ 1 + week + (1 + week|id),lead,Gamma(link="log"))
# glmmtmb_mod %>% ranef %>% as.data.frame
# glmmtmb_mod
# 
# 
# 
# 
# Penicillin = cbind(Penicillin,model.matrix( ~ sample, Penicillin)) 
# glmmtmb_mod = glmmTMB(diameter ~ 1 + sampleB + sampleC + sampleD + sampleE + sampleF +
#                         (1 + sampleB + sampleC + sampleD + sampleE + sampleF|plate),Penicillin,Gamma,
#                       control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3),
#                                                optimizer=optim,optArgs=list(method="BFGS")))
# glmmtmb_mod %>% ranef %>% as.data.frame
# glmmtmb_mod
# varcovar = VarCorr(glmmtmb_mod)$cond[[1]]
# X_1 = glmmtmb_mod$frame %>% filter(plate=="a") %>% select(sampleB,sampleC,sampleD,sampleE,sampleF)
# X_1 = rep(1,nrow(X_1)) %>% cbind(X_1) %>% as.matrix
# glmmtmb_mod$fit$parfull
# sigma(glmmtmb_mod)
# 
# #########---------- normal scale mixture quantiles ------------##########
# S = matrix(nrow=10000,ncol=ncol(X_1))
# for (i in 1:10000) {
#   b_1 = rmvnorm(n=1,sigma = varcovar) %>% t %>% as.vector
#   eta_1 = X_1%*%glmmtmb_mod$fit$par[1:ncol(X_1)] + X_1%*%b_1
#   W_1 = diag(as.vector(1/(eta_1^2*sigma(glmmtmb_mod)^2)))
#   S[i,] = rmvnorm(1,sigma = solve(crossprod(X_1,W_1)%*%X_1))
# }
# quantile(S[,1],probs = c(0.025,0.975))
# quantile(S[,2],probs = c(0.025,0.975))
# quantile(S[,3],probs = c(0.025,0.975))
# quantile(S[,4],probs = c(0.025,0.975))
# quantile(S[,5],probs = c(0.025,0.975))
# quantile(S[,6],probs = c(0.025,0.975))
# 
# 
# 
# glmmtmb_mod = glmmTMB(Reaction ~ 1 + Days + (1 + Days|Subject),sleepstudy,Gamma(link="log"))
# lme4_mod = glmer(Reaction ~ 1 + Days + (1 + Days|Subject),sleepstudy,Gamma)
# glmmtmb_mod %>% ranef %>% as.data.frame
# lme4_mod %>% ranef %>% as.data.frame
# glmmtmb_mod
# varcovar = VarCorr(glmmtmb_mod)$cond[[1]]
# X_1 = glmmtmb_mod$frame %>% filter(Subject=="308") %>% select(Days)
# X_1 = rep(1,nrow(X_1)) %>% cbind(X_1) %>% as.matrix
# glmmtmb_mod$fit$parfull
# sigma(glmmtmb_mod)
# 
# #########---------- normal scale mixture quantiles ------------##########
# S = matrix(nrow=10000,ncol=ncol(X_1))
# for (i in 1:10000) {
#   b_1 = rmvnorm(n=1,sigma = varcovar) %>% t %>% as.vector
#   eta_1 = X_1%*%glmmtmb_mod$fit$par[1:2] + X_1%*%b_1
#   W_1 = diag(nrow=nrow(X_1))
#   diag(W_1) = (exp(2*eta_1)/(sigma(glmmtmb_mod)^2*exp(2*eta_1)))
#   S[i,] = rmvnorm(1,sigma = solve(crossprod(X_1,W_1)%*%X_1))
#   }
# quantile(S[,1],probs = c(0.025,0.975))
# quantile(S[,2],probs = c(0.025,0.975))
# 
# #############
# 
# glmmtmb_mod = glmmTMB(angle  ~ 1 + recipe + temperature + (1 + recipe + temperature|replicate),cake,Gamma(link="log"))
# glmmtmb_mod %>% ranef %>% as.data.frame
# glmmtmb_mod
# 
# 
# glmmtmb_mod = glmmTMB(obesity  ~ 1 + year + (1 + year|id),muscatine[1:100,],binomial)
# glmmtmb_mod %>% ranef %>% as.data.frame
# glmmtmb_mod

##############---------epilepsy data-------------- (faraway) ###############


# data(epilepsy)
# epilepsy$id = as.factor(epilepsy$id)
# epilepsy$age_adj = rep(c(0,8,10,12,14))
# epilepsy[,'age_w'] = epilepsy$age*52 + epilepsy$age_adj
# epilepsy[,'age_w'] = (epilepsy[,'age_w'] - mean(epilepsy[,'age_w']))/sqrt(var(epilepsy[,'age_w']))
# epilepsy$treated = epilepsy$treat*epilepsy$expind
# epilepsy$log_tadj = log(epilepsy$timeadj)
# glmmtmb_mod = glmmTMB(seizures  ~ 1 + treated + age_w + (1 + treated + age_w|id),
#                       data = epilepsy,family = poisson, offset = log_tadj)
# glmmtmb_mod %>% ranef %>% as.data.frame
# glmmtmb_mod
# lme_mod = glmer(seizures  ~ 1 + treated + age_w + (1 + treated + age_w|id),
#                       data = epilepsy,family = poisson, offset = log_tadj)
# lme_mod %>% ranef %>% as.data.frame
# lme_mod
# 
# X_i = list()
# ind = 1
# for (i in levels(epilepsy$id)) {
#   X_i[[ind]] = epilepsy %>% filter(id==i) %>% select(age)
#   ind = ind + 1
# }
# matequal <- function(x, y){dim(x) == dim(y) && all(x == y) }
# 
# pairwise_comp = matrix(nrow = length(X_i), ncol = length(X_i))
# for (i in 1:length(X_i)) {
#   for (j in i:length(X_i)) {
#     pairwise_comp[i,j] = matequal(X_i[[i]], X_i[[j]])
#   }
# }
# diag(pairwise_comp) = FALSE
# all_comp = pairwise_comp[upper.tri(pairwise_comp)]
# sum(all_comp)
# pairwise_comp_treated = pairwise_comp[29:59,29:59]
# all_comp_treated = pairwise_comp_treated[upper.tri(pairwise_comp_treated)]
# sum(all_comp_treated)
# 
# for (i in 1:nrow(pairwise_comp_treated)) {
#   for (j in i:ncol(pairwise_comp_treated)) {
#     if (pairwise_comp_treated[i,j]==TRUE) {
#       print(c(i,j)+28)
#     }
#   }
# }
# ordered_epilepsy = epilepsy[order(epilepsy$age,epilepsy$id),]
# ordered_epilepsy_treat = ordered_epilepsy %>% filter(treat==1)
# 
# varcovar = VarCorr(glmmtmb_mod)$cond[[1]]
# varcovar_lme = VarCorr(lme_mod)$id
# X_1 = glmmtmb_mod$frame %>% filter(id=="29") %>% select(treated,age_w,"(offset)")
# X_1 = rep(1,nrow(X_1)) %>% cbind(X_1) %>% as.matrix
# glmmtmb_mod$fit$parfull
# sigma(glmmtmb_mod)
# lme_mod@beta
# lme_mod
# summary(lme_mod)
# 
# #########---------- normal scale mixture quantiles ------------##########
# S = matrix(nrow=10000,ncol=ncol(X_1)-1)
# for (i in 1:10000) {
#   b_1 = rmvnorm(n=1,sigma = varcovar_lme) %>% t %>% as.vector
#   eta_1 = X_1%*%c(lme_mod@beta,1) + X_1[,-ncol(X_1)]%*%b_1
#   W_1 = diag(nrow=nrow(X_1))
#   diag(W_1) = poisson()$var(poisson()$linkinv(eta_1))
#   S[i,] = rmvnorm(1,sigma = solve(crossprod(X_1[,-ncol(X_1)],W_1)%*%X_1[,-ncol(X_1)]))
#   }
# quantile(S[,1],probs = c(0.025,0.975))









