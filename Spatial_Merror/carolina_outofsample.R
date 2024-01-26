rm(list = ls())
library(raster)
library(tidyverse)
library(GGally)
library(mgcv)
library(sp)
library(TMB)
library(mvtnorm)
library(autoFRK)
library(spaMM)
library(RandomFields)
library(ggplot2)
library(forcats)
library(fields)
library(optimx)
library(FRK)
library(Matrix)
library(MASS)
library(glm2)
library(xtable)
library(plotly)
library(RColorBrewer)
library(rgdal)
library(simex)
library(ROCR)

setwd("D:/Study new/PhD/R code/Spatial_Merror")
#gdbsource("poisson2D_2step.R",interactive=TRUE)
TMB::compile("multX_bin.cpp","-O1 -g",DLLFLAGS="")
dyn.load(dynlib("multX_bin"))


##-----------------
## Two predictors
##-----------------
rmin = raster("us_tmin.txt")/100 # To celsius. Might be worth dividing by 1000 to examine changes in tens of a degree
rmax = raster("us_tmax.txt")/100 # To celsius. Might be worth dividing by 1000 to examine changes in tens of a degree
r <- raster::stack(rmin, rmax)

coords <- read.csv("2010_coord.csv") %>% 
  dplyr::select(Longitude, Latitude)
#coords = coords[-502,]
coordinates(coords) <- c("Longitude" , "Latitude")
proj4string(coords) <- CRS("+proj=longlat +datum=WGS84")
UTM <- spTransform(coords, CRS("+proj=utm +zone=16 ellps=WGS84"))



rasValue <- raster::extract(r, coords) 
combinePointValue <- cbind(UTM@coords, rasValue) %>% as.data.frame
combinePointValue %>% 
  head

##-----------------
## Response
##-----------------
resp <- read.csv("Wren_PA_2010.csv")
combinePointValue$resp <- resp$x #[-502]
combinePointValue$index <- 1:nrow(combinePointValue)


##-----------------
## Explore data & processing
##-----------------
combinePointValue %>% 
  head

#ggpairs(combinePointValue)


presences = combinePointValue %>% filter(resp == 1) %>% arrange(Latitude)
n_presences = nrow(presences)
presence_test.index = seq(1,n_presences,by = 5)

absences = combinePointValue %>% filter(resp == 0) %>% arrange(Latitude)
n_absences = nrow(absences)
absence_test.index = seq(1,n_absences,by = 5)

combinePointValue_train = rbind(presences[-presence_test.index,],absences[-absence_test.index,])
train_ind = combinePointValue_train$index
combinePointValue_test = combinePointValue[-train_ind,]

#Some plots

#pdf('temp_data_train.pdf')

par(mfrow = c(2,2), las = 1)
plot(rmin, main = "Min Temp.", xlim=c(-101,-65), xlab = "Longitude", ylab = "Latitude")
points(coords$Longitude[train_ind],coords$Latitude[train_ind],col=combinePointValue_train$resp, cex = 0.2, pch = 19)
legend(-77,32,legend = c("Absence","Presence"),col = 0:1,pch=19, bg = "green")
plot(rmax, main = "Max Temp.", xlim=c(-101,-65), xlab = "Longitude", ylab = "Latitude")
points(coords$Longitude[train_ind],coords$Latitude[train_ind],col=combinePointValue_train$resp, cex = 0.2, pch = 19)
legend(-77,32,legend = c("Absence","Presence"),col = 0:1,pch=19, bg = "green")

#rm(rmin, rmax, r, rasValue)

rmin_pi = raster("us_tmin_pi.txt")/20 #%>% sqrt
rmax_pi = raster("us_tmax_pi.txt")/20 #%>% sqrt
PI_raster = stack(rmin_pi,rmax_pi)
PI_obs = raster::extract(PI_raster,coords)
PI_obs[882,1] = mean(PI_obs[,1], na.rm = TRUE) #881 if outlier excluded, 882 if outlier included
plot(rmin_pi, main = "Min Temp. Error", xlim=c(-101,-65), xlab = "Longitude", ylab = "Latitude")
plot(rmax_pi, main = "Max Temp. Error", xlim=c(-101,-65), xlab = "Longitude", ylab = "Latitude")

#dev.off()

#############-------------- Fit models and compare ------------#############

n = length(combinePointValue_train$resp)

sizes = rep(1,n)

xbf_num = 200
xsqbf_num = 200
x2bf_num = 150
spbf_num = 15

x_basisfunctions <- mrts(UTM@coords[train_ind,], k = xbf_num, x = UTM@coords[train_ind,]) %>%
  as.matrix

x2_basisfunctions <- mrts(UTM@coords[train_ind,], k = x2bf_num, x = UTM@coords[train_ind,]) %>%
  as.matrix

x3_basisfunctions <- mrts(UTM@coords[train_ind,], k = xsqbf_num, x = UTM@coords[train_ind,]) %>%
  as.matrix

sp_basisfunctions <- mrts(UTM@coords[train_ind,], k = spbf_num, x = UTM@coords[train_ind,]) %>%
  as.matrix
sp_basisfunctions <- sp_basisfunctions[,-1] #no intercept


#Fit the model for W via autoFRK()

W = combinePointValue_train$us_tmin
W2 = combinePointValue_train$us_tmax
W3 = combinePointValue_train$us_tmin^2

fit_Wmodel <- autoFRK(Data = matrix(W,ncol=1), loc = UTM@coords[train_ind,], G = x_basisfunctions)
tmin_pred <- predict(fit_Wmodel)$pred

fit_W2model <- autoFRK(Data = matrix(W2,ncol=1), loc = UTM@coords[train_ind,], G = x2_basisfunctions)
tmax_pred <- predict(fit_W2model)$pred

fit_W3model <- autoFRK(Data = matrix(W3,ncol=1), loc = UTM@coords[train_ind,], G = x3_basisfunctions)
tmin2_pred <- predict(fit_W3model)$pred

combinePointValue_train$tmin_pred = tmin_pred
combinePointValue_train$tmax_pred = tmax_pred

mean((tmin_pred-W)^2)
mean((tmax_pred-W2)^2)
PI_obs %>% apply(2,function(x) mean(x,na.rm = TRUE))

sigmax_hat <- fit_Wmodel$M
alphax_hat <- fit_Wmodel$w
phi_hat <- fit_Wmodel$s

sigmax2_hat <- fit_W2model$M
alphax2_hat <- fit_W2model$w
phi_hat2 <- fit_W2model$s

sigmax3_hat <- fit_W3model$M
alphax3_hat <- fit_W3model$w
phi_hat3 <- fit_W3model$s


###################-------------Quadratic models with min temp only------------####################

#Naive


naive_mod = glm(resp ~ us_tmin + I(us_tmin^2) + sp_basisfunctions , family = binomial, data =  combinePointValue_train)

no.spatial_mod = glm(resp ~ tmin_pred + I(tmin_pred^2) , family = binomial, data = combinePointValue_train)

cholSigma = diag(log(abs(naive_mod$coefficients[-(1:3)])))[,1:1]
cholSigma_vec = cholSigma[lower.tri(cholSigma, diag = TRUE)]

inputdata=list(sp_basis = sp_basisfunctions , y = as.matrix(combinePointValue_train$resp), 
               predX = combinePointValue_train$us_tmin, predX2 = combinePointValue_train$us_tmin^2, true_rank = 1, 
               identity_matrix = diag(nrow=ncol(sp_basisfunctions)), sizes = sizes) 

parameterlist=list(beta = naive_mod$coefficients[1:3] , cholSigma = cholSigma_vec , 
                   alpha_rho = naive_mod$coefficients[-(1:3)])


obj <- MakeADFun(data=inputdata,DLL = "multX_bin",
                 parameters = parameterlist, random = c("alpha_rho"), silent=T)


fit1 <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr,
               control=list(trace = 0, iter.max = 1000, eval.max = 1000))


report1 = sdreport(obj,getJointPrecision = TRUE, getReportCovariance = FALSE, skip.delta.method = TRUE)
report_summary1=summary(report1)
params1 = report_summary1 %>% 
  as.data.frame %>% 
  filter(rownames(report_summary1)%in%c("alpha_rho","beta")) %>%
  dplyr::select("Estimate") %>% .$Estimate

report_summary1[1:3,]



#Corrected

inputdata=list(sp_basis = sp_basisfunctions , y = as.matrix(combinePointValue_train$resp), 
               predX = tmin_pred, predX2 = tmin_pred^2, true_rank = 1, 
               identity_matrix = diag(nrow=ncol(sp_basisfunctions)), sizes = sizes) 

parameterlist=list(beta = naive_mod$coefficients[1:3] , cholSigma = cholSigma_vec , 
                   alpha_rho = naive_mod$coefficients[-(1:3)])


obj <- MakeADFun(data=inputdata,DLL = "multX_bin",
                 parameters = parameterlist, random = c("alpha_rho"), silent=T)


fit2 <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr,
               control=list(trace = 0, iter.max = 1000, eval.max = 1000))

report2 = sdreport(obj,getJointPrecision = TRUE, getReportCovariance = FALSE, skip.delta.method = TRUE)
report_summary2=summary(report2)
params2 = report_summary2 %>% 
  as.data.frame %>% 
  filter(rownames(report_summary2)%in%c("alpha_rho","beta")) %>%
  dplyr::select("Estimate") %>% .$Estimate

sigmarho_hat = as.list(report2,"Est",TRUE)$Sigma

#Non-parametric estimate of variance of score

varscore_x = 0
varscore_y = 0

eta = params2[1] + params2[2]*tmin_pred + params2[3]*tmin_pred^2 + sp_basisfunctions%*%params2[-(1:3)] 
iterative_weights = sizes*binomial()$var(binomial()$linkinv(eta)) %>% as.vector %>% diag.spam

for (j in 1:n) {
  varscore_x = varscore_x + tcrossprod(x_basisfunctions[j,]*(W[j] - tmin_pred[j])/phi_hat ) 
  
  score_y0 = combinePointValue_train$resp[j] - sizes[j]*binomial()$linkinv(eta[j]) 
  
  score_y1 = tmin_pred[j]*(combinePointValue_train$resp[j] - sizes[j]*binomial()$linkinv(eta[j])) 
  
  score_y2 = tmin_pred[j]^2*(combinePointValue_train$resp[j] - sizes[j]*binomial()$linkinv(eta[j]))
  
  score_y_rho = sp_basisfunctions[j,]*(combinePointValue_train$resp[j] - sizes[j]*binomial()$linkinv(eta[j])) 
  score_y = c(score_y0,score_y1,score_y2,score_y_rho)
  varscore_y = varscore_y + tcrossprod(score_y)
}



varscore_y[-(1:3),-(1:3)] = varscore_y[-(1:3),-(1:3)] + ginv(sigmarho_hat,1*(.Machine$double.eps)^(1/4))
varscore_x = varscore_x + ginv(sigmax_hat)

#observed negative Hessians (joint)

Xstar = cbind(1, tmin_pred, tmin_pred^2, sp_basisfunctions)
Y_hess = crossprod(Xstar,iterative_weights)%*%Xstar + 
  bdiag.spam(matrix(0,nrow=3,ncol=3),ginv(sigmarho_hat,1*(.Machine$double.eps)^(1/4)))

X_hess = crossprod(x_basisfunctions)/phi_hat + ginv(sigmax_hat) #X part/block

XY_hess = fit2$par[2]*crossprod(Xstar,iterative_weights)%*%x_basisfunctions + 
  2*fit2$par[3]*crossprod(Xstar,iterative_weights)%*%diag(c(tmin_pred))%*%x_basisfunctions

jointPrecision = bdiag(Y_hess, X_hess)
jointPrecision[1:(3+ncol(sp_basisfunctions)),(4+ncol(sp_basisfunctions)):ncol(jointPrecision)] =
  XY_hess


varscore = bdiag.spam(varscore_y,varscore_x)
invnegHess = as.matrix(solve(jointPrecision))
varcovar = tcrossprod(invnegHess%*%varscore,invnegHess)
together_vars = c(varcovar[1,1],varcovar[2,2],varcovar[3,3])



#Corrected V2 (separate basis functions for X^2)


inputdata=list(sp_basis = sp_basisfunctions , y = as.matrix(combinePointValue_train$resp), 
               predX = tmin_pred, predX2 = tmin2_pred, true_rank = 1, 
               identity_matrix = diag(nrow=ncol(sp_basisfunctions)), sizes = sizes) 

parameterlist=list(beta = naive_mod$coefficients[1:3] , cholSigma = cholSigma_vec , 
                   alpha_rho = naive_mod$coefficients[-(1:3)])


obj <- MakeADFun(data=inputdata,DLL = "multX_bin",
                 parameters = parameterlist, random = c("alpha_rho"), silent=T)


fit3 <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr,
               control=list(trace = 0, iter.max = 1000, eval.max = 1000))

report3 = sdreport(obj,getJointPrecision = TRUE, getReportCovariance = FALSE, skip.delta.method = TRUE)
report_summary3=summary(report3)
params3 = report_summary3 %>% 
  as.data.frame %>% 
  filter(rownames(report_summary3)%in%c("alpha_rho","beta")) %>%
  dplyr::select("Estimate") %>% .$Estimate

sigmarho_hat = as.list(report3,"Est",TRUE)$Sigma

#Non-parametric estimate of variance of score

varscore_x = 0
varscore_x2 = 0
varscore_y = 0

eta = params2[1] + params2[2]*tmin_pred + params2[3]*tmin2_pred + sp_basisfunctions%*%params2[-(1:3)] 
iterative_weights = sizes*binomial()$var(binomial()$linkinv(eta)) %>% as.vector %>% diag.spam

for (j in 1:n) {
  varscore_x = varscore_x + tcrossprod(x_basisfunctions[j,]*(W[j] - tmin_pred[j])/phi_hat ) 
  
  varscore_x2 = varscore_x2 + tcrossprod(x3_basisfunctions[j,]*(W3[j] - tmin2_pred[j])/phi_hat3 )
  
  score_y0 = combinePointValue_train$resp[j] - sizes[j]*binomial()$linkinv(eta[j]) 
  
  score_y1 = tmin_pred[j]*(combinePointValue_train$resp[j] - sizes[j]*binomial()$linkinv(eta[j])) 
  
  score_y2 = tmin2_pred[j]*(combinePointValue_train$resp[j] - sizes[j]*binomial()$linkinv(eta[j]))
  
  score_y_rho = sp_basisfunctions[j,]*(combinePointValue_train$resp[j] - sizes[j]*binomial()$linkinv(eta[j])) 
  score_y = c(score_y0,score_y1,score_y2,score_y_rho)
  varscore_y = varscore_y + tcrossprod(score_y)
}



varscore_y[-(1:3),-(1:3)] = varscore_y[-(1:3),-(1:3)] + ginv(sigmarho_hat,1*(.Machine$double.eps)^(1/4))
varscore_x = varscore_x + ginv(sigmax_hat)
varscore_x2 = varscore_x2 + ginv(sigmax3_hat)

#observed negative Hessians (joint)

Xstar = cbind(1, tmin_pred, tmin2_pred, sp_basisfunctions)
Y_hess = crossprod(Xstar,iterative_weights)%*%Xstar + 
  bdiag.spam(matrix(0,nrow=3,ncol=3),ginv(sigmarho_hat,1*(.Machine$double.eps)^(1/4)))

X_hess = crossprod(x_basisfunctions)/phi_hat + ginv(sigmax_hat) #X part/block
X2_hess = crossprod(x3_basisfunctions)/phi_hat3 + ginv(sigmax3_hat)

XY_hess = cbind(fit3$par[2]*crossprod(Xstar,iterative_weights)%*%x_basisfunctions,
                fit3$par[3]*crossprod(Xstar,iterative_weights)%*%x3_basisfunctions)

jointPrecision = bdiag(Y_hess, X_hess, X2_hess)
jointPrecision[1:(3+ncol(sp_basisfunctions)),(4+ncol(sp_basisfunctions)):ncol(jointPrecision)] =
  XY_hess


varscore = bdiag.spam(varscore_y,varscore_x, varscore_x2)
invnegHess = as.matrix(solve(jointPrecision))
varcovar = tcrossprod(invnegHess%*%varscore,invnegHess)
separate_vars = c(varcovar[1,1],varcovar[2,2],varcovar[3,3])



################--------------Simex etc------------##################

uGAM = gam(resp ~ us_tmin + I(us_tmin^2) + s(Longitude,Latitude) , family = binomial, data =  combinePointValue_train)

set.seed(2464)
simexGAM <- simex(uGAM, SIMEXvariable = c("us_tmin"),
                  measurement.error = PI_obs[train_ind,'us_tmin_pi'], asymptotic = FALSE)

cGAM = gam(resp ~ tmin_pred + I(tmin_pred^2) + s(Longitude,Latitude) , family = binomial, data =  combinePointValue_train)

uGAMsumm = summary(uGAM)
simexGAM_summ = summary(simexGAM)
cGAMsumm = summary(cGAM)

dnaive = glm(resp ~ us_tmin + I(us_tmin^2) , family = binomial, data =  combinePointValue_train)

set.seed(2464)
simexNS <- simex(dnaive, SIMEXvariable = c("us_tmin"),
                 measurement.error = PI_obs[train_ind,'us_tmin_pi'], asymptotic = FALSE)

simexNS_summ = summary(simexNS)
summ_dnaive = summary(dnaive)
summ_no.spatial = summary(no.spatial_mod)


################-----------Out of sample prediction performance-----------------################


#Naive - uncorrected spatial

rho_preds <- as.matrix(mrts(UTM@coords[train_ind,], k = spbf_num, x = UTM@coords[-train_ind,]))[,-1]%*%params1[-(1:3)]
combinePointValue_test = combinePointValue_test %>% mutate(naive = 
   binomial()$linkinv(params1[1] + params1[2]*us_tmin + params1[3]*us_tmin^2 + rho_preds))


#Corrected

x_preds <- as.matrix(mrts(UTM@coords[train_ind,], k = xbf_num, x = UTM@coords[-train_ind,] ))%*%alphax_hat
rho_preds <- as.matrix(mrts(UTM@coords[train_ind,], k = spbf_num, x = UTM@coords[-train_ind,] ))[,-1]%*%params2[-(1:3)]
combinePointValue_test = combinePointValue_test %>% mutate(correct = 
    binomial()$linkinv(params2[1] + params2[2]*x_preds +params2[3]*x_preds^2 + rho_preds))


#Corrected V2

x_preds <- as.matrix(mrts(UTM@coords[train_ind,], k = xbf_num, x = UTM@coords[-train_ind,] ))%*%alphax_hat
x2_preds <- as.matrix(mrts(UTM@coords[train_ind,], k = xsqbf_num, x = UTM@coords[-train_ind,] ))%*%alphax3_hat
rho_preds <- as.matrix(mrts(UTM@coords[train_ind,], k = spbf_num, x = UTM@coords[-train_ind,] ))[,-1]%*%params3[-(1:3)]
combinePointValue_test = combinePointValue_test %>% mutate(correctV2 =  
  binomial()$linkinv(params3[1] + params3[2]*x_preds +params3[3]*x2_preds + rho_preds))


#dNaive


combinePointValue_test = combinePointValue_test %>% mutate(dNaive = binomial()$linkinv(
  dnaive$coefficients[1] + dnaive$coefficients[2]*us_tmin + dnaive$coefficients[3]*us_tmin^2))



#No spatial

x_preds <- as.matrix(mrts(UTM@coords[train_ind,], k = xbf_num, x = UTM@coords[-train_ind,] ))%*%alphax_hat
combinePointValue_test = combinePointValue_test %>% mutate(no.spatial = 
 binomial()$linkinv(no.spatial_mod$coefficients[1] +
 no.spatial_mod$coefficients[2]*x_preds +
 no.spatial_mod$coefficients[3]*x_preds^2))


#SimexNS

combinePointValue_test = combinePointValue_test %>% mutate(simexNS = 
 binomial()$linkinv(simexNS$coefficients[1] +
 simexNS$coefficients[2]*us_tmin +
 simexNS$coefficients[3]*us_tmin^2))


#uGAM


combinePointValue_test = combinePointValue_test %>% mutate(uGAM = 
  binomial()$linkinv(predict(uGAM,newdata = combinePointValue_test)))


#cGAM

newdat = combinePointValue_test %>% mutate(tmin_pred = x_preds)

combinePointValue_test = combinePointValue_test %>% mutate(cGAM = 
 binomial()$linkinv(predict(cGAM,newdata = newdat)))


#SimexGAM


combinePointValue_test = combinePointValue_test %>% mutate(simexGAM = 
 binomial()$linkinv(predict(simexGAM,newdata = combinePointValue_test)))


preds = combinePointValue_test %>% select(naive:simexGAM)

#MSE
MSE = apply(preds,2,function(x){
  (combinePointValue_test$resp - x) %>% .^2 %>% mean
})

#MAE

MAE = apply(preds,2,function(x){
  (combinePointValue_test$resp - x) %>% abs %>% mean
})

#zero-one error

ZOE = apply(preds,2,function(x){
  (combinePointValue_test$resp - round(x)) %>% abs %>% mean
})

#cross-entropy

CE = apply(preds,2,function(x){
  -(combinePointValue_test$resp*log(x) + (1-combinePointValue_test$resp)*log(1-x)) %>% mean
})


#AUC

AUC = apply(preds,2,function(x){
  pred <- prediction(x,combinePointValue_test$resp)
  perf <- performance(pred,"auc")
  perf@y.values[[1]]
})




rbind(MSE,AUC,CE) %>% as.data.frame %>% t %>% xtable(digits=4) %>% print(include.rownames=FALSE)







