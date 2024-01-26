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
library(gridExtra)

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
coords = coords #[-502,]
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
longlatdat = coords@coords %>% cbind(resp$x) %>% as_tibble %>% rename(resp = V3)
longlatdat$resp = as.factor(longlatdat$resp)
levels(longlatdat$resp) = c("Absence","Presence")

##-----------------
## Explore data
##-----------------
combinePointValue %>% 
  head

#ggpairs(combinePointValue)

#Some plots

#pdf('temp_data.pdf')

par(mfrow = c(2,2), las = 1)
plot(rmin, main = "Min Temp.", xlim=c(-101,-65), xlab = "Longitude", ylab = "Latitude")
points(coords$Longitude,coords$Latitude,col=combinePointValue$resp, cex = 0.2, pch = 19)
legend(-77,32,legend = c("Absence","Presence"),col = 0:1,pch=19, bg = "green")
plot(rmax, main = "Max Temp.", xlim=c(-101,-65), xlab = "Longitude", ylab = "Latitude")
points(coords$Longitude,coords$Latitude,col=combinePointValue$resp, cex = 0.2, pch = 19)
legend(-77,32,legend = c("Absence","Presence"),col = 0:1,pch=19, bg = "green")

#rm(rmin, rmax, r, rasValue)

rmin_pi = raster("us_tmin_pi.txt")/20 #%>% sqrt
rmax_pi = raster("us_tmax_pi.txt")/20 #%>% sqrt
PI_raster = stack(rmin_pi,rmax_pi)
PI_obs = raster::extract(PI_raster,coords)
PI_obs[882,1] = mean(PI_obs[,1], na.rm = TRUE) #881 if outlier excluded, 882 if outlier included
plot(rmin_pi, main = "Min Temp. Error", xlim=c(-101,-65), xlab = "Longitude", ylab = "Latitude")
plot(rmax_pi, main = "Max Temp. Error", xlim=c(-101,-65), xlab = "Longitude", ylab = "Latitude")

gg1 = rasterToPoints(rmin) %>% as_tibble %>%
  rename(Longitude = x, Latitude = y) %>%
  #left_join(longlatdat) %>%
  filter(Longitude > -101) %>%
  ggplot(aes(x=Longitude,y=Latitude,z=us_tmin)) +
  geom_contour_filled(bins = 9) + #breaks = seq(0,0.45,by=0.05)
  scale_fill_brewer(palette = "Blues") +
  #facet_wrap(. ~ names, nrow = 2) +
  theme_bw() +
  #labs(fill = "Value") +
  coord_cartesian(xlim=c(-100,-65)) +
  theme(
    #legend.position="right", 
    #axis.title=element_text(size=25), 
    #strip.text = element_text(size=12),
    #axis.text=element_text(size=25),aspect.ratio = 1,
    #axis.ticks.length=unit(.25, "cm"),
    #legend.key.size = unit(1.5, 'cm'), 
    #legend.key.height = unit(1.5, 'cm'), 
    #legend.key.width = unit(1.5, 'cm'),
    legend.title = element_blank(), 
    #legend.text = element_text(size=15),
    plot.title = element_text(hjust = 0.5)
  ) +
    ggtitle('Min Temp.') +
  geom_point(data = longlatdat, 
             aes(x = Longitude, y = Latitude, color = resp, z = NULL), size= 0.4)+
  scale_color_manual(values = c("Absence" = "peachpuff","Presence"="red"))


gg2 = rasterToPoints(rmax) %>% as_tibble %>%
  rename(Longitude = x, Latitude = y) %>% filter(Longitude > -101) %>%
  ggplot(aes(x=Longitude,y=Latitude,z=us_tmax)) +
  geom_contour_filled(bins = 9) + #breaks = seq(0,0.45,by=0.05)
  scale_fill_brewer(palette = "Blues") +
  #facet_wrap(. ~ names, nrow = 2) +
  theme_bw() +
  #labs(fill = "Value") +
  coord_cartesian(xlim=c(-100,-65)) +
  theme(
    legend.position="right", 
    #axis.title=element_text(size=25), 
    #strip.text = element_text(size=12),
    #axis.text=element_text(size=25),aspect.ratio = 1,
    #axis.ticks.length=unit(.25, "cm"),
    #legend.key.size = unit(1.5, 'cm'), 
    #legend.key.height = unit(1.5, 'cm'), 
    #legend.key.width = unit(1.5, 'cm'),
    legend.title = element_blank(), 
    #legend.text = element_text(size=15),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle('Max Temp.') +
  geom_point(data = longlatdat,
             aes(x = Longitude, y = Latitude, color = resp, z = NULL), size= 0.4)+
  scale_color_manual(values = c("Absence" = "peachpuff","Presence"="red"))


gg3 = rasterToPoints(rmin_pi) %>% as_tibble %>%
  rename(Longitude = x, Latitude = y) %>% filter(Longitude > -101) %>%
  ggplot(aes(x=Longitude,y=Latitude,z=us_tmin_pi)) +
  geom_contour_filled(binwidth = 0.25) + #breaks = seq(0,0.45,by=0.05)
  #scale_fill_brewer(palette = "YlGn") +
  #facet_wrap(. ~ names, nrow = 2) +
  theme_bw() +
  #labs(fill = "Value") +
  coord_cartesian(xlim=c(-100,-65)) +
  theme(
    legend.position="right", 
    #axis.title=element_text(size=25), 
    #strip.text = element_text(size=12),
    #axis.text=element_text(size=25),aspect.ratio = 1,
    #axis.ticks.length=unit(.25, "cm"),
    #legend.key.size = unit(1.5, 'cm'), 
    #legend.key.height = unit(1.5, 'cm'), 
    #legend.key.width = unit(1.5, 'cm'),
    legend.title = element_blank(), 
    #legend.text = element_text(size=15),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle('Min Temp. Error') 


gg4 = rasterToPoints(rmax_pi) %>% as_tibble %>%
  rename(Longitude = x, Latitude = y) %>% filter(Longitude > -101) %>%
  ggplot(aes(x=Longitude,y=Latitude,z=us_tmax_pi)) +
  geom_contour_filled(binwidth = 0.15) + #breaks = seq(0,0.45,by=0.05)
  #scale_fill_brewer(palette = "YlGn") +
  #facet_wrap(. ~ names, nrow = 2) +
  theme_bw() +
  #labs(fill = "Value") +
  coord_cartesian(xlim=c(-100,-65)) +
  theme(
    legend.position="right", 
    #axis.title=element_text(size=25), 
    #strip.text = element_text(size=12),
    #axis.text=element_text(size=25),aspect.ratio = 1,
    #axis.ticks.length=unit(.25, "cm"),
    #legend.key.size = unit(1, 'cm'), 
    #legend.key.height = unit(1.5, 'cm'), 
    #legend.key.width = unit(1.5, 'cm'),
    legend.title = element_blank(), 
    #legend.text = element_text(size=15),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle('Max Temp. Error') 

grid.arrange(gg1,gg2,gg3,gg4,ncol = 2)

ggsave("temp_data.pdf")


# temperature.dat = rmin_points %>% inner_join(rmax_points) %>% 
#   inner_join(rminpi_points) %>% inner_join(rmaxpi_points) %>%
#   pivot_longer(us_tmin:us_tmax_pi,names_to = 'names',values_to = 'values')



#dev.off()

#############-------------- Fit models and compare ------------#############

n = length(combinePointValue$resp)

sizes = rep(1,n)

xbf_num = 200
xsqbf_num = 200
x2bf_num = 150
spbf_num = 15

x_basisfunctions <- mrts(UTM@coords, k = xbf_num, x = UTM@coords) %>%
  as.matrix

x2_basisfunctions <- mrts(UTM@coords, k = x2bf_num, x = UTM@coords) %>%
  as.matrix

x3_basisfunctions <- mrts(UTM@coords, k = xsqbf_num, x = UTM@coords) %>%
  as.matrix

sp_basisfunctions <- mrts(UTM@coords, k = spbf_num, x = UTM@coords) %>%
  as.matrix
sp_basisfunctions <- sp_basisfunctions[,-1] #no intercept


#Fit the model for W via autoFRK()

W = combinePointValue$us_tmin
W2 = combinePointValue$us_tmax
W3 = combinePointValue$us_tmin^2

fit_Wmodel <- autoFRK(Data = matrix(W,ncol=1), loc = UTM@coords, G = x_basisfunctions)
tmin_pred <- predict(fit_Wmodel)$pred

fit_W2model <- autoFRK(Data = matrix(W2,ncol=1), loc = UTM@coords, G = x2_basisfunctions)
tmax_pred <- predict(fit_W2model)$pred

fit_W3model <- autoFRK(Data = matrix(W3,ncol=1), loc = UTM@coords, G = x3_basisfunctions)
tmin2_pred <- predict(fit_W3model)$pred

combinePointValue$tmin_pred = tmin_pred
combinePointValue$tmax_pred = tmax_pred

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

########----------Try choose number of basis functions----------########

# testMSE = rep(0,100)
# train = sample(n,800)
# 
# for (num_knots in 5:104) {
#   train_basisfunctions <- mrts(coords[train,], k = num_knots, x = coords[train,]) %>% as.matrix
#   train_model <- autoFRK(Data = matrix(combinePointValue$resp[train],ncol=1), 
#                         loc = coords[train,], G = train_basisfunctions)
#   train_alpha <- train_model$w
#   test_basisfunctions <- mrts(coords[train,], k = num_knots, x = coords[-train,]) %>% as.matrix
#   test_preds = test_basisfunctions%*%train_alpha
#   testMSE[num_knots-4] = sum((test_preds - combinePointValue$resp[-train])^2)/249
# }
# 
# testMSE
# plot(testMSE)




# combinePointValue %>% 
#   ggplot(aes(Latitude,us_tmin)) +
#   geom_point() +
#   theme_bw() +
#   labs(title = "Scatterplot", x = "Latitude" , y = 'Min Temp') +
#   scale_fill_brewer(palette = 1) +
#   geom_line(aes(y = tmin_pred), col = 2)
# 
# combinePointValue %>% 
#   ggplot(aes(Latitude,us_tmax)) +
#   geom_point() +
#   theme_bw() +
#   labs(title = "Scatterplot", x = "Latitude" , y = 'max Temp') +
#   scale_fill_brewer(palette = 1) +
#   geom_line(aes(y = tmax_pred), col = 2)
# 
# combinePointValue %>% 
#   ggplot(aes(Longitude,us_tmin)) +
#   geom_point() +
#   theme_bw() +
#   labs(title = "Scatterplot", x = "Longitude" , y = 'Min Temp') +
#   scale_fill_brewer(palette = 1) +
#   geom_line(aes(y = tmin_pred), col = 2)
# 
# combinePointValue %>% 
#   ggplot(aes(Longitude,us_tmax)) +
#   geom_point() +
#   theme_bw() +
#   labs(title = "Scatterplot", x = "Longitude" , y = 'max Temp') +
#   scale_fill_brewer(palette = 1) +
#   geom_line(aes(y = tmax_pred), col = 2)
# 
# 
# fig <- plot_ly(combinePointValue, x = ~Latitude, y = ~Longitude, z = ~us_tmin)
# fig <- fig %>% add_markers()
# fig <- fig %>% add_trace(z = ~tmin_pred, name = 'pred', mode = 'markers')
# fig
# 
# fig1 <- plot_ly(combinePointValue, x = ~Latitude, y = ~Longitude, z = ~us_tmax)
# fig1 <- fig1 %>% add_markers()
# fig1 <- fig1 %>% add_trace(z = ~tmax_pred, name = 'pred', mode = 'markers')
# fig1



#Fit using different methods/models

#From MEE paper

# mod1 = glm(resp ~ us_tmin + I(us_tmin^2) , family = binomial, data =  combinePointValue)
#
# mod2 = glm(resp ~ us_tmin + I(us_tmin^2) + us_tmax + I(us_tmax^2) , family = binomial, data =  combinePointValue)
#
# mod3 = glm(combinePointValue$resp ~ tmin_pred + I(tmin_pred^2) , family = binomial)
#
# mod4 = glm(combinePointValue$resp ~ tmin_pred + I(tmin_pred^2) + tmax_pred + I(tmax_pred^2) , family = binomial)
#
# mod5 = fitme(resp ~ us_tmin + I(us_tmin^2) + Matern(1|Longitude + Latitude),
#              data = combinePointValue, family = "binomial")
# 
#
# mod1_coef = mod1$coefficients
# pred_raster = rmin
# pred_raster@data@values = binomial()$linkinv(mod1$coefficients[1] +
#                                                mod1$coefficients[2]*rmin@data@values +
#                                                mod1$coefficients[3]*rmin@data@values^2)
#
# my.palette <- brewer.pal(n = 20, name = "YlGn")
# plot(pred_raster, main = "Predicted Probabilities", xlim=c(-101,-65), xlab = "Longitude",
#      ylab = "Latitude", col = my.palette)
#


# US_raster <- rasterToPoints(r) %>% as_tibble %>% select(x,y)
# 
# coordinates(US_raster) <- c("x" , "y")
# proj4string(US_raster) <- CRS("+proj=longlat +datum=WGS84")
# UTM_raster <- spTransform(US_raster, CRS("+proj=utm +zone=16 ellps=WGS84"))
# 
# x_preds <- as.matrix(mrts(UTM@coords, k = 100, x = UTM_raster@coords))%*%alphax_hat
# US_raster = US_raster %>% as_tibble %>%
#   mutate(probs = binomial()$linkinv(mod3$coefficients[1] +
#           mod3$coefficients[2]*x_preds +
#           mod3$coefficients[3]*x_preds^2))
# 
# probs_raster = rasterFromXYZ(US_raster)
# my.palette <- brewer.pal(n = 20, name = "YlGn")
# plot(probs_raster, main = "Predicted Probabilities", xlim=c(-101,-65), xlab = "Longitude",
#      ylab = "Latitude", col = my.palette)

#fitme model

# US_raster <- rasterToPoints(rmin) %>% as_tibble
# coordinates(US_raster) <- c("x" , "y")
# proj4string(US_raster) <- CRS("+proj=longlat +datum=WGS84")
# UTM_raster <- spTransform(US_raster, CRS("+proj=utm +zone=16 ellps=WGS84")) %>% as_tibble %>%
#   rename(Longitude = x, Latitude = y)
# US_raster = US_raster %>% as_tibble %>%
#   mutate(probs = predict(mod5,newdata = UTM_raster)) %>% select(x,y,probs)
# 
# probs_raster = rasterFromXYZ(US_raster)
# my.palette <- brewer.pal(n = 20, name = "YlGn")
# plot(probs_raster, main = "Predicted Probabilities", xlim=c(-101,-65), xlab = "Longitude",
#      ylab = "Latitude", col = my.palette)

#Deprecated
# US_raster <- rasterToPoints(rmin) %>% as_tibble
# raster_coords = US_raster %>% select(x,y)
# full_basisfunctions <- mrts(raster_coords, k = 100, x = raster_coords) %>% as.matrix
# full_model <- autoFRK(Data = matrix(US_raster[,3]%>%as.matrix,ncol=1), loc = raster_coords, G = full_basisfunctions)
# full_pred <- predict(full_model)$pred
#
# US_raster = US_raster %>% select(x,y) %>% mutate(probs =
#                                    binomial()$linkinv(mod3$coefficients[1] +
#                                                         mod3$coefficients[2]*full_pred +
#                                                         mod3$coefficients[3]*full_pred^2))




#Models with both min and max temp

# naive_mod = glm(resp ~ us_tmin + us_tmax + sp_basisfunctions , family = binomial, data =  combinePointValue)
# 
# no.spatial_mod = glm(resp ~ tmin_pred + tmax_pred , family = binomial, data = combinePointValue)
# 
# cholSigma = diag(log(abs(naive_mod$coefficients[-(1:3)])))[,1:1]
# cholSigma_vec = cholSigma[lower.tri(cholSigma, diag = TRUE)]
# 
# inputdata=list(sp_basis = sp_basisfunctions , y = as.matrix(combinePointValue$resp), 
#                predX = tmin_pred, predX2 = tmax_pred, true_rank = 1, 
#                identity_matrix = diag(nrow=ncol(sp_basisfunctions)), sizes = sizes) 
# 
# parameterlist=list(beta = naive_mod$coefficients[1:3] , cholSigma = cholSigma_vec , 
#                    alpha_rho = naive_mod$coefficients[-(1:3)])
# 
# 
# obj <- MakeADFun(data=inputdata,DLL = "multX_bin",
#                  parameters = parameterlist, random = c("alpha_rho"), silent=T)
# 
# 
# # fit <- optimr(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS", hessian = TRUE,
# #                       control=list(maxit=500, trace=0))
# 
# fit <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr,
#               control=list(trace = 0, iter.max = 1000, eval.max = 1000))
# 
# 
# ###-------Get SEs--------###
# 
# #Examine outputted SEs from TMB; *should* be correct when model is correct and true X is used
# 
# report = sdreport(obj,getJointPrecision = TRUE, getReportCovariance = FALSE, skip.delta.method = TRUE)
# report_summary=summary(report)
# params = report_summary %>% 
#   as.data.frame %>% 
#   filter(rownames(report_summary)%in%c("alpha_rho","beta")) %>%
#   dplyr::select("Estimate") %>% .$Estimate
# 
# sigmarho_hat = as.list(report,"Est",TRUE)$Sigma
# 
# #Non-parametric estimate of variance of score
# 
# 
# 
# varscore_x = 0
# varscore_x2 = 0
# varscore_y = 0
# 
# eta = params[1] + params[2]*tmin_pred + params[3]*tmax_pred + sp_basisfunctions%*%params[-(1:3)] 
# iterative_weights = sizes*binomial()$var(binomial()$linkinv(eta)) %>% as.vector %>% diag.spam
# 
# for (j in 1:n) {
#   varscore_x = varscore_x + tcrossprod(x_basisfunctions[j,]*(W[j] - tmin_pred[j])/phi_hat ) 
#   
#   varscore_x2 = varscore_x2 + tcrossprod(x2_basisfunctions[j,]*(W2[j] - tmax_pred[j])/phi_hat2 )
#   
#   score_y0 = combinePointValue$resp[j] - sizes[j]*binomial()$linkinv(eta[j]) 
#   
#   score_y1 = tmin_pred[j]*(combinePointValue$resp[j] - sizes[j]*binomial()$linkinv(eta[j])) 
#   
#   score_y2 = tmax_pred[j]*(combinePointValue$resp[j] - sizes[j]*binomial()$linkinv(eta[j]))
#   
#   score_y_rho = sp_basisfunctions[j,]*(combinePointValue$resp[j] - sizes[j]*binomial()$linkinv(eta[j])) 
#   score_y = c(score_y0,score_y1,score_y2,score_y_rho)
#   varscore_y = varscore_y + tcrossprod(score_y)
# }
# 
# 
# 
# varscore_y[-(1:3),-(1:3)] = varscore_y[-(1:3),-(1:3)] + ginv(sigmarho_hat,1*(.Machine$double.eps)^(1/4))
# varscore_x = varscore_x + ginv(sigmax_hat)
# varscore_x2 = varscore_x2 + ginv(sigmax2_hat)
# 
# #observed negative Hessians (joint)
# 
# Xstar = cbind(1, tmin_pred, tmax_pred, sp_basisfunctions)
# Y_hess = crossprod(Xstar,iterative_weights)%*%Xstar + 
#   bdiag.spam(matrix(0,nrow=3,ncol=3),ginv(sigmarho_hat,1*(.Machine$double.eps)^(1/4)))
# 
# X_hess = crossprod(x_basisfunctions)/phi_hat + ginv(sigmax_hat) #X part/block
# X2_hess = crossprod(x2_basisfunctions)/phi_hat2 + ginv(sigmax2_hat)
# 
# XY_hess = cbind(fit$par[2]*crossprod(Xstar,iterative_weights)%*%x_basisfunctions,
#                 fit$par[3]*crossprod(Xstar,iterative_weights)%*%x2_basisfunctions)
# 
# jointPrecision = bdiag(Y_hess, X_hess, X2_hess)
# jointPrecision[1:(3+ncol(sp_basisfunctions)),(4+ncol(sp_basisfunctions)):ncol(jointPrecision)] =
#   XY_hess
# 
# 
# varscore = bdiag.spam(varscore_y,varscore_x, varscore_x2)
# invnegHess = as.matrix(solve(jointPrecision))
# varcovar = tcrossprod(invnegHess%*%varscore,invnegHess)
# c(varcovar[1,1],varcovar[2,2],varcovar[3,3])
# 
# #compare the different methods
# 
# est_table = rbind(fit$par[1:3],naive_mod$coefficients[1:3],no.spatial_mod$coefficients)
# 
# summ_naive = summary(naive_mod)
# summ_no.spatial = summary(no.spatial_mod)
# 
# var_table = rbind(diag(varcovar)[1:3],diag(summ_naive$cov.scaled)[1:3],diag(summ_no.spatial$cov.scaled))
# 
# est_table = est_table %>% as_tibble %>% mutate(method = c("Correct","Naive","No Spatial"))
# var_table = var_table %>% as_tibble %>% mutate(method = c("Correct","Naive","No Spatial"))
# colnames(est_table) <- colnames(var_table) <- c("Intercept","tmin","tmax","Method")
# 
# est_table
# var_table
# 
# est_table = est_table %>% pivot_longer(Intercept:tmax, names_to = 'parameter')
# var_table = var_table %>% pivot_longer(Intercept:tmax, names_to = 'parameter')
# latex_table = inner_join(est_table,var_table,c('Method','parameter'))
# xtable(latex_table)

###################-------------Quadratic models with min temp only------------####################

#Naive


naive_mod = glm(resp ~ us_tmin + I(us_tmin^2) + sp_basisfunctions , family = binomial, data =  combinePointValue)

no.spatial_mod = glm(resp ~ tmin_pred + I(tmin_pred^2) , family = binomial, data = combinePointValue)

cholSigma = diag(log(abs(naive_mod$coefficients[-(1:3)])))[,1:1]
cholSigma_vec = cholSigma[lower.tri(cholSigma, diag = TRUE)]

inputdata=list(sp_basis = sp_basisfunctions , y = as.matrix(combinePointValue$resp), 
               predX = combinePointValue$us_tmin, predX2 = combinePointValue$us_tmin^2, true_rank = 1, 
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

inputdata=list(sp_basis = sp_basisfunctions , y = as.matrix(combinePointValue$resp), 
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
  
  score_y0 = combinePointValue$resp[j] - sizes[j]*binomial()$linkinv(eta[j]) 
  
  score_y1 = tmin_pred[j]*(combinePointValue$resp[j] - sizes[j]*binomial()$linkinv(eta[j])) 
  
  score_y2 = tmin_pred[j]^2*(combinePointValue$resp[j] - sizes[j]*binomial()$linkinv(eta[j]))
  
  score_y_rho = sp_basisfunctions[j,]*(combinePointValue$resp[j] - sizes[j]*binomial()$linkinv(eta[j])) 
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


inputdata=list(sp_basis = sp_basisfunctions , y = as.matrix(combinePointValue$resp), 
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
  
  score_y0 = combinePointValue$resp[j] - sizes[j]*binomial()$linkinv(eta[j]) 
  
  score_y1 = tmin_pred[j]*(combinePointValue$resp[j] - sizes[j]*binomial()$linkinv(eta[j])) 
  
  score_y2 = tmin2_pred[j]*(combinePointValue$resp[j] - sizes[j]*binomial()$linkinv(eta[j]))
  
  score_y_rho = sp_basisfunctions[j,]*(combinePointValue$resp[j] - sizes[j]*binomial()$linkinv(eta[j])) 
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

uGAM = gam(resp ~ us_tmin + I(us_tmin^2) + s(Longitude,Latitude) , family = binomial, data =  combinePointValue)

set.seed(2464)
simexGAM <- simex(uGAM, SIMEXvariable = c("us_tmin"),
                    measurement.error = PI_obs[,'us_tmin_pi'], asymptotic = FALSE)

cGAM = gam(resp ~ tmin_pred + I(tmin_pred^2) + s(Longitude,Latitude) , family = binomial, data =  combinePointValue)

uGAMsumm = summary(uGAM)
simexGAM_summ = summary(simexGAM)
cGAMsumm = summary(cGAM)

dnaive = glm(resp ~ us_tmin + I(us_tmin^2) , family = binomial, data =  combinePointValue)

set.seed(2464)
simexNS <- simex(dnaive, SIMEXvariable = c("us_tmin"),
                    measurement.error = PI_obs[,'us_tmin_pi'], asymptotic = FALSE)

simexNS_summ = summary(simexNS)
summ_dnaive = summary(dnaive)
summ_no.spatial = summary(no.spatial_mod)

#compare the different methods


est_table = rbind(fit2$par[1:3],fit3$par[1:3],fit1$par[1:3],dnaive$coefficients,no.spatial_mod$coefficients,
 simexNS$coefficients,uGAM$coefficients[1:3],cGAM$coefficients[1:3],simexGAM$coefficients[1:3])

var_table = rbind(together_vars,separate_vars,report_summary1[1:3,2]^2,diag(summ_dnaive$cov.scaled)[1:3],
                  diag(summ_no.spatial$cov.scaled),simexNS_summ$coefficients$jackknife[,2]^2,
                  diag(uGAMsumm$cov.scaled)[1:3],diag(cGAMsumm$cov.scaled)[1:3],simexGAM_summ$coefficients$jackknife[1:3,2]^2)

est_table = est_table %>% as_tibble %>% mutate(method = 
 c("dFRK","dFRK2","Naive",'dNaive',"NS",'SNS','uGAM','cGAM','sGAM'))
var_table = var_table %>% as_tibble %>% mutate(method =  
 c("dFRK","dFRK2","Naive",'dNaive',"NS",'SNS','uGAM','cGAM','sGAM'))

colnames(est_table) <- colnames(var_table) <- c("Intercept","tmin","tminsq","Method")

est_table
var_table

est_table = est_table %>% pivot_longer(Intercept:tminsq, names_to = 'Parameter', values_to = 'estimate')
var_table = var_table %>% pivot_longer(Intercept:tminsq, names_to = 'Parameter', values_to = 'SE')
var_table$SE = sqrt(var_table$SE)
latex_table = inner_join(est_table,var_table,c('Method','Parameter'))
latex_table
xtable(latex_table)

pivot_wider(est_table,names_from = Parameter,values_from = estimate) %>% 
  xtable(digits = c(0,0,3,3,4)) %>% print(include.rownames=FALSE)
pivot_wider(var_table,names_from = Parameter,values_from = SE) %>% 
  xtable(digits = c(0,0,3,4,4)) %>% print(include.rownames=FALSE)


###############----------------predicted probability plots-------------#####################

#Naive - uncorrected spatial

#pdf("predicted_probs.pdf",width=21,height = 21)

#par(mfrow=c(3,3))

US_raster <- rasterToPoints(rmin) %>% as_tibble
coordinates(US_raster) <- c("x" , "y")
proj4string(US_raster) <- CRS("+proj=longlat +datum=WGS84")
UTM_raster <- spTransform(US_raster, CRS("+proj=utm +zone=16 ellps=WGS84"))
rho_preds <- as.matrix(mrts(UTM@coords, k = spbf_num, x = UTM_raster@coords))[,-1]%*%params1[-(1:3)]
US_raster = US_raster %>% as_tibble %>% mutate(Anaive = binomial()$linkinv(params1[1] +
                                                        params1[2]*us_tmin +
                                                        params1[3]*us_tmin^2 + rho_preds))
all_probs = US_raster %>% select(x,y,Anaive) %>% rename(Longitude = x, Latitude = y)

# probs_raster = rasterFromXYZ(US_raster)
# my.palette <- brewer.pal(n = 20, name = "YlGn")
# plot(probs_raster, main = "Naive", xlim=c(-101,-65), xlab = "Longitude",
#      ylab = "Latitude", col = my.palette)



#Corrected

US_raster <- rasterToPoints(rmin) %>% as_tibble %>% select(x,y)
coordinates(US_raster) <- c("x" , "y")
proj4string(US_raster) <- CRS("+proj=longlat +datum=WGS84")
UTM_raster <- spTransform(US_raster, CRS("+proj=utm +zone=16 ellps=WGS84"))
x_preds <- as.matrix(mrts(UTM@coords, k = xbf_num, x = UTM_raster@coords ))%*%alphax_hat
rho_preds <- as.matrix(mrts(UTM@coords, k = spbf_num, x = UTM_raster@coords ))[,-1]%*%params2[-(1:3)]
US_raster = US_raster %>% as_tibble %>% mutate(probs =  binomial()$linkinv(params2[1] +
                                                                             params2[2]*x_preds +
                                                                             params2[3]*x_preds^2 + rho_preds))
all_probs$Bcorrect = US_raster$probs

# probs_raster = rasterFromXYZ(US_raster)
# my.palette <- brewer.pal(n = 20, name = "YlGn")
# plot(probs_raster, main = "dFRK", xlim=c(-101,-65), xlab = "Longitude",
#      ylab = "Latitude", col = my.palette)


#Corrected V2

US_raster <- rasterToPoints(rmin) %>% as_tibble %>% select(x,y)
coordinates(US_raster) <- c("x" , "y")
proj4string(US_raster) <- CRS("+proj=longlat +datum=WGS84")
UTM_raster <- spTransform(US_raster, CRS("+proj=utm +zone=16 ellps=WGS84"))
#x_preds <- as.matrix(mrts(UTM@coords, k = xbf_num, x = UTM_raster@coords ))%*%alphax_hat
x2_preds <- as.matrix(mrts(UTM@coords, k = xsqbf_num, x = UTM_raster@coords ))%*%alphax3_hat
rho_preds <- as.matrix(mrts(UTM@coords, k = spbf_num, x = UTM_raster@coords ))[,-1]%*%params3[-(1:3)]
US_raster = US_raster %>% as_tibble %>% mutate(probs =  binomial()$linkinv(params3[1] +
                                                        params3[2]*x_preds +
                                                        params3[3]*x2_preds + rho_preds))
all_probs$CcorrectV2 = US_raster$probs

# probs_raster = rasterFromXYZ(US_raster)
# my.palette <- brewer.pal(n = 20, name = "YlGn")
# plot(probs_raster, main = "dFRK2", xlim=c(-101,-65), xlab = "Longitude",
#      ylab = "Latitude", col = my.palette)


#dNaive


# pred_raster = rmin
# pred_raster@data@values = binomial()$linkinv(dnaive$coefficients[1] +
#                                                dnaive$coefficients[2]*rmin@data@values +
#                                                dnaive$coefficients[3]*rmin@data@values^2)

# my.palette <- brewer.pal(n = 20, name = "YlGn")
# plot(pred_raster, main = "dNaive", xlim=c(-101,-65), xlab = "Longitude",
#      ylab = "Latitude", col = my.palette)


US_raster <- rasterToPoints(rmin) %>% as_tibble %>%
mutate(probs = binomial()$linkinv(dnaive$coefficients[1] +
                                               dnaive$coefficients[2]*us_tmin +
                                               dnaive$coefficients[3]*us_tmin^2))

all_probs$Fdnaive = US_raster$probs

#No spatial


US_raster <- rasterToPoints(rmin) %>% as_tibble %>% select(x,y)

coordinates(US_raster) <- c("x" , "y")
proj4string(US_raster) <- CRS("+proj=longlat +datum=WGS84")
UTM_raster <- spTransform(US_raster, CRS("+proj=utm +zone=16 ellps=WGS84"))

#x_preds <- as.matrix(mrts(UTM@coords, k = xbf_num, x = UTM_raster@coords))%*%alphax_hat
US_raster = US_raster %>% as_tibble %>%
  mutate(probs = binomial()$linkinv(no.spatial_mod$coefficients[1] +
                                      no.spatial_mod$coefficients[2]*x_preds +
                                      no.spatial_mod$coefficients[3]*x_preds^2))

all_probs$GNS = US_raster$probs

# probs_raster = rasterFromXYZ(US_raster)
# my.palette <- brewer.pal(n = 20, name = "YlGn")
# plot(probs_raster, main = "NS", xlim=c(-101,-65), xlab = "Longitude",
#      ylab = "Latitude", col = my.palette)


#SimexNS

# pred_raster = rmin
# pred_raster@data@values = binomial()$linkinv(simexNS$coefficients[1] +
#                                                simexNS$coefficients[2]*rmin@data@values +
#                                                simexNS$coefficients[3]*rmin@data@values^2)
# 
# my.palette <- brewer.pal(n = 20, name = "YlGn")
# plot(pred_raster, main = "SNS", xlim=c(-101,-65), xlab = "Longitude",
#      ylab = "Latitude", col = my.palette)

US_raster <- rasterToPoints(rmin) %>% as_tibble %>%
  mutate(probs = binomial()$linkinv(simexNS$coefficients[1] +
                                      simexNS$coefficients[2]*us_tmin +
                                      simexNS$coefficients[3]*us_tmin^2))

all_probs$HSNS = US_raster$probs


#uGAM


US_raster <- rasterToPoints(rmin) %>% as_tibble
coordinates(US_raster) <- c("x" , "y")
proj4string(US_raster) <- CRS("+proj=longlat +datum=WGS84")
UTM_raster <- spTransform(US_raster, CRS("+proj=utm +zone=16 ellps=WGS84"))
UTM_raster = UTM_raster %>% as_tibble %>% rename(Longitude = x, Latitude = y)
US_raster = US_raster %>% as_tibble %>% select(x,y) %>% 
  mutate(probs = binomial()$linkinv(predict(uGAM,newdata = UTM_raster)))

all_probs$DuGAM = US_raster$probs


# probs_raster = rasterFromXYZ(US_raster)
# my.palette <- brewer.pal(n = 20, name = "YlGn")
# plot(probs_raster, main = "uGAM", xlim=c(-101,-65), xlab = "Longitude",
#      ylab = "Latitude", col = my.palette)


#cGAM


# US_raster <- rasterToPoints(rmin) %>% as_tibble %>% select(x,y)
# coordinates(US_raster) <- c("x" , "y")
# proj4string(US_raster) <- CRS("+proj=longlat +datum=WGS84")
# UTM_raster <- spTransform(US_raster, CRS("+proj=utm +zone=16 ellps=WGS84")) 
# x_preds <- as.matrix(mrts(UTM@coords, k = xbf_num, x = UTM_raster@coords ))%*%alphax_hat
# UTM_raster = UTM_raster %>% as_tibble %>% rename(Longitude = x, Latitude = y) %>% mutate(tmin_pred = x_preds)
# US_raster = US_raster %>% as_tibble %>% mutate(probs =  
#                                                  binomial()$linkinv(predict(cGAM,newdata = UTM_raster)))
# 
# probs_raster = rasterFromXYZ(US_raster)
# my.palette <- brewer.pal(n = 20, name = "YlGn")
# plot(probs_raster, main = "cGAM", xlim=c(-101,-65), xlab = "Longitude",
#      ylab = "Latitude", col = my.palette)


#SimexGAM


US_raster <- rasterToPoints(rmin) %>% as_tibble
coordinates(US_raster) <- c("x" , "y")
proj4string(US_raster) <- CRS("+proj=longlat +datum=WGS84")
UTM_raster <- spTransform(US_raster, CRS("+proj=utm +zone=16 ellps=WGS84"))
UTM_raster = UTM_raster %>% as_tibble %>% rename(Longitude = x, Latitude = y)
US_raster = US_raster %>% as_tibble %>% select(x,y) %>% mutate(probs =  
                binomial()$linkinv(predict(simexGAM,newdata = UTM_raster)))

# probs_raster = rasterFromXYZ(US_raster)
# my.palette <- brewer.pal(n = 20, name = "YlGn")
# plot(probs_raster, main = "sGAM", xlim=c(-101,-65), xlab = "Longitude",
#      ylab = "Latitude", col = my.palette)

all_probs$EsGAM = US_raster$probs

all_probs = all_probs %>% pivot_longer(Anaive:EsGAM,names_to = 'method',values_to = 'Probability')

all_probs$names = factor(all_probs$method,labels = c('FRK[naive]','dFRK','dFRK2','GAM[naive]','GAM[simex]','NS[naive]','NS[corr]','NS[simex]'))

my.palette <- brewer.pal(n = 20, name = "YlGn")


pdf("predicted_probs.pdf",width=28,height = 14)

#par(mar = c(4, 4, 0.1, 0.1)) 

ggplot(all_probs,aes(x=Longitude,y=Latitude,fill=Probability)) +
geom_raster() +
scale_fill_gradientn(colours = my.palette) +
facet_wrap(. ~ names, nrow = 2, labeller = label_parsed) +
theme_classic() +
coord_cartesian(xlim=c(-100,-65)) +
#scale_fill_distiller(type = 'seq', palette = 1) +
theme(legend.position="bottom", axis.title=element_text(size=25), 
      strip.text = element_text(size=25,face="bold"),
      axis.text=element_text(size=25),aspect.ratio = 1,
      axis.ticks.length=unit(.25, "cm")) +
theme(legend.key.size = unit(2, 'cm'), #change legend key size
      legend.key.height = unit(2, 'cm'), #change legend key height
      legend.key.width = unit(2, 'cm'), #change legend key width
      legend.title = element_text(size=25), #change legend title font size
      legend.text = element_text(size=25)) #change legend text font size

#ggsave("predicted_probs.pdf", width = 28, height = 14)

dev.off()

































