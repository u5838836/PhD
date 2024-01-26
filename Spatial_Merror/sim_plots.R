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
library(parallel)
library(foreach)
library(doParallel)
library(glm2)
library(xtable)
library(gridExtra)
library(colorspace)

setwd("D:/Study new/PhD/R code/Spatial_Merror")

true_coefs <- c(-0.1,1,1)
m <- 22500

set.seed(24)
xy_all <- data.frame(x = runif(m), y = runif(m))

X2_all <- RFsimulate(model = RMmatern(nu = 1.5, var = 0.5, scale = 0.1), x = xy_all$x, y = xy_all$y, n = 1)@data[,1]

X_all <- RFsimulate(model = RMmatern(nu = 1, var = 0.5, scale = 0.1), x = xy_all$x, y = xy_all$y, n = 1)@data[,1]

s = 2*(2*xy_all$x - 1)
r = 2*(2*xy_all$y - 1)
spranef_all <- 2*exp(-(s^2+r^2)) 

Matern(d = 0.2, range = 0.1, smoothness = 1.5)
Matern(d = 0.2, range = 0.1, smoothness = 1)

#observed X's and rho's

plot_data = xy_all %>% rename(longitude = y, latitude = x) %>% 
  mutate(X1 = X_all, X2 = X2_all, spranef = spranef_all) %>% as_tibble

p1 = ggplot(plot_data, aes(x = latitude, y = longitude, colour = X1)) +
  geom_point() + #shape = 1, size = 1
  scale_colour_viridis_c() +
  labs(title = "Covariate one", x = bquote(s[i1]) , y = bquote(s[i2]), colour = bquote(X[1])) +
  theme_bw() + 
  theme(legend.position="right",plot.title = element_text(size = 20, hjust = 0.5),
        axis.text=element_text(size=15), legend.text = element_text(size=20),
        axis.title=element_text(size=20),legend.title=element_text(size=20))

p2 = ggplot(plot_data, aes(x = latitude, y = longitude, colour = X2)) +
  geom_point() + #shape = 1, size = 1
  scale_colour_viridis_c() +
  labs(title = "Covariate two", x = bquote(s[i1]) , y = bquote(s[i2]), colour = bquote(X[2])) +
  theme_bw() + 
  theme(legend.position="right",plot.title = element_text(size = 20, hjust = 0.5),
        axis.text=element_text(size=15), legend.text = element_text(size=20),
        axis.title=element_text(size=20),legend.title=element_text(size=20))

p3 = ggplot(plot_data, aes(x = latitude, y = longitude, colour = spranef)) +
  geom_point() + #shape = 1, size = 1
  scale_colour_viridis_c() +
  labs(title = "Residual spatial effect", x = bquote(s[i1]) , y = bquote(s[i2]), colour = bquote(H[rho])) +
  theme_bw() + 
  theme(legend.position="right",plot.title = element_text(size = 20, hjust = 0.5),
        axis.text=element_text(size=15), legend.text = element_text(size=20),
        axis.title=element_text(size=20),legend.title=element_text(size=20))

pdf("observed_all.pdf",width = 21)
grid.arrange(p1,p2,p3, ncol=3)
dev.off()


#Matern correlations and true rho field

distances = seq(0,sqrt(2),length.out = 1000)

matern_true = tibble(distances = distances) %>% 
  mutate(field1 = Matern(d = distances, range = 0.1, smoothness = 1.5), 
         field2 = Matern(d = distances, range = 0.1, smoothness = 1))

p4 = matern_true %>% pivot_longer(cols = field1:field2, names_to = "name", values_to = 'correlation') %>%
  arrange(name) %>% mutate(nu = ifelse(name == "field1",'1.5','1.0')) %>%
  ggplot(aes(x = distances, y = correlation, col = nu)) +
  geom_line() +
  labs(title = "Matern fields", x = "Euclidean distance" , y = "Correlation", col = expression(nu)) +
  theme_bw() + 
  theme(legend.position="right",plot.title = element_text(size = 20, hjust = 0.5),
        axis.text=element_text(size=10), legend.text = element_text(size=10),
        axis.title=element_text(size=15))

x = 0:1000/1000
y = 0:1000/1000
s = 2*(2*x - 1)
r = 2*(2*y - 1)
true_rho =  expand.grid(s,r) %>% as_tibble %>% 
  rename(s = Var1, r = Var2) %>% mutate(spranef = 2*exp(-(s^2+r^2)) )
true_rho = expand.grid(x,y) %>% cbind(true_rho) %>% rename(x = Var1, y = Var2)


p5 = ggplot(true_rho, aes(x = x, y = y, z = spranef)) +
  geom_raster(aes(fill = spranef),hjust=0.5,vjust=0.5,interpolate=TRUE) +
  labs(title = "True residual spatial effect", x = bquote(s[i1]) , y = bquote(s[i2]), fill = bquote(H[rho])) +
  theme_bw() + 
  theme(legend.position="right",plot.title = element_text(size = 20, hjust = 0.5),
        axis.text=element_text(size=10), legend.text = element_text(size=10),
        axis.title=element_text(size=15)) +
  scale_fill_viridis_c()

pdf("true_all.pdf",width = 14)
grid.arrange(p4,p5, ncol=2)
dev.off()

#############-------------Boxplots, coverage for sims---------------###############

setwd("D:/Study new/Sim_Results.spatial")
setwd("D:/Study new/Sim_Results.spatial/0.5")
setwd("D:/Study new/Sim_Results.spatial/size5")
setwd("D:/Study new/Sim_Results.spatial/size10")
setwd("D:/Study new/Sim_Results.spatial/BF200")

resp_type = "pois"
big_table_dummy = NULL
num_runs = 100
for (i in 1:num_runs) {
  load(paste0(i,"_",resp_type,"_bigtable.Rdata"))
  big_table_dummy = rbind(big_table_dummy,big_table)
}
big_table = big_table_dummy
big_table$type[big_table$type=="naive"] = "Naive"
big_table$type[big_table$type=="dnaive"] = "dNaive"
big_table$type[big_table$type=="correct"] = "dFRK"
big_table$type[big_table$type=="no_spatial"] = "NS"

big_table = big_table %>% group_by(type,name,ssize,merror) %>% 
  mutate(true_sd = sqrt(var(estimate))) %>% 
  mutate(ratio = mean(sqrt(evar))/true_sd) %>% ungroup %>%
  mutate(tratio = paste0(type,'\n(',round(ratio,digits = 2),')'))


full_table = big_table %>% group_by(type,name,ssize,merror) %>% 
  summarise(mean = mean(estimate, na.rm = TRUE), sd = sd(estimate, na.rm = TRUE),
            se = sqrt(mean(evar, na.rm = TRUE))) %>% arrange(ssize,merror)

xtable(full_table,digits=c(0,0,0,0,2,3,3,3))

big_table = big_table %>% filter(name != "intercept") %>%
  mutate(lower = estimate - 1.96*sqrt(evar), upper = estimate + 1.96*sqrt(evar)) %>%
  mutate(covered = ifelse(upper > 1 & lower < 1 , 1 , 0 ))

coverage_table = big_table %>% group_by(type,name,ssize,merror) %>% summarise(coverage = mean(covered))
xtable(coverage_table,digits=3)


maxratio1 = big_table %>% filter(name == "slope") %>% select(ratio) %>% max
maxratio2 = big_table %>% filter(name == "slope2") %>% select(ratio) %>% max

big_table$type = factor(big_table$type,levels = c("dFRK","NS","Naive","dNaive"))
levels(big_table$type)

pdf(paste0('slope1_',resp_type,'.pdf'))

p6 = big_table %>% filter(name == "slope") %>%
    ggplot(aes(x = type, y = estimate)) +
    geom_boxplot(aes(fill = ratio)) +
    facet_grid(merror ~ ssize) +
    theme_bw() +
    geom_hline(aes(yintercept = 1),linetype="dotted",colour = "black") +
    labs(#title = ifelse(resp_type == "pois","Poisson first slope","Binomial first slope"),
         x = 'Procedure' , y = 'Estimate', fill = 'Ratio') +
    #paletteer_c("ggthemes::Red-Blue Diverging", 30) +
    #scale_fill_brewer(type = 'div', palette = 1) +
    scale_fill_continuous_divergingx(palette = 'RdBu', mid = 1, limits = c(0,1.5)) +
    #scale_fill_distiller(palette = "Spectral") + 
    theme(legend.position="right",plot.title = element_text(size = 22, hjust = 0.5),
          axis.text=element_text(size=11, vjust = 0.5, hjust=0.5, angle = 90), 
          axis.title=element_text(size=14,face="bold"), panel.spacing.y = unit(1.5, "lines")) +
    scale_x_discrete(labels=c('Naive' = expression(FRK[Naive]), 'NS' = expression(NS[Corr]),
                              "dNaive" = expression(NS[Naive]))) +
    theme(axis.title=element_blank())

p6

dev.off()



pdf(paste0('slope2_',resp_type,'.pdf'))

p7 = big_table %>% filter(name == "slope2") %>%
  ggplot(aes(x = type, y = estimate)) +
  geom_boxplot(aes(fill = ratio)) +
  facet_grid(merror ~ ssize) +
  theme_bw() +
  geom_hline(aes(yintercept = 1),linetype="dotted",colour = "black") +
  labs(#title = ifelse(resp_type == "pois","Poisson second slope","Binomial second slope"), 
       x = 'Procedure' , y = 'Estimate', fill = 'Ratio') + #
  #paletteer_c("ggthemes::Red-Blue Diverging", 30) +
  #scale_fill_brewer(type = 'div', palette = 1) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 1, limits = c(0,1.5)) +
  #scale_fill_distiller(palette = "Spectral") + 
  theme(legend.position="right",plot.title = element_text(size = 22, hjust = 0.5),
        axis.text=element_text(size=11, vjust = 0.5, hjust=0.5, angle = 90), 
        axis.title=element_text(size=14,face="bold")) +
  scale_x_discrete(labels=c('Naive' = expression(FRK[Naive]), 'NS' = expression(NS[Corr]),
                            "dNaive" = expression(NS[Naive]))) +
  theme(axis.title=element_blank())

p7


dev.off()



########### For ECSSC #############

pdf('ECSSC.pdf')

ssize.labs <- c("n = 250", "n = 500", "n = 1000")
names(ssize.labs) <- c(250,500,1000)

p8 = big_table %>% filter(name == "slope" & merror == 0.1) %>%
  ggplot(aes(x = type, y = estimate)) +
  geom_boxplot(aes(fill = type)) +
  facet_wrap(. ~ ssize, labeller = labeller(ssize = ssize.labs)) +
  theme_bw() +
  geom_hline(aes(yintercept = 1),linetype="dotted",colour = "black") +
  labs(#title = ifelse(resp_type == "pois","Poisson first slope","Binomial first slope"),
    x = 'Procedure' , y = 'Estimate', fill = 'Ratio') +
  #paletteer_c("ggthemes::Red-Blue Diverging", 30) +
  #scale_fill_brewer(type = 'div', palette = 1) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 1, limits = c(0,1.5)) +
  #scale_fill_distiller(palette = "Spectral") + 
  theme(legend.position="none",plot.title = element_text(size = 22, hjust = 0.5),
        axis.text=element_text(size=11, vjust = 0.5, hjust=0.5, angle = 90), 
        axis.title=element_text(size=14,face="bold"), panel.spacing.y = unit(1.5, "lines")) +
  scale_x_discrete(labels=c('Naive' = expression(FRK[Naive]), 'NS' = expression(NS[Corr]),
                            "dNaive" = expression(NS[Naive]))) +
  #theme(axis.title=element_blank())+
  scale_fill_manual(breaks = type,
                    values = c("cyan", "yellow", "red", "green"))

p8

dev.off()

