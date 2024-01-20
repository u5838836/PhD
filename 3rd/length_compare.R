###########------------ conditional vs unconditional intervals -----------------##############
rm(list=ls())
library(dplyr)
library(ggplot2)
setwd("D:/Study new/PhD/R code/3rd/Figures")

# sigma_b = 1
# quantiles = seq(0.5,0.99,0.01)
# mixn = rep(0,length(quantiles))
# for (q in 1:length(quantiles)) {
#   SSS = rep(0,50000)
#   for (i in 1:50000) {
#     b_dot = rnorm(1,0,sigma_b)
#     SSS[i] = pnorm(qnorm(quantiles[q])*exp(0.125*sigma_b^2+b_dot/2))
#   }
#   mixn[q] = mean(SSS)
# }
#plot(quantiles,mixn,xlab = "cond_quantile", ylab = "mixN_coverage")
#abline(a = 0,b = 1)
#can roughly interpret as y-axis is expected conditional interval length, x-axis is mixN interval length


all_dat = tibble()
alldat1 = tibble()
allcondlengths = tibble()

for (sigma_b in c(0.1,1,5)) {
  #actual lengths
  
  beta = 1
  alphas = seq(0.01,0.99,0.01)
  conditional_halflengths = rep(0,length(alphas))
  unconditional_halflengths = rep(0,length(alphas))
  clensq = rep(0,length(alphas))
  for (j in 1:length(alphas)) {
    SSSS = rep(0,50000)
    SSS = rep(0,50000)
    for (i in 1:50000) {
      b_dot = rnorm(1,0,sigma_b)
      SSSS[i] = rnorm(1,0,exp((- beta - b_dot)/2))
      SSS[i] = qnorm(1-alphas[j]/2)*exp((-beta-b_dot)/2)
    }
    unconditional_halflengths[j] = quantile(SSSS,probs=1-alphas[j]/2)
    conditional_halflengths[j] = qnorm(1-alphas[j]/2)*exp(- beta/2 + 0.125*sigma_b^2)
    clensq[j] = mean(SSS^2)
  }
  #plot(unconditional_halflengths,conditional_halflengths,main="Fixed sigma, changing alpha")
  #abline(a = 0,b = 1)
  dat = tibble(e_cond_len = conditional_halflengths, uncond_len = unconditional_halflengths, alpha = alphas,
               logcond = log(conditional_halflengths), loguncond = log(unconditional_halflengths),
               naive_norm = qnorm(1-alphas/2)*exp(- beta/2 + sigma_b^2/4), clensq = clensq, 
               uncond_lensq = unconditional_halflengths^2, sigma_b = sigma_b )
  
  all_dat = rbind(all_dat,dat)
  
  #fixed sigma and alpha
  
  #sigma_b = 0.1
  alpha = 0.05
  SSSS = rep(0,50000)
  cond_lengths = rep(0,50000)
  for (i in 1:50000) {
    b_dot = rnorm(1,0,sigma_b)
    SSSS[i] = rnorm(1,0,exp((-beta-b_dot)/2))
    cond_lengths[i] = qnorm(1-alpha/2)*exp((-beta-b_dot)/2)
  }
  unconditional_length = quantile(SSSS,probs=1-alpha/2)
  econditional_length = qnorm(1-alpha/2)*exp(-beta/2 + 0.125*sigma_b^2)
  
  dat1 = tibble(uncond_len = unconditional_length, econd_len = econditional_length, sigma_b = sigma_b)
  alldat1 = rbind(alldat1,dat1)
  allcondlength = tibble(cond_lengths, sigma_b = sigma_b)
  allcondlengths = rbind(allcondlengths,allcondlength)
}

all_dat$sigma_b = as.factor(all_dat$sigma_b)
levels(all_dat$sigma_b)
levels(all_dat$sigma_b) = c(expression(sigma[b] == 0.1),
                            expression(sigma[b] == 1), expression(sigma[b] == 5))

all_dat %>% 
  ggplot(aes(x = uncond_len,y = e_cond_len)) +
  facet_wrap(~sigma_b, scales = 'free') +
  geom_point(aes(col = alpha)) +
  theme_bw() +
  labs(x = "Unconditonal Length" , y = "Expected Conditional Length") +
  geom_abline(slope = 1,intercept = 0,col="red")
ggsave(filename = "cond_uncond.pdf",height = 14/3,width = 14)


all_dat %>% 
  ggplot(aes(x = uncond_len,y = naive_norm)) +
  facet_wrap(~sigma_b, scales = 'free') +
  geom_point(aes(col = alpha)) +
  theme_bw() +
  labs(x = "Unconditional Length" , y = "Naive Normal Length") +
  geom_abline(slope = 1,intercept = 0,col="red")
ggsave(filename = "naive_uncond.pdf",height = 14/3,width = 14)


# dat %>% 
#   ggplot(aes(x = uncond_lensq,y = clensq)) +
#   geom_point(aes(col = alpha)) +
#   theme_bw() +
#   geom_abline(slope = 1,intercept = 0,col="red")


# dat %>% 
#   ggplot(aes(x = loguncond, y = logcond)) +
#   geom_point(aes(col = alpha)) +
#   theme_bw() +
#   geom_abline(slope = 1,intercept = 0,col="red")


allcondlengths = allcondlengths %>% left_join(alldat1,by="sigma_b")

allcondlengths %>% 
  ggplot(aes(x = cond_lengths)) +
  facet_wrap(~sigma_b, scales = 'free') +
  geom_histogram() +
  theme_bw() +
  labs(x = "Conditional Interval Lengths") +
  geom_vline(aes(xintercept = uncond_len), col='red') +
  geom_vline(aes(xintercept = econd_len), col="blue")
ggsave(filename = "cond_hist.pdf",height = 14/3,width = 14)



###########----------- Bernoulli ----------------###############

#actual (half)lengths

# sigma_b = 5
# alphas = seq(0.01,0.99,0.01)
# conditional_halflengths = rep(0,length(alphas))
# unconditional_halflengths = rep(0,length(alphas))
# clensq = rep(0,length(alphas))
# var_uncond = rep(0,length(alphas))
# for (j in 1:length(alphas)) {
#   SSSS = rep(0,50000)
#   SSS = rep(0,50000)
#   for (i in 1:50000) {
#     b_dot = rnorm(1,0,sigma_b)
#     p = exp(b_dot)/(1+exp(b_dot))
#     SSSS[i] = rnorm(1,0,1/sqrt(p*(1-p)))
#     SSS[i] = qnorm(1-alphas[j]/2)/sqrt(p*(1-p))
#   }
#   unconditional_halflengths[j] = quantile(SSSS,probs=1-alphas[j]/2)
#   conditional_halflengths[j] = mean(SSS)
#   clensq[j] = mean(SSS^2)
#   var_uncond[j] = var(SSSS)
# }
# dat = tibble(e_cond_len = conditional_halflengths, uncond_len = unconditional_halflengths, alpha = alphas,
#              logcond = log(conditional_halflengths), loguncond = log(unconditional_halflengths),
#              naive_norm = qnorm(1-alphas/2)*sqrt(var_uncond[1]), clensq = clensq, uncond_lensq = unconditional_halflengths^2 )
# dat %>% 
#   ggplot(aes(x = uncond_len,y = e_cond_len)) +
#   geom_point(aes(col = alpha)) +
#   theme_bw() +
#   labs(x = "Unconditonal Length" , y = "Expected Conditional Length") +
#   geom_abline(slope = 1,intercept = 0,col="red")
# 
# dat %>% 
#   ggplot(aes(x = uncond_len,y = naive_norm)) +
#   geom_point(aes(col = alpha)) +
#   theme_bw() +
#   labs(x = "Unconditional Length" , y = "Naive Normal Length") +
#   geom_abline(slope = 1,intercept = 0,col="red")
# 
# #fixed sigma and alpha
# 
# sigma_b = 5
# alpha = 0.05
# SSSS = rep(0,50000)
# cond_lengths = rep(0,50000)
# for (i in 1:50000) {
#   b_dot = rnorm(1,0,sigma_b)
#   p = exp(b_dot)/(1+exp(b_dot))
#   SSSS[i] = rnorm(1,0,1/sqrt(p*(1-p)))
#   cond_lengths[i] = qnorm(1-alpha/2)/sqrt(p*(1-p))
# }
# unconditional_length = quantile(SSSS,probs=1-alpha/2)
# econditional_length = mean(cond_lengths)
# hist(cond_lengths,xlab='Conditonal Interval Lengths',probability = TRUE, main = '')
# abline(v = unconditional_length,col='red')
# abline(v = econditional_length, col = 'blue')
# 






