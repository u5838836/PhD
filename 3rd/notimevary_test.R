# rm(list = ls())
# library(tidyverse)
# library(lme4)
# library(mvtnorm)
# library(GGally)
# 
# num_individuals <- 10000
# num_timepoints <- 3
# dat <- data.frame(id = rep(1:num_individuals, each = num_timepoints),
#                   age_baseline = rep(rnorm(num_individuals, sd = 5), each = num_timepoints),
#                   gender = rep(sample(c("male","female"), num_individuals, replace = TRUE), each = num_timepoints)
# )
# dat
# 
# true_ranef <- rmvnorm(num_individuals, sigma = matrix(c(1,0.5,0.5,1), nrow = 2)) #cbind(rnorm(num_individuals, sd = 1), rnorm(num_individuals, sd = 0.5))
# true_eta <- 15 + 4*dat$age_baseline + rowSums(cbind(1,dat$age_baseline)*true_ranef[rep(1:num_individuals, each = num_timepoints),])
# dat$y <- rnorm(length(true_eta), mean = true_eta, sd = 1)
# 
# fit <- lmer(y ~ age_baseline + (age_baseline| id), data = dat, REML = FALSE) #' This is a misspecified model, and you can take out gender from the random slope if want to make things look comparable
# 
# summary(fit)
# 
# getME(fit, "Z") %>% head
# 
# compare_random_effects <- data.frame(est = ranef(fit)$id, true = true_ranef) 
# summary(compare_random_effects)
# 
# ggpairs(compare_random_effects, 
#         lower = "blank", 
#         upper = list(continuous = function(data, mapping) {
#           ggplot(data = data, mapping = mapping) +
#             geom_point() +
#             geom_abline(intercept = 0, slope = 1)}
#         ))


#time-varying test

rm(list = ls())
library(tidyverse)
library(lme4)
library(mvtnorm)
library(GGally)

num_individuals <- 20000
num_timepoints <- 3
dat <- data.frame(id = rep(1:num_individuals, each = num_timepoints),
                  age = rnorm(num_individuals*num_timepoints, sd = 5),
                  gender = rep(sample(c("male","female"), num_individuals, replace = TRUE), each = num_timepoints)
)
dat

true_ranef <- rmvnorm(num_individuals, sigma = matrix(c(1,0.5,0.5,1), nrow = 2)) #cbind(rnorm(num_individuals, sd = 1), rnorm(num_individuals, sd = 0.5))
true_eta <- 15 + 4*dat$age + rowSums(cbind(1,dat$age)*true_ranef[rep(1:num_individuals, each = num_timepoints),])
dat$y <- rnorm(length(true_eta), mean = true_eta, sd = 1)

fit <- lmer(y ~ age + (age| id), data = dat, REML = FALSE) #' This is a misspecified model, and you can take out gender from the random slope if want to make things look comparable

summary(fit)

getME(fit, "Z") %>% head

compare_random_effects <- data.frame(est = ranef(fit)$id, true = true_ranef)
summary(compare_random_effects)

ggpairs(compare_random_effects,
        lower = "blank",
        upper = list(continuous = function(data, mapping) {
          ggplot(data = data, mapping = mapping) +
            geom_point() +
            geom_abline(intercept = 0, slope = 1)}
        ))

