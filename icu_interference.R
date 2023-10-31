# Reading icu data
library(haven)
library(dplyr)
data <- read_dta("~/Downloads/icu_pseudo_data.dta")
set.seed(2023)
n <- nrow(data)

# Defining variables
Y <- 1 - data$dead7
Y <- 1 - data$dead28
Y <- 1 - data$dead90
data <- select(data, !c(dead7, dead28, dead90))

# Implementing 2016 Luedtke and van der Laan algo
library(SuperLearner)

## Estimating Q_0
### Data adaptative
Q_n <- SL.gam(Y, X, family = binomial)
Q_n <- SL.nnet(Y, X, family = binomial)
### Parametric
Q_n <- SL.bayesglm(Y, X, family = binomial)
Q_n <- SL.glm(Y, X, family = binomial)
Q_n <- SL.glm.interaction(Y, X, family = binomial)
Q_n <- SL.mean(Y, X, family = binomial)
Q_n <- SL.step(Y, X, family = binomial)
Q_n <- SL.step.interaction(Y, X, family = binomial)
Q_n <- SL.step.forward(Y, X, family = binomial)

## Estimating Q_{b,o}
Y_tilde <- (2*A - 1)*(Y - mean(Y))/g_n(A,W) + mean(Y) #should be g_0 instead of g_n, if it is known
### Data adaptative
Q_b <- SL.gam(Y_tilde, X)
Q_b <- SL.nnet(Y_tilde, X)
### Parametric
Q_b <- SL.bayesglm(Y_tilde, X)
Q_b <- SL.glm(Y_tilde, X)
Q_b <- SL.glm.interaction(Y_tilde, X)
Q_b <- SL.mean(Y_tilde, X)
Q_b <- SL.step(Y_tilde, X)
Q_b <- SL.step.interaction(Y_tilde, X)
Q_b <- SL.step.forward(Y_tilde, X)
