# Reading icu data
library(haven)
data <- read_dta("~/Downloads/icu_pseudo_data.dta")
set.seed(2023)
n <- nrow(data)

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
Y_tilde <- (2*A - 1)*(Y - mean(Y))/g_n(A,W) + mean(Y)