# Reading icu data
library(haven)
library(dplyr)
data <- read_dta("~/Downloads/icu_pseudo_data.dta")
set.seed(2023)
n <- nrow(data)

# Defining variables
Y <- 100*(1 - data$dead7)
data <- select(data, !c(dead7, dead28, dead90))

A <- data$icu_bed
data <- select(data, !c(icu_bed))

## define kappa (might change later)
kappa <- mean(icu_bed)

## Removing uninteresting variables
data <- select(data, !c(id))

#chosen variables, may change, follows Wang, Qi and Shi (2022)
L <- select(data, c(age, male, sofa_score, open_beds_cmp))


# Implementing 2016 Luedtke and van der Laan algo
library(SuperLearner)

## Estimating g_0 (Luedtke and van der Laan assume they know it)
### Data adaptative
g_n <- SL.gam(A, L, family = binomial)
g_n <- SL.nnet(A, L, family = binomial)
### Parametric
g_n <- SL.bayesglm(A, L, family = binomial)
g_n <- SL.glm(A, L, family = binomial)
g_n <- SL.glm.interaction(A, L, family = binomial)
g_n <- SL.mean(A, L, family = binomial)
g_n <- SL.step(A, L, family = binomial)
g_n <- SL.step.interaction(A, L, family = binomial)
g_n <- SL.step.forward(A, L, family = binomial)


## Estimating Q_0
X <- cbind(A,L)
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
Y_tilde <- (2*A - 1)*(Y - mean(Y))/(A*predict(g_n, L) + (1-A)*(1-predict(g_n,L))) + mean(Y) #should be g_0 instead of g_n, if it is known
### Data adaptative
Q_b <- SL.gam(Y_tilde, L)
Q_b <- SL.nnet(Y_tilde, L)
### Parametric
Q_b <- SL.bayesglm(Y_tilde, L)
Q_b <- SL.glm(Y_tilde, L)
Q_b <- SL.glm.interaction(Y_tilde, L)
Q_b <- SL.mean(Y_tilde, L)
Q_b <- SL.step(Y_tilde, L)
Q_b <- SL.step.interaction(Y_tilde, L)
Q_b <- SL.step.forward(Y_tilde, L)

## Estimating d_0 (might change in new setting)
d_n <- function(l){
  if(predict(Q_b,l) > 0){
    return(1)
  }else{
    return(0)
  }
}

## Estimating Psi
Psi_hat = mean(predict(Q_n,cbind(data.frame(A = rep(0,n)),L))*(1-d(L)) + predict(Q_n,cbind(data.frame(A = rep(1,n)),L))*d(L))

## Confidence interval
### Estimating sigma_0
eta_n <- 
tau_n <- max(0,eta_n)
sigma_n <- mean(I(A = d_n(L))*(Y - predict(Q_n,X))/((A*predict(g_n, L) + (1-A)*(1-predict(g_n,L)))) 
              + predict(Q_n, cbind(data.frame(A = d_n(L)), L))
              - mean(predict(Q_n, cbind(data.frame(A = d_n(L)), L))) - tau_n*(d_n(L) - kappa))

