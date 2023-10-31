# Reading icu data
library(haven)
library(dplyr)
data <- read_dta("~/Downloads/icu_pseudo_data.dta")
set.seed(2023)
n <- nrow(data)

# Defining variables
Y <- (1 - data$dead7)
data <- select(data, !c(dead7, dead28, dead90))

A <- data$icu_bed
data <- select(data, !c(icu_bed))

## define kappa (might change later)
kappa <- mean(A)
#kappa <- 0.2
#kappa <- 0.5

## Removing uninteresting variables
data <- select(data, !c(id))

#chosen variables, may change, follows Wang, Qi and Shi (2022)
L <- data[,c("age", "male", "sofa_score", "open_beds_cmp")]

# Implementing 2016 Luedtke and van der Laan algo
library(SuperLearner)
X <- cbind(A,L)
## Estimating g_0 (Luedtke and van der Laan assume they know it)
### Data adaptative
g_n <- SL.gam(A, L, family = binomial)
#g_n <- SL.nnet(A, L, family = binomial)
### Parametric
#g_n <- SL.bayesglm(A, L, family = binomial)
#g_n <- SL.glm(A, L, family = binomial)
#g_n <- SL.glm.interaction(A, L, family = binomial, obsWeights = NULL)
#g_n <- SL.mean(A, L, family = binomial)
#g_n <- SL.step(A, L, family = binomial)
#g_n <- SL.step.interaction(A, L, family = binomial)
#g_n <- SL.step.forward(A, L, family = binomial)

#g_n <- glm(A~., data = X, family = "binomial")

## Estimating Q_0

### Data adaptative
Q_n <- SL.gam(Y, X, family = binomial)
#Q_n <- SL.nnet(Y, X, family = binomial)
### Parametric
#Q_n <- SL.bayesglm(Y, X, family = binomial)
#Q_n <- SL.glm(Y, X, family = binomial)
#Q_n <- SL.glm.interaction(Y, X, family = binomial)
#Q_n <- SL.mean(Y, X, family = binomial)
#Q_n <- SL.step(Y, X, family = binomial)
#Q_n <- SL.step.interaction(Y, X, family = binomial)
#Q_n <- SL.step.forward(Y, X, family = binomial)
#Q_n <- glm(Y~., data = X, family = "binomial")

## Estimating Q_{b,o}
Y_tilde <- (2*A - 1)*(Y - mean(Y))/(A*predict(g_n, L) + (1-A)*(1-predict(g_n,L))) + mean(Y) #should be g_0 instead of g_n, if it is known
### Data adaptative
Q_b <- SL.gam(Y_tilde, L)
#Q_b <- SL.nnet(Y_tilde, L)
### Parametric
#Q_b <- SL.bayesglm(Y_tilde, L)
#Q_b <- SL.glm(Y_tilde, L)
#Q_b <- SL.glm.interaction(Y_tilde, L)
#Q_b <- SL.mean(Y_tilde, L)
#Q_b <- SL.step(Y_tilde, L)
#Q_b <- SL.step.interaction(Y_tilde, L)
#Q_b <- SL.step.forward(Y_tilde, L)
#Q_b <- lm(Y~., data = L)



## Estimating d_0 (might change in new setting)
eta_n <- -quantile(-predict(Q_b,L), probs = c(kappa)) #P_n(Q_n(L) > tau) = P_n(-Q_n(L) <= -tau)
tau_n <- max(0,eta_n)
d_n <- function(l){
  if(predict(Q_b,l)[1] > tau_n){
    return(1)
  }else{
    return(0)
  }
}

## TMLE procedure
H <- function(a,l){
  return(a*d_n(l) + (1-a)*(1-d_n(l)))/(a*predict(g_n, l)[1] + (1-a)*(1-predict(g_n,l))[1])
}
Cov <- H(A,L)
logit <- function(p){return(log(p/(1-p)))}
Off <- logit(predict(Q_n,X))
fm <- glm(Y~ Off + Cov -1, family = "binomial")
epsilon_n <- fm$coefficients[2]
expit <- function(x){return(1/(1 + exp(-x)))}
Q_star <- function(a,l){
  Off_st <- logit(predict(Q_n,cbind(data.frame(A = a), l)))
  Cov_st <- H(a,l)
  return(expit(Off_st + epsilon_n*Cov_st))
}

## Estimating Psi
Psi_hat = mean(predict(Q_n,cbind(data.frame(A = rep(0,n)),L))[1]*(1-d_n(L)) + predict(Q_n,cbind(data.frame(A = rep(1,n)),L))[1]*d_n(L))
Psi_hat
## Confidence interval
### Estimating sigma_0
sigma_n <- sqrt(mean(((A*d_n(L) + (1-A)*(1-d_n(L)))*(Y - predict(Q_n,X)[1])/((A*predict(g_n, L)[1] + (1-A)*(1-predict(g_n,L)[1]))) 
              + predict(Q_n, cbind(data.frame(A = d_n(L)), L))[1]
              - mean(predict(Q_n, cbind(data.frame(A = d_n(L)), L))[1]) - tau_n*(d_n(L) - kappa))^2))

Psi_hat + c(-qnorm(0.95)*sigma_n/sqrt(n), qnorm(0.95)*sigma_n/sqrt(n))

