# Reading icu data (semi-parametric estimators)
library(haven)
library(dplyr)
data <- read_dta("~/Downloads/icu_pseudo_data.dta")
set.seed(2023)
n <- nrow(data)

# Defining variables
Y <- (1 - data$dead90)
data <- select(data, !c(dead7, dead28, dead90))

A <- data$icu_accept
data <- select(data, !c(icu_accept))

## define kappa (might change later)
kappa <- mean(A)/2
#kappa <- mean(A)
#kappa <- mean(data$icu_recommend)
#kappa <- 1

## Removing uninteresting variables
data <- select(data, !c(id))

#chosen variables, may change, follows Wang, Qi and Shi (2022)
L <- data[,c("age", "male", "sofa_score")]

# Implementing 2016 Luedtke and van der Laan algo
library(SuperLearner)
X <- cbind(A,L)
## Estimating g_0 (Luedtke and van der Laan assume they know it)
### Data adaptative
g_n <- SuperLearner(A, L, family = binomial, SL.library = "SL.gam")
#g_n <- SuperLearner(A, L, family = binomial, SL.library = "SL.nnet")
### Parametric
#g_n <- SuperLearner(A, L, family = binomial, SL.library = "SL.bayesglm")
#g_n <- SuperLearner(A, L, family = binomial, SL.library = "SL.glm")
#g_n <- SuperLearner(A, L, family = binomial, SL.library = "SL.glm.interaction")
#g_n <- SuperLearner(A, L, family = binomial, SL.library = "SL.mean")
#g_n <- SuperLearner(A, L, family = binomial, SL.library = "SL.step")
#g_n <- SuperLearner(A, L, family = binomial, SL.library = "SL.step.interaction")
#g_n <- SuperLearner(A, L, family = binomial, SL.library = "SL.step.forward")


## Estimating Q_0
### Data adaptative
Q_n <- SuperLearner(Y, X, family = binomial, SL.library = "SL.gam")
#Q_n <- SuperLearner(Y, X, family = binomial, SL.library = "SL.nnet")
### Parametric
#Q_n <- SuperLearner(Y, X, family = binomial, SL.library = "SL.bayesglm")
#Q_n <- SuperLearner(Y, X, family = binomial, SL.library = "SL.glm")
#Q_n <- SuperLearner(Y, X, family = binomial, SL.library = "SL.glm.interaction)
#Q_n <- SuperLearner(Y, X, family = binomial, SL.library = "SL.mean")
#Q_n <- SuperLearner(Y, X, family = binomial, SL .library = "SL.step")
#Q_n <- SuperLearner(Y, X, family = binomial, SL.library = "SL.step.interaction")
#Q_n <- SuperLearner(Y, X, family = binomial, SL.library = "SL.step.forward)


## Estimating Q_{b,o}
Y_tilde <- (2*A - 1)*(Y - mean(Y))/(A*predict(g_n, L)$pred + (1-A)*(1-predict(g_n,L)$pred)) #+ mean(Y) #, don't need mean(Y), see reference in Luedtke and van der Laan (2016) which is not the same target #should be g_0 instead of g_n, if it is known
### Data adaptative
Q_b <- SuperLearner(Y_tilde, L, SL.library = "SL.gam")
#Q_b <- SuperLearner(Y_tilde, L, SL.library = "SL.nnet")
### Parametric
#Q_b <- SuperLearner(Y_tilde, L, SL.library = "SL.bayesglm")
#Q_b <- SuperLearner(Y_tilde, L, SL.library = "SL.glm")
#Q_b <- SuperLearner(Y_tilde, L, SL.library = "SL.glm.interaction")
#Q_b <- SuperLearner(Y_tilde, L, SL.library = "SL.mean")
#Q_b <- SuperLearner(Y_tilde, L, SL.library = "SL.step")
#Q_b <- SuperLearner(Y_tilde, L, SL.library = "SL.step.interaction")
#Q_b <- SuperLearner(Y_tilde, L, SL.library = "SL.step.forward")



## Estimating d_0 (might change in new setting)
eta_n <- quantile(predict(Q_b,L)$pred, probs = c(1-kappa)) #P_n(Q_n(L) > tau) = P_n(-Q_n(L) <= -tau)
tau_n <- max(0,eta_n)
d_n <- function(l){
  return(as.integer(predict(Q_b,l)$pred > tau_n))
}

## TMLE procedure
H <- function(a,l){
  return(a*d_n(l) + (1-a)*(1-d_n(l)))/(a*predict(g_n, l)$pred + (1-a)*(1-predict(g_n,l)$pred))
}
Cov <- H(A,L)
logit <- function(p){return(log(p/(1-p)))}
Off <- logit(predict(Q_n,X)$pred)
fm <- glm(Y~ Off + Cov -1, family = "binomial")
epsilon_n <- fm$coefficients[2]
expit <- function(x){return(1/(1 + exp(-x)))}
Q_star <- function(a,l){
  Off_st <- logit(predict(Q_n,cbind(data.frame(A = a), l))$pred)
  Cov_st <- H(a,l)
  return(expit(Off_st + epsilon_n*Cov_st))
}

## Estimating Psi
Psi_hat = mean(Q_star(rep(0,n),L)*(1-d_n(L)) + Q_star(rep(1,n),L)*d_n(L))
Psi_hat
## Confidence interval
### Estimating sigma_0
sigma_n <- sqrt(mean(((A*d_n(L) + (1-A)*(1-d_n(L)))*(Y - predict(Q_n,X)$pred)/((A*predict(g_n, L)$pred + (1-A)*(1-predict(g_n,L)$pred))) 
              + predict(Q_n, cbind(data.frame(A = d_n(L)), L))$pred
              - mean(predict(Q_n, cbind(data.frame(A = d_n(L)), L))$pred) - tau_n*(d_n(L) - kappa))^2))

Psi_hat + c(-qnorm(0.975)*sigma_n/sqrt(n), qnorm(0.975)*sigma_n/sqrt(n))

#Survival curve
kappa_Ao2 <- mean(A)/2
kappa_A <- mean(A)
kappa_R <- mean(data$icu_recommend)
kappa_1 <- 1

eta_Ao2 <- quantile(predict(Q_b,L)$pred, probs = c(1-kappa_Ao2))
tau_Ao2 <- max(0, eta_Ao2)
eta_A <- quantile(predict(Q_b,L)$pred, probs = c(1-kappa_A))
tau_A <- max(0, eta_A)
eta_R <- quantile(predict(Q_b,L)$pred, probs = c(1-kappa_R))
tau_R <- max(0, eta_R)
eta_1 <- quantile(predict(Q_b,L)$pred, probs = c(1-kappa_1))
tau_1 <- max(0, eta_1)

Survival <- function(x){
  return(mean(as.integer(predict(Q_b,L)$pred > x)))
}

curve(Survival, from = -1, to = 1, xlab = "x", ylab = "y")