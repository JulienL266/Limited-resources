# Reading icu data (semi-parametric estimators)
library(haven)
library(dplyr)
data <- read_dta("~/Downloads/icu_pseudo_data.dta")
set.seed(2023)
n <- nrow(data)

# Defining variables
Y <- (1 - data$dead90)
data <- select(data, !c(dead7, dead28, dead90))

A <- data$icu_bed
data <- select(data, !c(icu_bed))

## define kappa (might change later)
kappa <- mean(A)
#kappa <- 0.2
#kappa <- 0.6
#kappa <- 0.05

## Removing uninteresting variables
data <- select(data, !c(id))

#chosen variables, may change, follows Wang, Qi and Shi (2022)
L <- data[,c("age", "male", "sofa_score", "open_beds_cmp")]

#Estimation 
l_n <- ceiling(sqrt(n))
  
D_tilde <- function(d, Q, g, y,a,w,j){
  D_1 <- ((a*d(w,j) + (1-a)*(1-d(w,j)))/(a*g(1,w,j) + (1-a)*g(0,w,j)))*(y - Q(a,w,j))
  return(D_1 + Q(d(w,j),w,j))
}
g_n <- function(a,w,j){
  A_j <- A[1:j]
  Y_j <- Y[1:j]
  L_j <- L[1:j,]
  A_jw <- A_j[which(L_j == w)]
  return(a*mean(A_jw) + (1-a)*mean(1-A_jw))
}

Q_n <- function(a,w,j){
  A_j <- A[1:j]
  Y_j <- Y[1:j]
  L_j <- L[1:j,]
  Y_ja <- Y_j[which(A_j == a)]
  L_ja <- L_j[which(A_j == a),]
  Y_jaw <- Y_ja[which(L_ja == w)]
  return(mean(Y_jaw))
}

d_n <- function(w,j){
  return(as.integer(Q_n(1,w,j) - Q_n(0,w,j) > 0))
}

D_n <- function(y,a,w,j){
  return(D_tilde(d_n, Q_n, g_n, y, a, w, j))
}

sigma_n <- function(j){
  res <- c()
  for(i in l_n:j-1){
    res <- c(res, D_n(Y[i],A[i],L[i,],j))
  }
  return(sd(res))
}

Psi_hat <- 0
Gamma_n <- 0
for(j in (l_n + 1):n){
  Psi_hat <- Psi_hat + D_n(Y[j], A[j], L[j,],j)/sigma_n(j)
  Gamma_n <- Gamma_n + 1/sigma_n(j)
}
Gamma_n <- Gamma_n/(n - l_n)
Psi_hat <- Psi_hat/(Gamma_n*(n - l_n))
#value function
Psi_hat
#two-sided CI
Psi_hat + c(-qnorm(0.975)*1/(Gamma_n*sqrt(n - l_n)), qnorm(0.975)*1/(Gamma_n*sqrt(n - l_n)))