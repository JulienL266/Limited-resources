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

## Categorizing sofa score as in the guidelines
data$sofa_score <- cut(data$sofa_score, breaks = c(-1,7,11,14))

#chosen variables, may change, follows Wang, Qi and Shi (2022)
L <- data[,c("age", "male", "sofa_score")]

# Selecting finite sample 
n_samp <- 20
samp <- sample(1:n, size = n_samp)
A_samp <- A[samp]
L_samp <- L[samp,]
Y_samp <- Y[samp]

# Ranking function Gamma
red <- which(L_samp$sofa_score == "(-1,7]")
order_red <- sample(red, size = length(red))

yellow <- which(L_samp$sofa_score == "(7,11]")
order_yellow <- sample(yellow, size = length(yellow))

blue <- which(L_samp$sofa_score == "(11,14]")
order_blue <- sample(blue, size = length(blue))

order <- c(order_red, order_yellow, order_blue)

Gamma <- function(i_samp){
  return(which(order == i_samp))
}

# IPW estimator
q_n <- function(a,l){
  B_n <- 0
  L_n <- 0
  for(i in 1:n){
    if(L[i,] == l){
      L_n <- L_n + 1
      if(A[i] == a){
        B_n <- B_n + 1
      }
    }
  }
  return(B_n/L_n)
}

q_star <- function(a,l){
  B_n_samp <- 0
  L_n_samp <- 0
  for(i in 1:n_samp){
    if(L_samp[i,] == l){
      L_n_samp <- L_n_samp + 1
      if(A_samp[i] == a){
        B_n_samp <- B_n_samp + 1
      }
    }
  }
  if(L_n_samp == 0){
    return(0)
  }else{
    return(B_n_samp/L_n_samp)
  }
}

Val.IPW <- 0 
for(i in 1:n){
  Val.IPW <- Val.IPW + Y[i]*q_star(A[i], L[i])/q_n(A[i], L[i])
}

## Variance estimation with bootstrap

# Parametric g-formula estimator

## Variance estimation with bootstrap