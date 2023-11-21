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
n <- 20
samp <- sample(1:n, size = 20)
A <- A[samp]
L <- L[samp,]
Y <- Y[samp]

# Ranking function Gamma
red <- which(L$sofa_score == "(-1,7]")
order_red <- sample(red, size = length(red))

yellow <- which(L$sofa_score == "(7,11]")
order_yellow <- sample(yellow, size = length(yellow))

blue <- which(L$sofa_score == "(11,14]")
order_blue <- sample(blue, size = length(blue))

order <- c(order_red, order_yellow, order_blue)

Gamma <- function(i){
  return(which(order == i))
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

## Variance estimation with bootstrap

# Parametric g-formula estimator

## Variance estimation with bootstrap