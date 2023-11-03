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
D_tilde <- function(d, Q, g, y,a,w,j){
  D_1 <- ((a*d(w,j) + (1-a)*(1-d(w,j)))/(a*g(1,w,j) + (1-a)*g(0,w,j)))*(y - Q(a,w,j))
  return(D_1 + Q(d(w,j),w,j))
}
g <- function(a,w,j){
  A_j <- A[1:j]
  Y_j <- Y[1:j]
  L_j <- L[1:j,]
  A_jw <- A_j[which(L_j == w)]
  return(a*mean(A_jw) + (1-a)*mean(1-A_jw))
}

Q(a,w,j){
  A_j <- A[1:j]
  Y_j <- Y[1:j]
  L_j <- L[1:j,]
  Y_ja <- Y_j[which(A_j == a)]
  L_ja <- L_j[which(A_j == a)]
  Y_jaw <- Y_ja[which(L_ja == w)]
  return(mean(Y_jaw))
}

