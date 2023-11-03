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
D_tilde <- function(d, Q, g, y,a,w){
  D_1 <- (a*d(w) + (1-a)*(1-d(w)))/(a*g(1,w) + (1-a)*g(0,w))
}