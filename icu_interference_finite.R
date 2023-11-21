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
L <- data[,c("sofa_score")] 

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
    if(L[i,]$sofa_score == l$sofa_score){
      L_n <- L_n + 1
      if(A[i] == a){
        B_n <- B_n + 1
      }
    }
  }
  return(B_n/L_n)
}

#q_star <- function(a,l){ #may need to implement ranking here, right now it's not there
 # B_n_samp <- 0
  #L_n_samp <- 0
  #for(i in 1:n_samp){
   # if(L_samp[i,]$sofa_score == l$sofa_score){
    #  L_n_samp <- L_n_samp + 1
     # if(A_samp[i] == a){
      #  B_n_samp <- B_n_samp + 1
      #}
    #}
  #}
  #if(L_n_samp == 0){
   # return(0)
  #}else{
   # return(B_n_samp/L_n_samp)
  #}
#}
q_star <- function(a,l){
  treated <- order[1:floor(kappa*n_samp)]
  if(l$sofa_score == "(-1,7]"){
    a_count <- 0
    for(i in 1:length(red)){
      if(i <= length(treated)){
        if(a == 1){
          a_count <- a_count + 1
        }
      }else{
        if(a == 0){
          a_count <- a_count + 1
        }
      }
    }
    return(a_count/length(red))
  }else if(l$sofa_score == "(7,11]"){
    a_count <- 0
    for(i in 1:length(yellow)){
      if(i + length(red) <= length(treated)){
        if(a == 1){
          a_count <- a_count + 1
        }
      }else{
        if(a == 0){
          a_count <- a_count + 1
        }
      }
    }
    return(a_count/length(yellow))
  }else if(l$sofa_score == "(11,14]"){
    a_count <- 0
    for(i in 1:length(blue)){
      if(i + length(red) + length(yellow) <= length(treated)){
        if(a == 1){
          a_count <- a_count + 1
        }
      }else{
        if(a == 0){
          a_count <- a_count + 1
        }
      }
    }
    return(a_count/length(blue))
  }
}

Val.IPW <- 0
pb <- txtProgressBar(min = 0, max = n, initial = 0, style = 3)
for(i in 1:n){
  setTxtProgressBar(pb,i)
  Val.IPW <- Val.IPW + Y[i]*q_star(A[i], L[i,])/q_n(A[i], L[i,])
}
Val.IPW <- Val.IPW/n

## Variance estimation with bootstrap[NOT RUN YET]
### bootstrapping elements in super cluster, as finite cluster propensity is known
### n goes to infinity, but our finite sample is fixed
B <- 1000
Val.IPW.boot <- rep(NA,B)
pb <- txtProgressBar(min = 0, max = B, initial = 0, style = 3)
for(b in 1:B){
  setTxtProgressBar(pb,b)
  boot_samp <- sample(1:n, size = n, replace = TRUE)
  A_boot <- A[boot_samp]
  L_boot <- L[boot_samp,]
  Y_boot <- Y[boot_samp]
  
  q_boot <- function(a,l){
    B_n <- 0
    L_n <- 0
    for(i in 1:n){
      if(L_boot[i,]$sofa_score == l$sofa_score){
        L_n <- L_n + 1
        if(A_boot[i] == a){
          B_n <- B_n + 1
        }
      }
    }
    return(B_n/L_n)
  }
  
  Val.IPW.boot[b] <- 0 
  for(i in 1:n){
    Val.IPW.boot[b] <- Val.IPW.boot[b] + Y_boot[i]*q_star(A_boot[i], L_boot[i,])/q_boot(A_boot[i], L_boot[i,])
  }
  Val.IPW.boot[b] <- Val.IPW.boot[b]/n
}
sigma <- sd(Val.IPW.boot)

Val.IPW + c(-qnorm(0.975)*sd/sqrt(n), qnorm(0.975)*sd/sqrt(n))

# Parametric g-formula estimator

## Variance estimation with bootstrap