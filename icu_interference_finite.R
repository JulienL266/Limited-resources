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
#precalculate these!!
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
#pre-calculate these!! (see below)
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
    if(length(red) == 0){
      if(length(treated) > 0){return(1)}else{return(0)}
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
    if(length(yellow) == 0){
      if(length(treated) > length(red)){return(1)}else{return(0)}
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
    if(length(yellow) == 0){
      if(length(treated) > length(red) + length(yellow)){return(1)}else{return(0)}
    }
    return(a_count/length(blue))
  }
}

Val.IPW <- 0
pb <- txtProgressBar(min = 0, max = n, initial = 0, style = 3)
for(i in 1:n){
  if(q_n(A[i], L[i,]) == 0){
    print("problem!")
  }
  if(is.nan(Val.IPW)){
    print("NaN problem")
    break
    }
  Val.IPW <- Val.IPW + Y[i]*q_star(A[i], L[i,])/(n*q_n(A[i], L[i,]))
  setTxtProgressBar(pb,i)
}


## Variance estimation with bootstrap
### bootstrapping elements in super cluster, as finite cluster propensity is known
### n goes to infinity, but our finite sample is fixed
B <- 100
Val.IPW.boot <- rep(NA,B)
pb <- txtProgressBar(min = 0, max = B, initial = 0, style = 3)
for(b in 1:B){
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
  setTxtProgressBar(pb,b)
}
### Normal approximation
sigma <- sd(Val.IPW.boot)
Val.IPW + c(-qnorm(0.975)*sigma/sqrt(n), qnorm(0.975)*sigma/sqrt(n))
### Bootstrap CI
2*Val.IPW - c(quantile(Val.IPW.boot, 0.975), quantile(Val.IPW.boot, 0.025))

# Parametric g-formula estimator[NOT RUN YET]
### see Theorem 1 and pages 17-18 for an equation describing the g-formula
X <- L$sofa_score
Q_Y <- glm(Y~ A + X)

f <- function(y,a,l){
  return((predict(Q_Y, list(A = a, X = l$sofa_score))*y + (1-y)*(1-predict(Q_Y, list(A = a, X = l$sofa_score))))*q_star(a,l)*length(which(L$sofa_score == l$sofa_score))/n)
}

Val.g <- 0
for(i in 1:n){
  Val.g <- Val.g + Y[i]*f(Y[i], A[i], L[i,])
}
## Variance estimation with bootstrap

B <- 100
Val.g.boot <- rep(NA,B)
pb <- txtProgressBar(min = 0, max = B, initial = 0, style = 3)
for(b in 1:B){
  boot_samp <- sample(1:n, size = n, replace = TRUE)
  A_boot <- A[boot_samp]
  L_boot <- L[boot_samp,]
  Y_boot <- Y[boot_samp]
  X.boot <- L_boot$sofa_score
  Q_Y.boot <- glm(Y_boot~ A_boot + X.boot)
  p_red <- length(which(L_boot$sofa_score == "(-1,7]"))/n
  p_yellow <- length(which(L_boot$sofa_score == "(7,11]"))/n
  p_blue <- length(which(L_boot$sofa_score == "(11,14]"))/n
  f.boot <- function(y,a,l){
    if(L_boot$sofa_score == "(-1,7]"){
      return(predict(Q_Y.boot, data.frame(A_boot = a, X.boot = l$sofa_score))*q_star(a,l)*p_red)
    }else if(L_boot$sofa_score == "(7,11]"){
      return(predict(Q_Y.boot, data.frame(A_boot = a, X.boot = l$sofa_score))*q_star(a,l)*p_yellow)
    }else if(L_boot$sofa_score == "(11,14]"){
      return(predict(Q_Y.boot, data.frame(A_boot = a, X.boot = l$sofa_score))*q_star(a,l)*p_blue)
    }
}
  
  
  Val.g.boot[b] <- 0
  for(i in 1:n){
    Val.g.boot[b] <- Val.g.boot[b] + Y_boot[i]*f.boot(Y_boot[i], A_boot[i], L_boot[i,])
  }
  setTxtProgressBar(pb,b)
}
### Normal approximation
sigma.g <- sd(Val.g.boot)
Val.g + c(-qnorm(0.975)*sigma.g/sqrt(n), qnorm(0.975)*sigma.g/sqrt(n))
### Bootstrap CI
2*Val.IPW - c(quantile(Val.g.boot, 0.975), quantile(Val.g.boot, 0.025))
