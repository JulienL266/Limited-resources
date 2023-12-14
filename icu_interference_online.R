# Reading icu data (semi-parametric estimators)
# Have issues with certain values not existing
library(haven)
library(dplyr)
library(SuperLearner)
data <- read_dta("~/Downloads/icu_pseudo_data.dta")
set.seed(2023)
n <- nrow(data)

# Defining variables
Y <- (1 - data$dead90)
data <- select(data, !c(dead7, dead28, dead90))

A <- data$icu_bed
data <- select(data, !c(icu_bed))

## define kappa (might change later)
#kappa <- mean(A)/2
kappa <- mean(A)
#kappa <- mean(data$icu_recommend)
#kappa <- 1


## Removing uninteresting variables
data <- select(data, !c(id))

#chosen variables, may change, needs to be low-dimensional
L <- data[,c("age", "male", "sofa_score","sepsis_dx", "winter", "periarrest", "out_of_hours", "news_score", "icnarc_score","site")]

#Estimation 
l_n <- ceiling(sqrt(n))
  
S <- (n-l_n)/l_n

cutoffs <- l_n + ceiling(S*(0:((n - l_n)/S)))

D_tilde <- function(d, Q, g, y,a,w,j){
  D_1 <- ((a*d(w,j) + (1-a)*(1-d(w,j)))/(a*g(1,w,j) + (1-a)*g(0,w,j)))*(y - Q(a,w,j))
  return(D_1 + Q(d(w,j),w,j))
}
#precompute models for each cutoff
#g_n <- function(a,w,j){
 # if(j == l_n + 1){
  #  k <- j
  #}else{
   # k <- cutoffs[max(which(cutoffs < j))] + 1
  #}
  #A_j <- A[1:k-1]
  #Y_j <- Y[1:k-1]
  #L_j <- L[1:k-1,]
  #empirical approach
  #A_jw <- c()
  #for(i in 1:(k-1)){
   # w_i <- L_j[i,]
    #if(sum(w_i == w) == length(w)){
     # A_jw <- c(A_jw, A_j[i])
    #}
  #}
  #if(length(A_jw) == 0){
   # return(0)
  #}
  #return(a*mean(A_jw) + (1-a)*mean(1-A_jw))
  #model approach
  #g_j <- SuperLearner(A_j, L_j, family = binomial, SL.library = "SL.gam")
  #return(a*predict(g_j, w)$pred + (1-a)*(1-predict(g_j,w)$pred))
#}
#precomputing to make function better
g_cutoffs <- list()
pb <- txtProgressBar(min = 1, max = n, initial = 1, style = 3)
for(j in (cutoffs + 1)){
  if(j == l_n + 1){
    k <- j
  }else{
    k <- cutoffs[max(which(cutoffs < j))] + 1
  }
  A_j <- A[1:k-1]
  Y_j <- Y[1:k-1]
  L_j <- L[1:k-1,]
  g_cutoffs[[j]] <- SuperLearner(A_j, L_j, family = binomial, SL.library = "SL.gam")
  setTxtProgressBar(pb,j)
}
g_n <- function(a,w,j){
  return(a*as.numeric(predict(g_cutoffs[[cutoffs[max(which(cutoffs < j))] + 1]], w)$pred) + (1-a)*(1-as.numeric(predict(g_cutoffs[[cutoffs[max(which(cutoffs < j))] + 1]],w)$pred)))
}
#Q_n <- function(a,w,j){
 # if(j == l_n + 1){
  #  k <- j
  #}else{
   # k <- cutoffs[max(which(cutoffs < j))] + 1
  #}
  #A_j <- A[1:k-1]
  #Y_j <- Y[1:k-1]
  #L_j <- L[1:k-1,]
  #empirical approach
  #Y_ja <- Y_j[which(A_j == a)]
  #L_ja <- L_j[which(A_j == a),]
  #if(length(Y_ja == 0)){
   # return(0)
  #}
  #Y_jaw <- c()
  #for(i in 1:length(L_ja)){
   # if(sum(L_ja[i,] == w) == length(w)){
    #  Y_jaw <- c(Y_jaw, Y_ja[i])
    #}
  #}
  #if(length(Y_jaw) == 0){
    #return(0)
  #}
  #return(mean(Y_jaw))
  #X_j <- cbind(A_j, L_j)
  #Q_j <- SuperLearner(Y_j, X_j, family = binomial, SL.library = "SL.gam")
  #return(predict(Q_j, cbind(data.frame(A = a), w))$pred)
#}
Q_cutoffs <- list()
pb <- txtProgressBar(min = 1, max = n, initial = 1, style = 3)
for(j in (cutoffs + 1)){
  if(j == l_n + 1){
    k <- j
  }else{
    k <- cutoffs[max(which(cutoffs < j))] + 1
  }
  A_j <- A[1:k-1]
  Y_j <- Y[1:k-1]
  L_j <- L[1:k-1,]
  X_j <- cbind(A_j, L_j)
  Q_cutoffs[[j]] <- SuperLearner(Y_j, X_j, family = binomial, SL.library = "SL.gam")
  setTxtProgressBar(pb,j)
}
Q_n <- function(a,w,j){
  cov <- cbind(data.frame(A_j = a), as.data.frame(w))
  return(as.numeric(predict(Q_cutoffs[[cutoffs[max(which(cutoffs < j))] + 1]], cov)$pred))
}
## need to change to our optimal regime
eta_n <- rep(NA,n)
tau_n <- rep(NA,n)
pb <- txtProgressBar(min = l_n + 1, max = 13011, initial = l_n + 1, style = 3)
for(j in (l_n + 1):n){
  Y_tilde <- (2*A[1:(j-1)] - 1)*(Y[1:(j-1)] - mean(Y[1:(j-1)]))/(A[1:(j-1)]*predict(g_n, L[1:(j-1),])$pred + (1-A[1:(j-1)])*(1-predict(g_n,L[1:(j-1),])$pred))
  ### Data adaptative
  Q_b <- SuperLearner(Y_tilde, L[1:(j-1),], SL.library = "SL.gam")
  ## Estimating d_0 (might change in new setting)
  eta_n[j] <- quantile(predict(Q_b,L[1:(j-1),])$pred, probs = c(1-kappa)) 
  tau_n[j] <- max(0,eta_n[j])
  setTxtProgressBar(pb,j)
}
d_n <- function(w,j){
  return(as.integer(Q_n(1,w,j) - Q_n(0,w,j) > tau_n[j]))
}

D_n <- function(y,a,w,j){
  return(D_tilde(d_n, Q_n, g_n, y, a, w, j))
}

sigma_n <- function(j){
  #res <- c()
  #for(i in 1:(j-1)){
   # res <- c(res, D_n(Y[i],A[i],L[i,],j))
  #}
  res <- D_n(Y[1:(j-1)], A[1:(j-1)], L[1:(j-1),],j)
  return(sd(res))
}

Psi_hat <- 0
Gamma_n <- 0
pb <- txtProgressBar(min = 1, max = 13011, initial = 1, style = 3)
for(j in (l_n + 1):n){
  setTxtProgressBar(pb,j)
  Psi_hat <- Psi_hat + D_n(Y[j], A[j], L[j,],j)/sigma_n(j)
  if(is.na(Psi_hat)){
    print(j)
  }
  if(is.na(D_n(Y[j], A[j], L[j,],j))){
    print("D_n")
  }
  if(is.na(sigma_n(j))){
    print("sigma is NA")
  }
  if(sigma_n(j) == 0){
    print("sigma = 0")
  }
  Gamma_n <- Gamma_n + 1/sigma_n(j)
  setTxtProgressBar(pb,j)
}
Gamma_n <- Gamma_n/(n - l_n)
Psi_hat <- Psi_hat/(Gamma_n*(n - l_n))
#value function
Psi_hat
#two-sided CI
Psi_hat + c(-qnorm(0.975)*1/(Gamma_n*sqrt(n - l_n)), qnorm(0.975)*1/(Gamma_n*sqrt(n - l_n)))
