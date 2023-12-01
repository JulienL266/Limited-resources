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
#kappa <- mean(A)/2
kappa <- mean(A)
#kappa <- mean(data$icu_recommend)
#kappa <- 1

## Removing uninteresting variables
data <- select(data, !c(id))

## Categorizing sofa score as in the guidelines
data$sofa_score <- cut(data$sofa_score, breaks = c(-1,7,11,14))

#chosen variables, may change, follows Wang, Qi and Shi (2022)
L <- data[,c("age", "male", "sofa_score", "sepsis_dx", "winter", "periarrest", "out_of_hours", "news_score", "icnarc_score")]
L$age <- cut(L$age, breaks = c(17,quantile(L$age, c(1/3, 2/3)),104)) #categorizing age into 3 quartiles
L$news_score <- cut(L$news_score, breaks = c(-1,quantile(L$news_score, c(1/3,2/3)),20))
L$icnarc_score <- cut(L$icnarc_score, breaks = c(-1,quantile(L$icnarc_score, c(1/3, 2/3)),53))

# Selecting finite sample size
n_samp <- 20


# IPW estimator
## precalculation step
q_n <- function(a,l){
  B_n <- 0
  L_n <- 0
  for(i in 1:n){
    l_eq <- TRUE
    for(j in 1:ncol(L)){
      if(L[i,j] != l[1,j]){
        l_eq <- FALSE
      }
    }
    if(l_eq){
      L_n <- L_n + 1
      if(A[i] == a){
        B_n <- B_n + 1
      }
    }
  }
  return(B_n/L_n)
}

dims <- c(3,2,3,2,2,2,2,3,3,2)
q_n.image <- array(rep(NA, prod(dims)), dim = dims)
for(i_age in 1:length(levels(L$age))){
  for(i_male in 1:2){
    for(i_sofa_score in 1:length(levels(L$sofa_score))){
      for(i_sepsis_dx in 1:2){
        for(i_winter in 1:2){
          for(i_periarrest in 1:2){
            for(i_out_of_hours in 1:2){
              for(i_news_score in 1:length(levels(L$news_score))){
                for(i_icnarc_score in 1:length(levels(L$icnarc_score))){
                  for(i_a in 1:2){
                  q_n.image[i_age, i_male, i_sofa_score, i_sepsis_dx, i_winter, i_periarrest, i_out_of_hours, i_news_score, i_icnarc_score, i_a] <- q_n(i_a - 1, data.frame(age = levels(L$age)[i_age], male = i_male - 1, 
                                               sofa_score = levels(L$sofa_score)[i_sofa_score],
                                               sepsis_dx = i_sepsis_dx - 1, winter = i_winter - 1,
                                               periarrest = i_periarrest - 1,
                                               out_of_hours = i_out_of_hours - 1,
                                               news_score = levels(L$news_score)[i_news_score],
                                               icnarc_score = levels(L$icnarc_score)[i_icnarc_score]))
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
q_n <- function(a,l){
  return(q_n.image[which(levels(L$age) == l$age), l$male + 1, which(levels(L$sofa_score) == l$sofa_score),
                   l$sepsis_dx + 1, l$winter + 1, l$periarrest + 1, l$out_of_hours + 1,
                   which(levels(L$news_score) == l$news_score),
                   which(levels(L$icnarc_score) == l$icnarc_score), a + 1])
}



## Pre-calculation
n_red <- length(which(L$sofa_score == "(-1,7]"))
n_yellow <- length(which(L$sofa_score == "(7,11]"))
n_blue <- length(which(L$sofa_score == "(11,14]"))

### check that probs make sense
p_1.red <- 0
for(r in 0:n_red){
  for(y in 0:n_yellow){
    for(b in 0: n_blue){
      if(r + y + b == n_samp){
        if(r <= kappa*n_samp){
          p_1.red <- p_1.red + (choose(n_red,r)*choose(n_yellow,y)*choose(n_blue,b)/choose(n, n_samp))
        }else{
            p_1.red <- p_1.red + (floor(kappa*n_samp)/r)*(choose(n_red,r)*choose(n_yellow,y)*choose(n_blue,b)/choose(n, n_samp))
        }
      }
    }
  }
}
p_1.yellow <- 0
for(r in 0:n_red){
  for(y in 0:n_yellow){
    for(b in 0: n_blue){
      if(r + y + b == n_samp){
        if(r + y <= kappa*n_samp){
          p_1.yellow <- p_1.yellow + (choose(n_red,r)*choose(n_yellow,y)*choose(n_blue,b)/choose(n, n_samp))
        }else{
          if(floor(kappa*n_samp) > r){
            if(y > 0){
              p_1.yellow <- p_1.yellow + ((floor(kappa*n_samp) - r)/y)*(choose(n_red,r)*choose(n_yellow,y)*choose(n_blue,b)/choose(n, n_samp))
            }else{
              p_1.yellow <- p_1.yellow + (choose(n_red,r)*choose(n_yellow,y)*choose(n_blue,b)/choose(n, n_samp))
            }
          }
        }
      }
    }
  }
}
p_1.blue <- 0
for(r in 0:n_red){
  for(y in 0:n_yellow){
    for(b in 0: n_blue){
      if(r + y + b == n_samp){
        if(r + y + b <= kappa*n_samp){
          p_1.blue<- p_1.blue + (choose(n_red,r)*choose(n_yellow,y)*choose(n_blue,b)/choose(n, n_samp))
        }else{
          if(floor(kappa*n_samp) > r + y){
            if(b > 0){
              p_1.blue<- p_1.blue + ((floor(kappa*n_samp) - (r+y))/b)*(choose(n_red,r)*choose(n_yellow,y)*choose(n_blue,b)/choose(n, n_samp))
            }else{
              p_1.blue<- p_1.blue + (choose(n_red,r)*choose(n_yellow,y)*choose(n_blue,b)/choose(n, n_samp))
            }
          }
        }
      }
    }
  }
}


q_star <- function(a,l){
  if(l$sofa_score == "(-1,7]"){
    if(a == 0){
      return(1 - p_1.red)
    }else{
      return(p_1.red)
    }
  }
  if(l$sofa_score == "(7,11]"){
    if(a == 0){
      return(1 - p_1.yellow)
    }else{
      return(p_1.yellow)
    }
  }
  if(l$sofa_score == "(11,14]"){
    if(a == 0){
      return(1 - p_1.blue)
    }else{
      return(p_1.blue)
    }
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
      l_eq <- TRUE
      for(j in 1:ncol(L)){
        if(L_boot[i,j] != l[1,j]){
          l_eq <- FALSE
        }
      }
      if(l_eq){
        L_n <- L_n + 1
        if(A_boot[i] == a){
          B_n <- B_n + 1
        }
      }
    }
    return(B_n/L_n)
  }
  
  dims <- c(3,2,3,2,2,2,2,3,3,2)
  q_boot.image <- rep(NA,n)
  Val.IPW.boot[b] <- 0 
  for(i in 1:n){
    if(is.na(q_boot.image[boot_samp[i]])){
    q_boot.image[boot_samp[i]] <- q_boot(A_boot[i], L_boot[i,])
    }
    Val.IPW.boot[b] <- Val.IPW.boot[b] + Y_boot[i]*q_star(A_boot[i], L_boot[i,])/q_boot.image[boot_samp[i]]
  }
  Val.IPW.boot[b] <- Val.IPW.boot[b]/n
  setTxtProgressBar(pb,b)
}
### Normal approximation
sigma <- sd(Val.IPW.boot)
Val.IPW + c(-qnorm(0.975)*sigma/sqrt(n), qnorm(0.975)*sigma/sqrt(n))
### Bootstrap CI
2*Val.IPW - c(quantile(Val.IPW.boot, 0.975), quantile(Val.IPW.boot, 0.025))

# Parametric g-formula estimator
### see Theorem 1 and pages 17-18 for an equation describing the g-formula
X <- data[,c("age", "male", "sofa_score", "sepsis_dx", "winter", "periarrest", "out_of_hours", "news_score", "icnarc_score","site")]
Q_Y <- glm(Y~., data = cbind(A,X), family = "binomial")

f <- function(y,a,l){
  return((predict(Q_Y, cbind(data.frame(A = a), l))*y + (1-y)*(1-predict(Q_Y, cbind(data.frame(A = a), l))))*q_star(a,l))
}

Val.g <- 0
for(i in 1:n){
  Val.g <- Val.g + Y[i]*f(Y[i],A[i],X[i,])
}
Val.g <- Val.g/n
## Variance estimation with bootstrap
set.seed(2023)
B <- 100
Val.g.boot <- rep(NA,B)
pb <- txtProgressBar(min = 0, max = B, initial = 0, style = 3)
for(b in 1:B){
  boot_samp <- sample(1:n, size = n, replace = TRUE)
  A_boot <- A[boot_samp]
  L_boot <- L[boot_samp,]
  Y_boot <- Y[boot_samp]
  X.boot <- data[boot_samp,c("age", "male", "sofa_score", "sepsis_dx", "winter", "periarrest", "out_of_hours", "news_score", "icnarc_score","site")]
  Q_Y.boot <- glm(Y_boot~ ., data = cbind(A_boot, X.boot), family = "binomial")

  f.boot <- function(y,a,l){
    return((y*predict(Q_Y.boot, cbind(data.frame(A_boot = a), l)) + (1-y)*(1 - predict(Q_Y.boot, cbind(data.frame(A_boot = a), l))))*q_star(a,l))
}
  
  Val.g.boot[b] <- 0
  for(i in 1:n){
    Val.g.boot[b] <- Val.g.boot[b] + Y_boot[i]*f.boot(Y_boot[i],A_boot[i],X.boot[i,])
  }
  Val.g.boot[b] <- Val.g.boot[b]/n
  setTxtProgressBar(pb,b)
}
### Normal approximation
sigma.g <- sd(Val.g.boot)
Val.g + c(-qnorm(0.975)*sigma.g/sqrt(n), qnorm(0.975)*sigma.g/sqrt(n))
### Bootstrap CI
2*Val.g - c(quantile(Val.g.boot, 0.975), quantile(Val.g.boot, 0.025))
