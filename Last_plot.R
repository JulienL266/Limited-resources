# Reading icu data
library(haven)
library(dplyr)
data <- read_dta("~/Downloads/icu_pseudo_data.dta")
set.seed(2023)
n <- nrow(data)

## Defining variables
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
data$site <- as.factor(data$site)
data$news_score <- as.factor(data$news_score)

## chosen variables, may change, follows Wang, Qi and Shi (2022)
L <- data[,c("age", "male", "sofa_score", "sepsis_dx", "winter", "periarrest", "out_of_hours", "news_score", "icnarc_score","site")]
n_samp <- 20
#simulations are done according to SEM in following -> change to param. g-formula and IPW


# Intervention density
## Pre-calculation
n_red <- length(which(L$sofa_score == "(-1,7]"))
n_yellow <- length(which(L$sofa_score == "(7,11]"))
n_blue <- length(which(L$sofa_score == "(11,14]"))

### computing probs for different color groups
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

# Dependencies part
library(dplyr)
library(pracma)
library(mc2d)
library(dplyr)
library("multicool")
library(xtable)
library("latex2exp")
library(parallel)
expit<-function(x){1/(1+exp(-x))}

# Functions.R part
fullsim<-function(n, kappa,  p_Y_i, p_L_i=NULL, seed=1, conditional=F, comp_in=NULL,  estimate=F,   B){
  
  Start<-Sys.time()
  gcomp<-NULL
  if(estimate==F){
    gcomp<-gcomp(n, kappa, conditional, comp_in, p_L_i, p_Y_i)
  }	
  gcomp_time<-Sys.time()-Start
  
  Start<-Sys.time()
  sim<-sim(n, kappa, B, conditional, comp_in, seed, gcomp[["p_A_gplus_l"]], p_L_i, p_Y_i, estimate)
  sim_time<-Sys.time()-Start
  
  Means<-vector(4, mode="list")
  names(Means)<-c("Y", "Y_g", "Y_gprime", "TrueY_g")
  Means[1]<-sim[["Y"]]        %>% {rowSums(.)/n} %>% mean(.)
  Means[2]<-sim[["Y_g"]]      %>% {rowSums(.)/n} %>% mean(.)
  Means[3]<-sim[["Y_gprime"]] %>% {rowSums(.)/n} %>% mean(.)
  Means[4]<-gcomp[["E_Y_gplus"]]
  
  Counts<-vector(3, mode="list")
  names(Counts)<-c("Y", "Y_g", "Y_gprime")
  Counts[[1]]<-table_counts(sim[["Y"]]       , n)
  Counts[[2]]<-table_counts(sim[["Y_g"]]     , n)
  Counts[[3]]<-table_counts(sim[["Y_gprime"]], n)
  
  results<-vector(5, mode="list")
  names(results)<-c("gform", "sim", "times", "means", "counts")
  if(estimate==F){
    results[[1]]<-gcomp
  }
  results[[2]]<-sim
  results[[3]]<-c(gcomp_time, sim_time)
  results[[4]]<-Means
  results[[5]]<-Counts
  return(results)
}

gcomp<-function(n, kappa, conditional=F, comp_in=NULL, p_L_i=NULL, p_Y_i){
  
  #Preliminaries
  kappa_n=ceiling(kappa*n)
  if(conditional){
    comp<-{c(round(comp_in*n)[1:3], n-sum(round(comp_in*n)[1:3]))} %>%  {sapply(FUN=function(x){rep(x, .[x])}, 1:4)} %>% unlist
  }
  
  ##Computing g-formula#########################################################################################
  #$Q_{L_i}$:    
  ##Creating PMF for $L_i$ 
  #p_L_i<-rep(0.25, 4)
  #p_L_i<-c(1/5, 1/10, 3/5, 1/10)
  
  if(conditional){
    #$\mathcal{P}_L(l_i)$:    
    ##Creating $\mathcal{P}_L(l_i)$ for each $l_i \in \{1,2,3,4\}$
    L_grid_l<- rbind(
      as.vector(table(comp)) - c(1,0,0,0),
      as.vector(table(comp)) - c(0,1,0,0),
      as.vector(table(comp)) - c(0,0,1,0),
      as.vector(table(comp)) - c(0,0,0,1)
    )
    p_L_bb<-1
    #$P^*_{n, \mathbb{l^*}}(\Lambda(L_i^*)=c)$
    ##For each element in $\mathcal{P}_L^*(l^*_i)$, computing fractions of individuals in each rank-group (that is $P^*_{n, \mathbb{l^*}}(\Lambda(L_i^*)=c)$), for each $c \in \{1,2,3,4\}$, and for each $l^*_i \in \{1,2,3,4\}$.
    P_star_l_bb_lambda_grid<-as.vector(table(comp)[order(sapply(1:4, f_LambdaL_i))])/n
    
    
    #$\omega_{n, \mathbb{l^*}}$
    ##For each element in $\mathcal{P}_L^*(l^*_i)$, computing border group (that is $\omega_{n, \mathbb{l^*}}$), for each $l^*_i \in \{1,2,3,4\}$.
    omega_l_bb<-{{P_star_l_bb_lambda_grid %>% t %>% cumsum %>% t} >(kappa_n/n)} %>% which %>% min
    
    
    #$g^*$
    ##Compute probability of treatment under regime, for each $l^*_i \in \{1,2,3,4\}$.
    p_A_gplus_l_bb<-NULL
    p_A_gplus_l<-rep(NA, 4)
    for(i in 1:4){
      p_A_gplus_l[i]<-
        if        (f_LambdaL_i(i)<omega_l_bb){
          1
        } else if (f_LambdaL_i(i)>omega_l_bb){
          0
        } else if (omega_l_bb==1) {
          ((kappa_n/n)) / 
            P_star_l_bb_lambda_grid[omega_l_bb]
        } else  											  {
          ((kappa_n/n) - {P_star_l_bb_lambda_grid %>% t %>% cumsum %>% t}[omega_l_bb - 1]) / 
            P_star_l_bb_lambda_grid[omega_l_bb]
        }
    }
  }else{
    #$\mathcal{P}_L(l_i)$:    
    ##Creating $\mathcal{P}_L(l_i)$ for each $l_i \in \{1,2,3,4\}$
    L_grid_l<-{expand.grid(0:(n-1) ,0:(n-1), 0:(n-1), 0:(n-1))} %>% {.[which(rowSums(.)==(n-1)), ]}
    
    #$\mathcal{P}_L^*(l^*_i)$:     
    ##Creating PMF for each element in $\mathcal{P}_L^*(l^*_i)$ (that is P(\mathbb{L}=\mathbb{l} \mid L^*_i=l_i)), for each $l^*_i \in \{1,2,3,4\}$ 
    p_L_bb<-{{ {p_L_i^{L_grid_l %>% t}} %>% t} %>% log %>% rowSums %>% exp}*	
      mcmapply(X=1:nrow(L_grid_l), FUN=function(X){factorial(n-1)/(factorial(L_grid_l[X, 1])*factorial(L_grid_l[X, 2])*factorial(L_grid_l[X, 3])*factorial(L_grid_l[X, 4]))}, mc.cores = 11)
    
    
    #$P^*_{n, \mathbb{l^*}}(\Lambda(L_i^*)=c)$
    ##For each element in $\mathcal{P}_L^*(l^*_i)$, computing fractions of individuals in each rank-group (that is $P^*_{n, \mathbb{l^*}}(\Lambda(L_i^*)=c)$), for each $c \in \{1,2,3,4\}$, and for each $l^*_i \in \{1,2,3,4\}$.
    P_star_l_bb_lambda_grid<-vector(4, mode="list")
    P_star_l_bb_lambda_grid[[1]]<-{{expand.grid(0:n ,0:n, 0:n, 0:n)} %>% {.[which(rowSums(.)==n & .[, 1]>=1), ]}}[, order(sapply(1:4, f_LambdaL_i))]/n
    P_star_l_bb_lambda_grid[[2]]<-{{expand.grid(0:n ,0:n, 0:n, 0:n)} %>% {.[which(rowSums(.)==n & .[, 2]>=1), ]}}[, order(sapply(1:4, f_LambdaL_i))]/n
    P_star_l_bb_lambda_grid[[3]]<-{{expand.grid(0:n ,0:n, 0:n, 0:n)} %>% {.[which(rowSums(.)==n & .[, 3]>=1), ]}}[, order(sapply(1:4, f_LambdaL_i))]/n
    P_star_l_bb_lambda_grid[[4]]<-{{expand.grid(0:n ,0:n, 0:n, 0:n)} %>% {.[which(rowSums(.)==n & .[, 4]>=1), ]}}[, order(sapply(1:4, f_LambdaL_i))]/n
    
    #$\omega_{n, \mathbb{l^*}}$
    ##For each element in $\mathcal{P}_L^*(l^*_i)$, computing border group (that is $\omega_{n, \mathbb{l^*}}$), for each $l^*_i \in \{1,2,3,4\}$.
    omega_l_bb<-matrix(NA, nrow(P_star_l_bb_lambda_grid[[1]]), 4)
    for(i in 1:4){
      omega_l_bb[, i]<-mcmapply(X=1:nrow(P_star_l_bb_lambda_grid[[i]]), FUN=
                                  function(X){
                                    {{P_star_l_bb_lambda_grid[[i]][X, ] %>% t %>% cumsum %>% t} >(kappa_n/n)} %>% which %>% min
                                  }
                                , mc.cores = 11)
    }
    
    
    #$g^*_{n, \maathbb{l^*}}$
    ##For each element in $\mathcal{P}_L^*(l^*_i)$, compute probability of treatment under regime, for each $l^*_i \in \{1,2,3,4\}$.
    p_A_gplus_l_bb<-matrix(NA, nrow(L_grid_l), 4)
    for(i in 1:4){
      p_A_gplus_l_bb[, i]<-mcmapply(X=1:nrow(P_star_l_bb_lambda_grid[[i]]), FUN=
                                      function(X){
                                        if        (f_LambdaL_i(i)<omega_l_bb[X, i]){
                                          1
                                        } else if (f_LambdaL_i(i)>omega_l_bb[X, i]){
                                          0
                                        } else if (omega_l_bb[X, i]==1) {
                                          ((kappa_n/n)) / 
                                            P_star_l_bb_lambda_grid[[i]][X, omega_l_bb[X, i]]
                                        } else  											  {
                                          ((kappa_n/n) - {P_star_l_bb_lambda_grid[[i]][X, ] %>% t %>% cumsum %>% t}[omega_l_bb[X, i] - 1]) / 
                                            P_star_l_bb_lambda_grid[[i]][X, omega_l_bb[X, i]]
                                        }
                                      }
                                    , mc.cores = 11)
    }
    
    #$g^*$
    ##Compute probability of treatment under regime, for each $l^*_i \in \{1,2,3,4\}$.
    p_A_gplus_l<-{p_A_gplus_l_bb*p_L_bb} %>% colSums
  }
  #$Q_{Y}
  ##Compute conditional expectations of Y
  p_Y_a_gplus_l<-matrix(p_Y_i, 2, 4)
  
  
  #$\Psi^{g, AB}_n$
  ##Compute g-formula.	
  E_Y_gplus<-0
  if(conditional){
    for(i in 0:1){
      for(j in 1:4){
        E_Y_gplus<- E_Y_gplus + 
          p_Y_a_gplus_l[i+1, j]*
          ((p_A_gplus_l[j])^(i)*(1-p_A_gplus_l[j])^(1-i))*
          as.vector(table(comp)/n)[j]
      }
    }		
  }else{
    for(i in 0:1){
      for(j in 1:4){
        E_Y_gplus<- E_Y_gplus + p_Y_a_gplus_l[i+1, j]*
          ((p_A_gplus_l[j])^(i)*(1-p_A_gplus_l[j])^(1-i))*
          p_L_i[j]
        
      }
    }		
  }
  
  
  #Storing results
  gcomp<-vector(9, mode="list")
  names(gcomp)<-c("p_L_i","L_grid_l","p_L_bb","P_star_l_bb_lambda_grid","omega_l_bb","p_A_gplus_l_bb","p_A_gplus_l","p_Y_a_gplus_l","E_Y_gplus")
  gcomp[["p_L_i"]]<-p_L_i
  gcomp[["L_grid_l"]]<-L_grid_l
  gcomp[["p_L_bb"]]<-p_L_bb
  gcomp[["P_star_l_bb_lambda_grid"]]<-P_star_l_bb_lambda_grid
  gcomp[["omega_l_bb"]]<-omega_l_bb
  gcomp[["p_A_gplus_l_bb"]]<-p_A_gplus_l_bb
  gcomp[["p_A_gplus_l"]]<-p_A_gplus_l
  gcomp[["p_Y_a_gplus_l"]]<-p_Y_a_gplus_l
  gcomp[["E_Y_gplus"]]<-E_Y_gplus	
  return(gcomp)
}


sim<-function(n, kappa, B, conditional=F, comp_in=NULL, seed=1, p_A_gplus_l_in2, p_L_i=NULL, p_Y_i, estimate_in=F){
  
  #Preliminaries
  kappa_n=ceiling(kappa*n)
  if(estimate_in==F){
    p_A_gplus_l<-p_A_gplus_l_in2
  }
  if(conditional){
    comp<-{c(round(comp_in*n)[1:3], n-sum(round(comp_in*n)[1:3]))} %>%  {sapply(FUN=function(x){rep(x, .[x])}, 1:4)} %>% unlist
  }
  
  set.seed(seed)
  
  #Initializing simulatiion matrices
  epsilon_L<-matrix(NA, B, n)
  epsilon_Z<-matrix(NA, B, n)
  epsilon_R<-matrix(NA, B, n)
  epsilon_A<-matrix(NA, B, n)
  epsilon_Y<-matrix(NA, B, n)
  delta_A<-matrix(NA, B, n)
  L<-matrix(NA, B, n)
  LambdaL<-matrix(NA, B, n)
  Z<-matrix(NA, B, n)
  R<-matrix(NA, B, n)
  R_gplus<-matrix(NA, B, n)
  A<-matrix(0, B, n)
  A_gplus<-matrix(0, B, n)
  A_gprimeplus<-matrix(0, B, n)
  Y<-matrix(NA, B, n)
  Y_g<-matrix(NA, B, n)
  Y_gprime<-matrix(NA, B, n)
  
  #Simulating error terms as all independent uniforms on [0,1]
  epsilon_L<-matrix(runif(B*n, 0,1), B, n)
  epsilon_Z<-matrix(runif(B*n, 0,1), B, n)
  epsilon_R<-matrix(runif(B*n, 0,1), B, n)
  epsilon_A<-matrix(runif(B*n, 0,1), B, n)
  epsilon_Y<-matrix(runif(B*n, 0,1), B, n)
  delta_A  <-matrix(runif(B*n, 0,1), B, n)
  
  #Simulating factual data according to structural equations
  if(conditional){
    L      <- {rep(comp, B)} %>% {matrix(., B, n, byrow=T)}
  }else{
    L      <-	{mapply(f_L_i, epsilon_L_i=as.vector(epsilon_L), MoreArgs = list(p_L_i_in=p_L_i))} %>% {matrix(., B, n)}
  }
  LambdaL<-	{sapply(as.vector(L), f_LambdaL_i)} %>% {matrix(., B, n)}
  Z      <-	epsilon_Z
  
  if(conditional){
    R      <- f_R(L[1, ], epsilon_R[1, ]) %>%  rep(B) %>% matrix(B, n, byrow=T)
    A      <- mclapply(1:B, 
                       function(b){
                         for(i in 1:n){
                           A[b, which(R[b, ]==i)]<-f_A_i(R[b, ],  R[b, which(R[b, ]==i)], A[b, ], L[b, ], L[b, which(R[b, ]==i)], epsilon_A[b, which(R[b, ]==i)], kappa_n, n)
                         }
                         return(A[b, ])
                       } 
                       , mc.cores = 11) %>% unlist %>% {matrix(., B, n, byrow=T)}
    
    R_gplus<- f_R_gplus(LambdaL[1, ], Z[1, ]) %>%  rep(B) %>% matrix(B, n, byrow=T)
    for(i in 1:n){ A_gplus[1, which(R_gplus[1, ]==i)]<-f_A_i_gplus(i, kappa_n)}
    A_gplus<-A_gplus[1, ] %>%  rep(B) %>% matrix(B, n, byrow=T) 
    
  } else{
    R      <-{lapply(1:B, function(b){f_R(L[b, ], epsilon_R[b, ])})} %>% unlist %>% {matrix(., B, n, byrow=T)}
    A      <- mclapply(1:B, 
                       function(b){
                         for(i in 1:n){
                           A[b, which(R[b, ]==i)]<-f_A_i(R[b, ],  R[b, which(R[b, ]==i)], A[b, ], L[b, ], L[b, which(R[b, ]==i)], epsilon_A[b, which(R[b, ]==i)], kappa_n, n)
                         }
                         return(A[b, ])
                       }
                       , mc.cores = 11) %>% unlist %>% {matrix(., B, n, byrow=T)}
    R_gplus<- {lapply(1:B, function(b){f_R_gplus(LambdaL[b, ], Z[b, ])})} %>% unlist %>% {matrix(., B, n, byrow=T)}
    A_gplus<- mclapply(1:B, 
                       function(b){
                         for(i in 1:n){
                           A_gplus[b, which(R_gplus[b, ]==i)]<-f_A_i_gplus(i, kappa_n)
                         }
                         return(A_gplus[b, ])
                       }
                       , mc.cores = 11) %>% unlist %>% {matrix(., B, n, byrow=T)}		
  }
  Y      <- {mcmapply(f_Y_i, L_i=as.vector(L), A_i=as.vector(A), epsilon_Y_i=as.vector(epsilon_Y), MoreArgs=list(p_Y_i_in=p_Y_i), mc.cores = 11) +0} %>% {matrix(., B, n)}
  #Simulating counterfactual data for regime $g$ according to structural equations
  Y_g         <- {mcmapply(f_Y_i, L_i=as.vector(L), A_i=as.vector(A_gplus)     , epsilon_Y_i=as.vector(epsilon_Y), MoreArgs=list(p_Y_i_in=p_Y_i), mc.cores = 11) +0} %>% {matrix(., B, n)}
  
  if(estimate_in==F){
    #Simulating counterfactual data for regime $g$ according to structural equations
    A_gprimeplus<- {mcmapply(f_A_i_gprimeplus, delta_A_i=as.vector(delta_A), L_i=as.vector(L), MoreArgs = list(p_A_gplus_l_in = p_A_gplus_l), mc.cores = 11) +0} %>% {matrix(., B, n)}
    Y_gprime    <- {mcmapply(f_Y_i, L_i=as.vector(L), A_i=as.vector(A_gprimeplus), epsilon_Y_i=as.vector(epsilon_Y), MoreArgs=list(p_Y_i_in=p_Y_i), mc.cores = 11) +0} %>% {matrix(., B, n)}
  }	
  #Storing results
  simresults<-vector(17, mode="list")
  names(simresults)<-c("epsilon_L","epsilon_Z","epsilon_R","epsilon_A","epsilon_Y","delta_A", "L","LambdaL","Z","R","R_gplus","A","A_gplus","A_gprimeplus","Y","Y_g","Y_gprime")
  simresults[["epsilon_L"]]<-epsilon_L
  simresults[["epsilon_Z"]]<-epsilon_Z
  simresults[["epsilon_R"]]<-epsilon_R
  simresults[["epsilon_A"]]<-epsilon_A
  simresults[["epsilon_Y"]]<-epsilon_Y
  simresults[["delta_A"]]<-delta_A
  simresults[["L"]]<-L
  simresults[["LambdaL"]]<-LambdaL
  simresults[["Z"]]<-Z
  simresults[["R"]]<-R
  simresults[["R_gplus"]]<-R_gplus
  simresults[["A"]]<-A
  simresults[["A_gplus"]]<-A_gplus
  simresults[["A_gprimeplus"]]<-A_gprimeplus
  simresults[["Y"]]<-Y
  simresults[["Y_g"]]<-Y_g
  simresults[["Y_gprime"]]<-Y_gprime
  return(simresults)
}


table_counts<-function(matrix, n){
  counttab<- {matrix %>% rowSums(.) %>% table %>% prop.table}
  fincounttab<-rep(0, n+1)
  names(fincounttab)<-as.character(0:n)
  fincounttab[names(counttab)]<-counttab
  return(fincounttab)
}


############################################################################################################
############################################################################################################
#2. Structural equations
############################################################################################################
############################################################################################################

#<self-explanatory>

############################################################################################################
#Descriptions
############################################################################################################

#<self-explanatory>

############################################################################################################
#Functions
############################################################################################################

f_L_i<-function(epsilon_L_i, p_L_i_in){
  (1)*(                                  epsilon_L_i<    p_L_i_in[1]   ) + 
    (2)*(epsilon_L_i>=    p_L_i_in[1]    & epsilon_L_i<sum(p_L_i_in[1:2])) + 
    (3)*(epsilon_L_i>=sum(p_L_i_in[1:2]) & epsilon_L_i<sum(p_L_i_in[1:3])) +
    (4)*(epsilon_L_i>=sum(p_L_i_in[1:3]))  
}

f_LambdaL_i<-function(L_i){
  (1)*(L_i==2) + 
    (2)*(L_i==3) + 
    (3)*(L_i==4) + 
    (4)*(L_i==1)  
}


f_Z<-function(epsilon_Z){
  epsilon_Z
}

f_R<-function(L, epsilon_R){
  order(order(L + epsilon_R))
}


f_R_gplus<-function(LambdaL, Z){
  order(order(LambdaL + Z))
}

f_A_i<-function(R, R_i, A, L, L_i, epsilon_A_i, kappa_n, n){
  (kappa_n > sum(A[which(R<R_i)]))
  
  (epsilon_A_i< (5-L_i)*.2 )	
  
}



f_A_i_gplus<-function(R_i_gplus, kappa_n){
  kappa_n >= R_i_gplus
}

f_A_i_gprimeplus<-function(delta_A_i, L_i, p_A_gplus_l_in){
  delta_A_i<p_A_gplus_l_in[L_i]
}


f_Y_i<-function(L_i, A_i, epsilon_Y_i,  p_Y_i_in){
  if(A_i==0 & L_i==1){ epsilon_Y_i<matrix(p_Y_i_in, 2,4)[1,1]} 
  else if(A_i==1 & L_i==1){ epsilon_Y_i<matrix(p_Y_i_in, 2,4)[2,1]} 
  else if(A_i==0 & L_i==2){ epsilon_Y_i<matrix(p_Y_i_in, 2,4)[1,2]} 
  else if(A_i==1 & L_i==2){ epsilon_Y_i<matrix(p_Y_i_in, 2,4)[2,2]} 
  else if(A_i==0 & L_i==3){ epsilon_Y_i<matrix(p_Y_i_in, 2,4)[1,3]} 
  else if(A_i==1 & L_i==3){ epsilon_Y_i<matrix(p_Y_i_in, 2,4)[2,3]} 
  else if(A_i==0 & L_i==4){ epsilon_Y_i<matrix(p_Y_i_in, 2,4)[1,4]} 
  else                    { epsilon_Y_i<matrix(p_Y_i_in, 2,4)[2,4]} 
}



############################################################################################################
############################################################################################################
#3. Estimation functions
############################################################################################################
############################################################################################################

#A. estimate
#B. estim
#C. estim_props

############################################################################################################
#Descriptions
############################################################################################################


#A. full-sim: this function is the primary function called for simulations to evaluate MC g-formula estimators. It essentially runs in three steps. First, a large MC simulation is conducted to approximate the true value of the parameter (used as a benchmark to evaluate bias); then a number of 'observed data' clusters are simulated, each of which, on its own, might constitute the observed data in a single study. Then, a parametric g-formula estimator is implemented for each of these sample replicates, except instead of using the known parameter values in the truth, the MC simulation is run using inputs estimated from that replicate.
## (i) It outputs a results object (a list) with: 
### (a) a benchmark simulation that approximates the "truth"
### (b) a list of replication samples (to simulate the frequentist paradigm) - the estimator will be implemented for each of these samples
### (c) a list of replications of the estimator on each of the replication samples
### (d) summaries across the replication samples to evaluate estimator performance 
## (ii) There are three secondary functions it calls:
### (a) fullsim
### (b) estim
### (c) estim_props
## (iii) Arguments
### (a) n_star_sim: number of units in observed data cluster (larger than n_sim, for computational purposes)
### (b) n_sim: number of units in cluster defining the estimand
### (c) kappa_sim: <see full sim>
### (d) p_Y_i_sim: <see full sim>
### (e) p_L_i_sim: <see full sim>
### (f) seed_sim: <see full sim>
### (g) B_bench: number of simulation replication for approximating the truth
### (h) B1: number of replications
### (i) B2: number of simulations for each MC g-formula estimate


#B. estim: this is the function that computes the nuisance parameter estimates for each of the B1 replications of the estimator

#C. estim_props: this is function that estimates the properties of the estimator by computing averages across replicates.

############################################################################################################
#Functions
############################################################################################################

estimate<-function(n_star_sim, n_sim, kappa_sim, p_Y_i_sim, p_L_i_sim, seed_sim, B_bench, B1, B2){
  benchmark<-fullsim(
    B=B_bench, #Cluster size for estimand
    n=n_sim, #Cluster size for estimand
    kappa=kappa_sim, 
    p_L_i=p_L_i_sim,
    p_Y_i=p_Y_i_sim,
    estimate=T,
    seed=seed_sim
  )
  
  samples<-	fullsim(
    seed=seed_sim+1,
    B=B1, #Number of replications
    n=n_star_sim, #Number of samples for estimation in each replication (cluster size in ''study'' data)
    kappa=kappa_sim,
    p_L_i=p_L_i_sim,
    p_Y_i=p_Y_i_sim,
    estimate=T
  )[["sim"]]
  
  replications <- lapply(X=1:B1, #For each replication, compute MC g-formula estimator
                         FUN=function(b){
                           estims<-estim(b, samples)
                           fullsim(
                             seed=seed_sim+1+b,
                             B=B2, 
                             n=n_sim,
                             kappa=kappa_sim,
                             p_L_i=estims[["p_L_i"]],
                             p_Y_i=estims[["p_Y_i"]],
                             estimate=T
                           )
                         })
  
  estim_proportions<-estim_props(replications)
  
  results<-vector(mode="list")
  names(results)<-c("Benchmark", "Samples", "Replications", "Properties")
  results[[1]]<-benchmark
  results[[2]]<-samples
  results[[3]]<-replications
  results[[4]]<-estim_proportions
  return(results)
}

estim<-function(b, sim){
  estims<-vector(2, mode="list")
  names(estims)<-c("p_L_i", "p_Y_i")
  estims[["p_L_i"]]<-prop.table(table(sim[["L"]][b, ]))
  estims[["p_Y_i"]]<-c(t(prop.table(table(sim[["L"]][b, ], sim[["sim"]][["A"]][b, ], sim[["Y"]][b, ]), margin=1:2)[,,2]))
  return(estims)
}

estim_props<-function(sims){
  estims_sums<-vector(2, mode="list")
  names(estims_sums)<-c("counts", "means")
  estims_sums[[1]]<-matrix(NA, length(sims), length(sims[[1]][["counts"]][["Y_g"]]))
  estims_sums[[2]]<-rep(NA, length(sims))
  for(i in 1:length(sims)){
    estims_sums[[1]][x, ]<- sims[[x]][["counts"]][["Y_g"]]
    estims_sums[[2]][x  ]<- sims[[x]][["means" ]][["Y_g"]]
  }
  
  return(estims_sums)
}



############################################################################################################
############################################################################################################
#4. Graphing functions
############################################################################################################
############################################################################################################

#A. graph
#B. graphforEstim

############################################################################################################
#Descriptions
############################################################################################################

#A. graph: See Figures 12 and 13
#B. graphforEstim: See Figure 14

############################################################################################################
#Functions
############################################################################################################



graph<-function(sim_in, n, conditional=F){
  par(oma=c(2,2,1,1))
  
  myhist   <-list(breaks=(0:n)-0.5, counts=sim_in[["counts"]][["Y_g"]]     , density=	sim_in[["counts"]][["Y_g"]])
  myhist2  <-list(breaks=(0:n)-0.5, counts=sim_in[["counts"]][["Y_gprime"]], density=	sim_in[["counts"]][["Y_gprime"]])
  class(myhist)  <- "histogram"
  class(myhist2) <- "histogram"
  
  if(conditional==F){
    plot(myhist, col=rgb(100,0,0,max = 100, alpha = 50, names = "lt.blue"),  border=rgb(100,0,0,max = 100, alpha = 10, names = "lt.blue"), main=" ", xlab=" ", ylab=" ", ylim=c(0, max(max(sim_in[["counts"]][["Y"]]), max(sim_in[["counts"]][["Y_g"]]), max(sim_in[["counts"]][["Y_gprime"]]))
    ))
  } else{
    plot(myhist, col=rgb(100,0,0,max = 100, alpha = 50, names = "lt.blue"),  border=rgb(100,0,0,max = 100, alpha = 10, names = "lt.blue"), main=" ", xlab=" ", ylab=" ", ylim=c(0, max(max(sim_in[["counts"]][["Y_g"]]), max(sim_in[["counts"]][["Y_gprime"]]))
    ))	
  }
  
  plot(myhist2, add=T, col=rgb(0,0,100, max = 100, alpha = 40, names = "lt.pink"),  border=rgb(0,0,100, max = 100, alpha = 10, names = "lt.pink"))
  
  
  if(conditional==F){
    myhist3  <-list(breaks=(0:n)-0.5, counts=sim_in[["counts"]][["Y"]]       , density=	sim_in[["counts"]][["Y"]])
    class(myhist3) <- "histogram"
    plot(myhist3, add=T, col=rgb(0,100,0, max = 100, alpha = 70, names = "lt.pink"),  border=rgb(0,100,1, max = 100, alpha = 10, names = "lt.pink"))
    abline(v=sim_in[["means"]][["Y"]]       *n, col="green", lty=1, lwd=2)
  }
  
  abline(v=sim_in[["means"]][["Y_g"]]     *n, col="Black", lty=1, lwd=2)
  abline(v=sim_in[["means"]][["Y_gprime"]]*n, col="grey" , lty=2, lwd=2)
  abline(v=sim_in[["means"]][["TrueY_g"]] *n, col="Red"  , lty=3, lwd=2)
  
  if(conditional==F){
    legend("topright", legend=c("Deterministic (g)", TeX('Stochastic ($g^*$)'), "Deterministic (factual)"),  fill=c("red", "blue", "green"), title=paste0("n=", n))
    legend("topleft", legend=c(TeX('$\\widetilde{E}(\\sum_{i=1}^nY_{i}^g$) '), TeX('$\\widetilde{E}(\\sum_{i=1}^nY_{i}^{g^*}$) '), TeX('$\\psi^{g,AB}_n$'), TeX('$\\widetilde{E}(\\sum_{i=1}^nY_{i}$)')), lty=c(1:3, 1), lwd=c(2,2,2,2), col=c("black", "grey", "red", "green"), bty="n")
    legend("topleft", legend=c(TeX('$\\widetilde{E}(\\sum_{i=1}^nY_{i}^g$) '), TeX('$\\widetilde{E}(\\sum_{i=1}^nY_{i}^{g^*}$) '), TeX('$\\psi^{g,AB}_n \\times n$'), TeX('$\\widetilde{E}(\\sum_{i=1}^nY_{i}$)')), lty=c(1:3, 1), lwd=c(2,2,2,2), col=c("black", "grey", "red", "green"), bty="n")
    
    mtext(TeX('$\\sum_{i=1}^nY_{i}=x$'), at=.5, side=1, line=-1, cex.lab=1, las=1, col="black", outer=T)
    mtext(TeX('P($\\sum_{i=1}^nY_{i}=x$)'), at=.5, side=2, line=-1, cex.lab=1, las=3, col="black", outer=T)
    
  }else{
    legend("topright", legend=c("Deterministic (g)", TeX('Stochastic ($g^*$)')),  fill=c("red", "blue"), title=paste0("n=", n))
    legend("topleft", legend=c(TeX('$\\widetilde{E}(\\sum_{i=1}^nY_{i}^g$ | \\mathbf{L}=\\mathbf{l}) '), TeX('$\\widetilde{E}(\\sum_{i=1}^nY_{i}^{g^*}$| \\mathbf{L}=\\mathbf{l}) '), TeX('$\\psi^{g,AB}_{n, \\mathbf{l}} \\times n$')), lty=c(1:3), lwd=c(2,2,2), col=c("black", "grey", "red"), bty="n")
    
    mtext(TeX('$\\sum_{i=1}^nY_{i}=x$'), at=.5, side=1, line=-1, cex.lab=1, las=1, col="black", outer=T)
    mtext(TeX('P($\\sum_{i=1}^nY_{i}=x | \\mathbf{L}=\\mathbf{l}$)'), at=.5, side=2, line=-1, cex.lab=1, las=3, col="black", outer=T)
  }
  
}

graphforEstim<-function(sim_in, n, estims){
  par(oma=c(2,2,1,1))
  
  myhist   <-list(breaks=(0:n)-0.5, counts=sim_in[["counts"]][["Y_g"]]     , density=	sim_in[["counts"]][["Y_g"]])
  class(myhist)  <- "histogram"
  
  ylims<-c(0, max(
    max(sim_in[["counts"]][["Y_g"]]),
    max(apply(estims[[1]], 2, FUN=function(x){quantile(x, c(0.025, 0.975))}))
  ))
  
  plot(myhist, col=rgb(100,0,0,max = 100, alpha = 50, names = "lt.blue"),  border=rgb(100,0,0,max = 100, alpha = 10, names = "lt.blue"), main=" ", xlab=" ", ylab=" ", ylim=ylims
  )	
  
  abline(v=sim_in[["means"]][["Y_g"]]     *n, col="Black", lty=1, lwd=2)
  
  for(i in 0:n){
    lines(c(i,i), apply(estims[[1]], 2, FUN=function(x){quantile(x, c(0.025, 0.975))})[,i+1], col="blue", lty=1)
    lines(c(i-n/100, i+n/100), rep(apply(estims[[1]], 2, FUN=function(x){quantile(x, c(0.025, 0.975))})[1,i+1], 2),  col="blue", lty=1)
    lines(c(i-n/100, i+n/100), rep(apply(estims[[1]], 2, FUN=function(x){quantile(x, c(0.025, 0.975))})[2,i+1], 2),  col="blue", lty=1)
  }
  points(0:10, apply(estims[[1]], 2, FUN=mean), col="blue", bg="red", pch=23)
  
  lines(quantile(estims[[2]], c(0.025, 0.975))*n, rep(mean(ylims), 2), col="red", lty=1)
  lines(rep(quantile(estims[[2]], c(0.025, 0.975))[1], 2)*n, c(mean(ylims)-ylims[2]/100, mean(ylims)+ylims[2]/100), col="red", lty=1)
  lines(rep(quantile(estims[[2]], c(0.025, 0.975))[2], 2)*n, c(mean(ylims)-ylims[2]/100, mean(ylims)+ylims[2]/100), col="red", lty=1)
  points(mean(estims[[2]])*n, mean(ylims), col="red", bg="blue", pch=23)
  
  legend("topright", 
         legend=c(
           TeX('$\\widetilde{\\beta}_{x, n}$'),
           TeX('$\\widetilde{E}(\\hat{\\beta}_{n, x})$, '),
           TeX('$n \\times \\widetilde{E}(\\hat{\\alpha}_{n})$'),
           TeX('$n \\times \\tilde{\\alpha}_{n}$'),
           TeX('$\\widetilde{Q}_{2.5, 97.5}(\\hat{\\beta}_{x, n})$'),
           TeX('$n \\times \\widetilde{Q}_{2.5, 97.5}(\\hat{\\alpha}_{n})$')
         ),  
         col=c(
           rgb(100,0,0,max = 100, alpha = 10, names = "lt.blue"),
           "blue",
           "red", 
           "black",
           "blue",
           "red" 
         ),
         pt.bg=c(
           rgb(100,0,0,max = 100, alpha = 50, names = "lt.blue"),
           "red",
           "blue",
           NA,
           NA,
           NA
         ),
         pch = c(22, 23, 23, NA, NA, NA),
         lty=c(
           NA, NA, NA,1,
           1,1
         ),
         lwd=c(
           NA, NA, NA,2,
           1,1
         ),
         title=paste0("n=", n), bty="n"
  )
  
  #mtext(TeX('$\\sum_{i=1}^nY_{i}=x$'), at=.5, side=1, line=-1, cex.lab=1, las=1, col="black", outer=T)
  #mtext(TeX('P($\\sum_{i=1}^nY_{i}=x$)'), at=.5, side=2, line=-1, cex.lab=1, las=3, col="black", outer=T)
  mtext(TeX('x'), at=.5, side=1, line=-1, cex.lab=1, las=1, col="black", outer=T)
  
  
}








# Main.R part
##Marginal
sim1<-fullsim(
  B=1,
  n=5, 
  kappa=6/15,
  p_L_i=c(1/4, 1/4, 1/4, 1/4),
  p_Y_i=c(0.55, 0.85, 0.75, 0.15, 0.95, 0.15, 0.05, 0.15)
)


sim3<-fullsim(
  B=100000,
  n=50, 
  kappa=6/15,
  p_L_i=c(1/5, 1/10, 3/5, 1/10),
  p_Y_i=c(0.55, 0.85, 0.75, 0.15, 0.95, 0.15, 0.05, 0.15)
)

sim5<-fullsim(
  B=100000,
  n=50, 
  kappa=6/15,
  p_L_i=c(1/6, 1/6, 1/2, 1/6),
  p_Y_i=c(0.55, 0.85, 0.75, 0.15, 0.95, 0.15, 0.05, 0.15)
)
sim6<-fullsim(
  B=1000,
  n=50, 
  kappa=6/15,
  p_L_i=c(1/5, 1/10, 3/5, 1/10),
  p_Y_i=c(0.55, 0.85, 0.75, 0.15, 0.95, 0.15, 0.05, 0.15)
)
sim7<-fullsim(
  B=10000,
  n=20, 
  kappa=6/15,
  p_L_i=c(.1, .1, .7, .1),
  p_Y_i=c(0.55, 0.85, 0.75, 0.15, 0.95, 0.15, 0.05, 0.15)
)

##Conditional
sim4<-fullsim(
  B=100000,
  n=50, 
  kappa=6/15,
  conditional=T,
  comp_in=c(1/5, 1/10, 3/5, 1/10),
  p_Y_i=c(0.55, 0.85, 0.75, 0.15, 0.95, 0.15, 0.05, 0.15)
)




##Organizing Table 1 for Appendix F#########################################################################################

gformtable<-rbind(unlist(cbind(sim1[["gform"]][["L_grid_l"]], sim1[["gform"]][["omega_l_bb"]], sim1[["gform"]][["p_L_bb"]], sim1[["gform"]][["p_A_gplus_l_bb"]])) %>% matrix(35, 13),
                  c(rep(NA, 9), sim1[["gform"]][["p_L_i"]]),
                  c(rep(NA, 9), sim1[["gform"]][["p_A_gplus_l"]]),
                  cbind(matrix(NA,2, 9), sim1[["gform"]][["p_Y_a_gplus_l"]]),
                  c(rep(NA, 12), sim1[["gform"]][["E_Y_gplus"]])
)


print(xtable(gformtable, type = "latex", digits=c(rep(0, 9), rep(3, 5))), file = "/Users/asarvet/Dropbox/Harvard PhD/Classes/Dissertation/FinitePopID_gformtable.tex")



##Figures 1 & 2#########################################################################################

##Marginal
graph(sim3, 50)
graph(sim5, 50)

##Conditional
graph(sim4, 50, conditional=T)



##Evaluating estimator performance, Figure 3#####################################################################################
Start<-Sys.time()	
estimation_sim1<-estimate(
  n_star_sim = , 
  n_sim = 10, 
  kappa_sim = 6/15, 
  p_Y_i_sim = rep(.25,4),
  p_L_i_sim = c(0.55, 0.85, 0.75, 0.15, 0.95, 0.15, 0.05, 0.15), 
  seed_sim = 1, 
  B_bench = 100000, 
  B1 = 500, 
  B2 = 1000
)
Sys.time()-Start	



graphforEstim(estimation_sim1[["Benchmark"]], 10, estimation_sim1[["Properties"]])


gformestable<-rbind(
  c(estimation_sim1[["Benchmark"]][["means"]][["Y_g"]],estimation_sim1[["Benchmark"]][["counts"]][["Y_g"]]),
  c(mean(stimation_sim1[["Properties"]][[2]]),apply(stimation_sim1[["Properties"]][[1]], 2, FUN=mean)),
  cbind(	quantile(stimation_sim1[["Properties"]][[2]], c(0.025, 0.975)),apply(stimation_sim1[["Properties"]][[1]], 2, FUN=function(x){quantile(x, c(0.025, 0.975))}))
)

print(xtable(gformestable, type = "latex", digits=rep(3, 13)), file = "/Users/asarvet/Dropbox/Harvard PhD/Classes/Dissertation/FinitePopID_gformesttable.tex")
