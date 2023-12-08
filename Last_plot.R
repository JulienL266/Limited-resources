
#Main.R part
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
