#reading icu data
library(haven)
data <- read_dta("~/Downloads/icu_pseudo_data.dta")
set.seed(2023)
n <- nrow(data)

#implementing 2016 Luedtke and van der Laan algo
library(SuperLearner)
#estimating Q_0

#estimating Q_{b,o}
