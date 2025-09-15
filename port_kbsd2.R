
# Simulating causal settings with 1,2 and 3 confounders with sparse data for some strata
# and checking if different diagnostics detect these strata with positivity violations.
library(simcausal)
set.seed(15082025)

sem1 <- DAG.empty() +
  node("L", distr = "rnorm", mean = 15, sd = 5) +
  node("A", distr = "rbern", prob = ifelse(L >= 10 & L <= 15, 0.02,         # if L in [10,15], then most are A=0 -> "cluster"
                                           ifelse(L > 15 & L <= 25, 0.98,   # if L in [15,25], most are A=1
                                                  ifelse(L < 10 | L > 25, 0.1, 0))))
dag1 <- set.DAG(sem1)
obs1 <- sim(dag1, rndseed = 05082025, n = 1000)
plot(obs1[-1])

# PoRT ---
setwd("C:/Users/victo/OneDrive/Desktop/Uni/Statistik/BA_thesis/positivity-violations-diag")
source("data/port_utils.R")
lst1 <- list(port = NULL, port_risca = NULL)
a_values <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.1)
gruber <- 5/(sqrt(nrow(obs1))*log(nrow(obs1)))
b_values <- c(0.001, 0.01, gruber, 0.05, 0.1)
port(A = "A", cov.quanti = "L", cov.quali = NULL, data = obs1, alpha = 0.05, beta = 0.01, gamma = 1)

# KBSD ---
source("kbsd.R")
o1 <- obs1[-1]
o1_1 <- o1
o1_1["A"] <- 1
o1_2 <- o1
o1_2["A"] <- 0
res1 <- kbsd(data = o1, int_data_list = list(o1_1, o1_2), disthalf_vec=c(L=4.9, A=0.5*0.5), # HALF of SD for A
             plot.out = F)
res1_plot <- kbsd(data = o1,
                  int_data_list = list(o1_1, o1_2),
                  disthalf_vec=c(L=4.9, A=0.5*0.5))
res1_plot
# makes sense that for IV=1 (which equals A=0) have lower EDP overall, as points are more scattered -> more "singles" without support

# check what L values the low-support points in A=0 have
subset1 <- res1[res1$shift == 1,]
subset1 <- subset1$diagnostic < 50 # for outliers
l_values1 <- obs1[subset1, "L"]  # the obs (with L, A values) which have low EDP (among A=0)
diag_values1 <- res1[subset1,] %>% filter(shift == 1)
diag_values1 <- diag_values1[, "diagnostic"]
plot(obs1[subset1,"L"])
plot(l_values1, diag_values1)
# 1. all outliers are obs with L < 10 | L > 30
# 2. the more extreme the L values, the less support
# WEIRD: why points in [20,30] not detected as with low support? ~~~~~~~~~~~~~~~~~~ TUE ~~~~~~~~~~~~~~~~~~~~

# check what L values the low-support points in A=1 have
subset2 <- res1[res1$shift == 2,]
subset2 <- subset2$diagnostic < 50 # for outliers
l_values2 <- obs1[subset2, "L"]  # the obs (with L, A values) which have low EDP (among A=1)
diag_values2 <- res1[subset2,] %>% filter(shift == 2)
diag_values2 <- diag_values2[, "diagnostic"]
plot(obs1[subset2,"L"])
plot(l_values2, diag_values2)
# WEIRD: why returns EDP values for points that don't exist in A=1? A=1 in obs1 does not have L<0 | L >28 ~~~~~~~~ TUE ~~~~~~

