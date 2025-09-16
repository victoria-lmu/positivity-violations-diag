
# Simulating causal settings with 1,2 and 3 confounders with sparse data for some strata
# and checking if different diagnostics detect these strata with positivity violations.
library(simcausal)
set.seed(15082025)


# 1) Clustered Observations ----

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
port(A = "A", cov.quanti = "L", cov.quali = NULL, data = obs1, alpha = 0.02, beta = gruber, gamma = 1)
# best choice: a = 0.02, beta = gruber -> detects 3 groups, smaller a allows to detect more &smaller groups

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
# for IV=1 (which equals A=1) overall lower EDP, but should have higher as points are closer together?
# IV = 2 (A=0) should have lower EDP as points are more scattered -> more "singles" without support

# check what L values the low-support points in A=0 have
shift1 <- res1[res1$shift == 1,]
outliers1 <- shift1$diagnostic < 50 # for outliers
l_values1 <- obs1[outliers1, "L"]  # the obs (with L, A values) which have low EDP (among A=0)
diag_values1 <- shift1[outliers1, ] # extract diag values of the outliers
plot(l_values1, diag_values1$diagnostic)
# 1. all outliers are obs with L < 10 | L > 30
# 2. the more extreme the L values, the less support
# ??: why returns EDP values for points that don't exist in A=1? A=1 in obs1 does not have L<0 | L >28
#     because created IV for all values! i.e. created o_int fr all values for A=0 AND A=1
#     and if we create o_int = (L=0, A=1), it has very low support as far from other obs! same for o_int = (L=34, A=1)
# wenn wir für jede unserer obs auf beide IV levels intervenieren (Konstruktion der cf), d.h. jede obs
# bekommt einmal A=0 und A=1, dann hätten die Punkte mit L nahe 0, nahe 30 mit A=1 sehr niedrige EDP

# check what L values the low-support points in A=0 have
shift2 <- res1[res1$shift == 2,]
outliers2 <- shift2$diagnostic < 50 # for outliers
l_values2 <- obs1[outliers2, "L"]  # to which original obs (&L values) do these outliers belong?
diag_values2 <- shift2[outliers2,] # what diag values to these outliers have
plot(l_values2, diag_values2$diagnostic)
# 1. again outliers are obs with L < 10 | L > 30
# 2. points in [20,30] detected as with low support

# effectively vergleicht kbsd die Verteilung von L innerhalb A=0 und A=1 -> wenn ähnlich, dann sollten für beide IV levels gleiche EDP
# wenn nicht ähnlich, zeigt es welche L-Werte im jeweiligen A "fehlen" (dort wo niedrige EDP)


# among A=0: interveniere auf jeden Punkt (auch aus A=1 -> einfach "runterziehen" auf A=0 Ebene)
#            und schaue wie viele andere obs in Reichweite
#            für intervenierten Punkt mit L-Wert = 24 aus A=1 im Bereich L = [20, 30] nur wenige andere Punkte in Reichweite
#            (die nächsten wären in A=0, aber da sind sehr wenige) (in A=1 sind viele in L=[20,30] aber zu weit weg bzgl A Wert
#            -> sind strenger bei Abstand in A weil disthalf Wert da kleiner (i.e. wenn wir schon kleinen US haben ist kernel-Info halbiert))

# among A=1: interveniere auf jeden Punkt (auch aus A=0) -> schon klar, dass für Punkt mit L<0 und L>30 aus A=0 "gezogen" nur wenige obs in Nähe
#            deshalb haben diese Punkte in künstlicher Intervention A=1 nur wenige EDP, generell für L > 20 immer weniger EDP da wenig obs in A=1

# essence: both PoRT & kbsd detect the critical obs, whereas kbsd is preciser as really on obs-level by providing individual EDPs,
#          whereas detecting "all critical obs" in PoRT only possible if alpha very small (to consider small groups)



# 2) Scattered Observations ----

sem2 <- DAG.empty() +
  node("L", distr = "rnorm", mean = 12, sd = 5) +
  node("A", distr = "rbern", prob = ifelse(L >= 0 & L <= 15, 0.1,  # if L in [0,15], most should be A=0, few in A=1
                                           ifelse(L > 15 & L <= 20, 0.8,        # if L in [15,20], then most are A=0, well scattered with sd=10
                                                  ifelse(L > 20 & L <= 25, 0.98,   # if L in [15,25], most are A=1
                                                         ifelse(L < 0 | L > 25, 0.1, 0)))))
dag2 <- set.DAG(sem2)
obs2 <- sim(dag2, rndseed = 05082025, n = 1000)
plot(obs2[-1])

# PoRT ---
setwd("C:/Users/victo/OneDrive/Desktop/Uni/Statistik/BA_thesis/positivity-violations-diag")
source("data/port_utils.R")
lst1 <- list(port = NULL, port_risca = NULL)
a_values <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.1)
gruber <- 5/(sqrt(nrow(obs1))*log(nrow(obs1)))
b_values <- c(0.001, 0.01, gruber, 0.05, 0.1)
port(A = "A", cov.quanti = "L", cov.quali = NULL, data = obs2, alpha = 0.03, beta = gruber, gamma = 1)
# best choice: a = 0.03, beta = gruber -> detects 5 groups that cover the 2 main violating groups

# KBSD ---
source("kbsd.R")
o2 <- obs2[-1]
o2_1 <- o2
o2_1["A"] <- 1
o2_2 <- o2
o2_2["A"] <- 0
res2 <- kbsd(data = o2, int_data_list = list(o2_1, o2_2), disthalf_vec=c(L=4.9, A=0.4*0.5), # HALF of SD for A
             plot.out = F)
res2_plot <- kbsd(data = o2,
                  int_data_list = list(o2_1, o2_2),
                  disthalf_vec=c(L=4.9, A=0.4*0.5))
res2_plot
# interesting: many obs in IV=1 (A=1) have much fewer support now
# but overall also in A=0 more outliers now (as wanted bc obs more scattered)


# assumption: points around L~10 should also have high EDP as in first simulation above: kbsd
# -> check EDP value for L~10 in obs1 (both for A=0 & A=1 ?)
# -> compare with EDP value for L~10 in obs2: prob will be similar which emphasises the limitation of kbsd
#    to distinguish informative support and uninformative support -> actually crucial for quality of extrapolation

around10_1 <- obs1[obs1$L > 9.7 & obs1$L <= 10,]
res1[res1$observation %in% around10_1$ID & res1$shift == 2,]  # all have EDP around 300

around10_2 <- obs2[obs2$L > 9.7 & obs2$L <= 10,]
res2[res2$observation %in% around10_2$ID & res2$shift == 2,]  # most have EDP around 500, but 2 have ~90

# but maybe not good example, I think distribution among A=0 needs to be even more scattered around 10, then compare again
