
source("setup.R")

# 2 Confounders But Balanced ----
# Baseline confounders: Age (L1) ~ N(50, 10) and Fitness (L2) dep on age -> most people in sample are not fit
# Treatment: BP medication (A) ~ Ber with logit(Age, Fitness)
set.seed(15082025)
sem2 <- DAG.empty() +
  node("L1", distr = "rnorm", mean = 50, sd = 10) +
  node("L2", distr = "rbern", prob = ifelse(L1 > 60, 0.2, 0.7)) +  # L2ness depends on L1
  node("A", distr = "rbern", prob = plogis(4*(L1 > 60) + 2.5*(1-L2) - 3*L2))
# if unfit -> P(A=1)=0.92, if old -> P(A=1) = 0.98, if unfit + old -> P(A=1)~1, if L2 -> P(A=1)=0.047
# if old & unfit -> def treated, i.e. P(A=1|old & unfit) should be high
dag2 <- set.DAG(sem2)
plotDAG(dag2)
obs2 <- sim(dag2, rndseed = 30072025, n = 1000)
table(obs2$A)  # balanced

obs2 %>% filter(L1 > 60 & A==1) %>% nrow()/obs2 %>% filter(L1 > 60) %>% nrow()  # viol #1: for beta=0.1

obs2 %>% filter(L2 == 0 & A==1) %>% nrow()/obs2 %>% filter(L2 ==0) %>% nrow()  # viol #2: for beta=0.05,0.1
obs2 %>% filter(L1 > 60 & L2 == 0 & A==1) %>% nrow()/obs2 %>% filter(L1 > 60 & L2 ==0) %>% nrow()  # viol #3: combo of the 2 above with beta = 0.01!

obs2 %>% filter(L2 == 1 & A==1) %>% nrow()/obs2 %>% filter(L2 ==1) %>% nrow()  # viol #4: for beta=0.1


# PoRT ---
source('data/port_utils.R')
a_values <- c(0.01, 0.025, 0.05, 0.1)
gruber2 <- 5/(sqrt(nrow(obs2))*log(nrow(obs2)))
b_values <- c(0.01, gruber2, 0.05, 0.1)
g_values <- 1:2
lst2 <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      # uncategorised L1
      #lst2[[paste0("gamma = ",g, ", alpha = ", a, ", beta = ", b)]] <-
      #  port(A = "A", cov.quanti = "L1", cov.quali = "L2", data = obs2, alpha = a, beta = b, gamma = g) 
      # categorised L1
      lst2[[paste0("gamma = ",g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A", cov.quanti =NULL, cov.quali = c("L1", "L2"), data = obs2, alpha = a, beta = b, gamma = g)
    }
  }
}
lst2
# uncategorised:
# gamma = 1, 2:
#   a=0.01/0.02 & b<=gruber: only viol #1, bc L1 as only cont var that can be split into smaller subgroups
#   from a=0.03, all detected (or say viol #3 covered by detection of 1,2) for b=0.1 (viol #2 from 0.05) as they should be!
#   viol #3 (intersecting group) poss to detect for gamma =2, but only covered bc of greedy cat for L1
#   also, young pp being treated rarely is detected for a=0.01! for that it's useful to do sensitivity ana

# categorised:
obs2$L1 <- cut(obs2$L1, breaks = c(0,20,40,60,90))
# gamma = 1: same as above, always detected from beta = 0.1 (viol #2 from 0.05) as should
# gamma = 2: now also intersecting group detected! bc no greedy cat by L1, rest always det for b=0.1 as wanted


# KBSD ---
source("kbsd.R")
# revert categorisation, i.e need cont L1, L2 again!
o2 <- obs2[-1]
o2_1 <- o2
o2_1$A <- 1
o2_2 <- o2
o2_2$A <- 0
res2 <- kbsd(data = o2, int_data_list = list(o2_1, o2_2), disthalf_vec=c(L1=10, L2=0.5, A=0.5*0.4),  # half of the SD for A!
             plot.out = F)
res2_plot <- kbsd(data = o2,
                  int_data_list = list(o2_1, o2_2),
                  disthalf_vec=c(L1=10, L2 = 0.5, A=0.2))
res2_plot
# visibly more support for IV=2 (A=0)
table(obs2$A)

# expected viol (confirmed by PoRT): age > 60, fit = 1, fit = 0

# what strata are those with low EDP in IV = 1 (A=1)
shift1 <- res2[res2$shift == 1,]
outliers1 <- shift1$diagnostic < quantile(res2[res2$shift==1, "diagnostic"], probs = 0.25)  # outliers (below whiskers' ends)
l_values1 <- obs2[outliers1, "L1"] 
diag_values1 <- shift1[outliers1,]
plot(l_values1, diag_values1$diagnostic)
# most are v young -> makes sense bc young ppl rarely treated but probs not extreme enough (only for a=0.01) relevant,
# but also some v old ppl! but NB: 4 outliers with high EDP?
l_values1 <- obs2[outliers1, "L2"]
diag_values1 <- shift1[outliers1,]
plot(l_values1, diag_values1$diagnostic)
# most with few support for A=1 have L2=1, i.e. are fit -> makes sense bc when trying to
# intervene fit ppl on A=1, only few other obs as fit ppl rarely treated

# what strata are those with low EDP in IV = 2 (A=0)
shift2 <- res2[res2$shift == 2,]
outliers2 <- shift2$diagnostic < quantile(res2[res2$shift==2, "diagnostic"], probs = 0.05)
l_values2 <- obs2[outliers2, "L1"] 
diag_values2 <- shift2[outliers2,]
plot(l_values2, diag_values2$diagnostic)
# most with few EDP are v old pp (but also some v young ones!) -> sensible bc few old ppl
# that got A=0 so that few support for such obs if IV02 (A00)
l_values2 <- obs2[outliers2, "L2"] 
diag_values2 <- shift2[outliers2,]
plot(l_values2, diag_values2$diagnostic)
# most with few support for A=0 are those with L2=0, i.e. unfit ppl bc unfit ppl are
# usually treated so that if intervene on A=0 for them, few unfit obs in neighbourhood


# 3 Confounders Version 1 ----
set.seed(15082025)
sem3 <- DAG.empty() +
  node("L1", distr = "rbern", prob = 0.3) +
  node("L2", distr = "rbern", prob = 0.1) + 
  node("L3", distr = "rbern", prob = 0.6) +
  node("A", distr = "rbern", prob = 0.26*L1 + 0.31*L2 + 0.42*L3)  # if L1, L2, L3=1 -> prob A=1
dag3 <- set.DAG(sem3)
plotDAG(dag3)
obs3 <- sim(dag3, rndseed = 12082025, n = 1000)
table(obs3$A)  

# check through all strata
# gamma = 1
table(obs3$L1, obs3$A)/rowSums(table(obs3$L1, obs3$A))
table(obs3$L2, obs3$A)/rowSums(table(obs3$L2, obs3$A))
table(obs3$L3, obs3$A)/rowSums(table(obs3$L3, obs3$A))  # low P(A=1|L3=0) -> viol #1 for b=0.1 (sample large enough with 386/1000)
# gamma =2
table(obs3[obs3$L1 == 0 & obs3$L2 == 0,"A"])
table(obs3[obs3$L1 == 0 & obs3$L3 == 0,"A"])  # expected viol #2: P(A=1)~0 for b=gruber (sample large enough with 27.8%)
table(obs3[obs3$L2 == 0 & obs3$L3 == 0,"A"])  # expected viol #3: P(A=1)~0 for b=0.1 (sample large enough with 34.6%)

table(obs3[obs3$L1 == 0 & obs3$L2 == 1,"A"])
table(obs3[obs3$L1 == 1 & obs3$L2 == 0,"A"])

table(obs3[obs3$L1 == 0 & obs3$L3 == 1,"A"])
table(obs3[obs3$L1 == 1 & obs3$L3 == 0,"A"])

table(obs3[obs3$L2 == 0 & obs3$L3 == 1,"A"])
table(obs3[obs3$L2 == 1 & obs3$L3 == 0,"A"])

table(obs3[obs3$L1 == 1 & obs3$L2 == 1,"A"])
table(obs3[obs3$L1 == 1 & obs3$L3 == 1,"A"])
table(obs3[obs3$L2 == 1 & obs3$L3 == 1,"A"])
# gamma = 3
obs3 %>% filter(L1==1 & L2==1 & L3==1 & A==1) %>% nrow()/
  obs3 %>% filter(L1==1 & L2==1 & L3==1) %>% nrow()  # expected viol #4: P(A=1)~1 for beta = 0.05/0.1



# PoRT ---
source('data/port_utils.R')
a_values <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.1)
gruber3 <- 5/(sqrt(nrow(obs3))*log(nrow(obs3)))
b_values <- c(0.01, gruber3, 0.05, 0.1)
g_values <- 1:3
lst3 <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst3[[paste0("gamma = ", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A", cov.quanti = NULL, cov.quali = c("L3", "L2", "L1"),
             data = obs3, alpha = a, beta = b, gamma = g) 
    }
  }
}
lst3
# gamma = 1:
# only for b = 0.1 the one wanted subgroup with L3=0 was detected; makes sense bc proba.exposure=0.073 which is only <0.1, not <0.05 etc.
# gamma = 2:
# all 3 viol detected for appropriate beta values
# BEST COMBO: any alpha, b = 0.05 to find all 3 viol (beta=gruber/0.1 only find 1 critical subgroup each)
# gamma = 3:
# for a=0.01, viol #4 with L1==1 & L2==1 & L3==1 -> P(A=1)~1 should have been found
# but instead, its complement with L1=0&L2=0&L3=0 -> P(A=1)~0 was found for b=0.01!
# reason: only for tree with split by 3 vars, a stratum with such an extreme beta could be found
#         proba.exposure not extreme enough for split by 1 var/by 2 vars, so continued to build tree

# essence: all violations except #4 were found, but #4 was found indirectly (L3=0 & L1=0 & L2=0) -> check if really so scarce
obs3 %>% filter(L3==0 & L1==0 & L2==0 & A==1) %>% nrow()/obs3 %>% filter(L3==0 & L1==0 & L2==0) %>% nrow() # it is


# KBSD ---
source("kbsd.R")
o3 <- obs3[-1]
o3_1 <- o3
o3_1$A <- 1
o3_2 <- o3
o3_2$A <- 0
res3 <- kbsd(data = o3,
             int_data_list = list(o3_1, o3_2),
             disthalf_vec=c(L1=0.5, L2 = 0.3, L3 = 0.5, A=0.5*0.5),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
res3_plot <- kbsd(data = o3,
                  int_data_list = list(o3_1, o3_2),
                  disthalf_vec=c(L1=0.5, L2 = 0.3, L3 = 0.5, A=0.5*0.5))  # use 1 SD for L_i, 0.5 SD for A
res3_plot

# what strata are those with low support for treated (IV = 1)
# according to port, the subgroups with low probs for treatment are: L3=0 & L1=0, L3=0 & L2=0, L3=0
subset_5 <- res3[res3$diagnostic < 50 & res3$shift == 1,]
obs3[subset_5$observation, ]
table(obs3[subset_5$observation, ][, c("L1", "L2")]) # many L1=0 & L2=0 among treated, but among L1=0 & L2=0 still enough treated&untreated (proba.exposure not extreme)
table(obs3[subset_5$observation, ][, c("L1", "L3")]) # confirmed that many from L3=0 & L1=0, ACTUALLY from L3=0
table(obs3[subset_5$observation, ][, c("L2", "L3")]) # confirmed that many from L3=0 & L2=0

# check additional group that port returned: how many among L3=0 & L1=0 & L2=0 treated/untreated:
table(obs3[obs3$L3==0 & obs3$L1==0 & obs3$L2==0, "A"])  # no one treated -> valid 
table(obs3[subset_5$observation, ][, "L1"],
      obs3[subset_5$observation, ][, "L2"],
      obs3[subset_5$observation, ][, "L3"]) # kbsd shows that they most small EDP from there -> det by kbsd as det by port!


# what strata are those with low support for untreated (IV = 2)
subset_6 <- res3[res3$diagnostic < 50 & res3$shift == 2,]
obs3[subset_6$observation, ]
table(obs3[subset_6$observation, ][, "L1"],
      obs3[subset_6$observation, ][, "L2"],
      obs3[subset_6$observation, ][, "L3"])  # split by L3: most of those with few EDP are L3=0 & L1=1 & L2=0 (for threshold=20, most are L3=1 & L1=0 & L2=1)
# check among L3=0 & L1=1 & L2=0, how many treated/untreated
table(obs3[obs3$L3==0 & obs3$L1==1 & obs3$L2==0, "A"])  # not extreme enough: proba.exposure = 0.17
# check among L3=1 & L1=0 & L2=1, how many treated/untreated
table(obs3[obs3$L1 == 0 & obs3$L2 == 1 & obs3$L3 == 1, "A"])  # not extreme enough: proba.exposure = 0.75
table(o3$A)  # also: way more untreated compared to treated in general!

# but generally, most EDP between 50 & 300 for both IV rules which is good compared to ex.2 

# essence: KBSD found all 3 strata among A=1 in which pos is violated (<=> low EDP), 
#          also found the additional stratum that PoRT had detected (L3=0 & L1=0 & L2=0),
#          even flagged 2 more strata among A=0 that do not present problem actually

### alternative metrics for disthalf_vec ----
# prob similar results as above
mad(o3$A)  # mad unsuitable
IQR(o3$A)  # using IQR would be suitable here, treatment distr not as imbalanced as in 2 Conf setting



# 3 Confounders Version 2 ----

# concatenate L1 + L2, so that low density in centre (bimodal distr) 
# → adjust A accordingly (sides with higher P(A), centre with lower P(A))
a <- rnorm(1000, 1, 1)
plot(density(a))
b <- rnorm(1000, 4, 1)
plot(density(b))
z <- c(a,b)
plot(density(z))
rug(z, col = "grey")

# first idea *******************************************************************
# L3_1 <- rnorm(500, 0, 1)
# L3_2 <- rnorm(500, 4, 1)
# L3_1_2 <- c(L3_1, L3_2)
# sem1 <- DAG.empty() +
#   node("L1", distr = "rbern", prob = 0.3) +
#   node("L2", distr = "rbern", prob = 0.1) +
#   node("L3", distr = "rconst", const = L3_1_2) +
#   #node("L4", distr = "rnorm", mean = 4, sd = 1) +
#   #node("L5", distr = "rnorm", mean = 5, sd = 1) +
#   #node("A", distr = "rbern",  prob = plogis(L2- 4*L1*((L3>=1)&(L3<=3))))
# dag1 <- set.DAG(sem1)
# plotDAG(dag1)
# data1 <- sim(dag1, n = 1000)
# ******************************************************************************

set.seed(27092025)
L3_1 <- rnorm(500, 3, 1)
L3_2 <- rnorm(500, 7, 1)
L3_1_2 <- c(L3_1, L3_2)
plot(density(L3_1_2))
sem1 <- DAG.empty() +
  node("L1", distr = "rbern", prob = 0.1) +
  node("L2", distr = "rbern", prob = 0.2) +
  node("L3", distr = "rconst", const = L3_1_2) +  # L3 follows a bimodal distribution,
  #        where vals around 5 have v low prob & vals to left & right have higher prob
  node("A", distr = "rbern", prob = plogis(L1 + L2 - 0.5*L3*(L3 < 6)))
dag1 <- set.DAG(sem1)
plotDAG(dag1)
data1 <- sim(dag1, n = 1000)

# visualise distribution of L3
plot(density(data1$L3), main = "bimodal L3 distribution from mixture")
# i.e. few values of L3 that are close to 5 -> will have viol there

data1 %>% filter(L3 < 4 & A==1) %>% nrow()/
  data1 %>% filter(L3 < 4) %>% nrow()

data1 %>% filter(L3 > 6 & A==1) %>% nrow()/
  data1 %>% filter(L3 > 6) %>% nrow()

data1 %>% filter(L3 >= 4 & L3 <= 6 & A==1) %>% nrow()/
  data1 %>% filter(L3 >= 4 & L3 <= 6) %>% nrow()  # viol #1: if L3 in [4,6], sample prop = 14.1% -> should find for g>=1, b=0.1

data1 %>% filter(L3 >= 4 & L3 <= 6 & L1==0 & A==1) %>% nrow()/
  data1 %>% filter(L3 >= 4 & L3 <= 6 & L1 == 0) %>% nrow()  # viol #2 (subviol of #1): if L3 in [4,6] AND L1=0, sample prop = 12.9% -> g>=2, b=0.1

data1 %>% filter(L3 >= 4 & L3 <= 6 & L2==0 & A==1) %>% nrow()/
  data1 %>% filter(L3 >= 4 & L3 <= 6 & L2 == 0) %>% nrow()  # viol #3 (subviol of #1): if L3 in [4,6] AND L2=0, sample prop = 11.3% -> g>=2, b=0.1


# distr of A within L1 and L2 is not too extreme
table(data1$L1, data1$A)
table(data1$L2, data1$A)


## PoRT: continuous var L3 uncategorised ----
source("data/port_utils.R")
lst5 <- list()
a_values <- c(0.01, 0.025, 0.05, 0.1)
b_values <- c(0.01, 5/(sqrt(nrow(dat))*log(nrow(data1))), 0.05, 0.1)
g_values <- 1:3
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst5[[paste0("gamma=", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A", cov.quanti = c("L3"), cov.quali = c("L1", "L2"),
             data = data1, alpha = a, beta = b, gamma = g)
    }
  }
}
lst5
#sink("port_3_bimodal_uncat.txt")
# gamma =1: detected the one subgroup to find (L3< 6.004 & L3>=3.956),
#           which was split into small strata for smaller a, b bc continuous
# gamma =2: viol with combo of L2=0 and L3< 5.626 & L3>=5.199 only for a=0.025 & b=0.01/gruber bc smaller stratum for L3 so that proba.exp extremer there
#           this viol #2 is covered by viol #1 for larger a=0.05/0.1 & b=0.1 (bc only detectable for b=0.1!)
# gamma =3: always identified where appropriate -> but one wrong subgroup? see below
data1 %>% filter( L3< 3.423 & L3>=3.169 & A==1) %>% nrow()/
  data1 %>% filter( L3< 3.423 & L3>=3.169) %>% nrow()  # error? subgroup size & proba is wrong!?


## PoRT: continuous var L3 categorised ----
data1_cat <- data1
source("data/port_utils.R")
data1_cat$L3 <- cut(data1_cat$L3, breaks = c(-Inf, 2, 4, 6, 8, Inf))
lst5_cat <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst5_cat[[paste0("gamma=", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A", cov.quanti = NULL, cov.quali = c("L1", "L2", "L3"),
             data = data1_cat, alpha = a, beta = b, gamma = g)
    }
  }
}
lst5_cat
#sink("port_3_bimodal_cat.txt")
# gamma = 1, 2, 3: det viol #1 but best example for info loss that can occur due to categorisation!!
#                  stratum where intersection with L2=0 (L2=0 & L3< 5.626 & L3>=5.199)
#                  is not included anymore bc categorisation too broad so that extreme
#                  proba that is in L2=0 & (L3< 5.626 & L3>=5.199), does not exist in (4,6]&L2=0 anymore


# essence: in both categorised & uncat both viol are found, but without categorisation,
#          can return preciser subgroup that includes L2=0, whereas after categorisation
#          the broader subgroup only is returned -> so emphasises importance of cat thresholds
#          appropriate/meaningful for the context as Chatton et al also emphasised -> to limit info loss


# KBSD ---
source("kbsd.R")
set.seed(27092025)
o5 <- data1[-1]
o5_1 <- o5
o5_1$A <- 1
o5_2 <- o5
o5_2$A <- 0
res5 <- kbsd(data = o5,
             int_data_list = list(o5_1, o5_2),
             disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3=sd(o5$L3), A=0.5*sd(o5$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
res5_plot <- kbsd(data = o5,
                  int_data_list = list(o5_1, o5_2),
                  disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3=sd(o5$L3), A=0.5*sd(o5$A)))  # use 1 SD for L_i, 0.5 SD for A
res5_plot
table(data1$A)  # see that alr fewer support among IV=1 bc fewer obs in A=1

# acc to viol, few support for IV=1 (A=1) if would estimate Y|A=1 further, for subgroup L3=(4,6]:
shift1 <- res5[res5$shift == 1,]
outliers1 <- shift1$diagnostic < quantile(shift1$diagnostic, probs = 0.25)  # create indices for the "outliers"
l_values1 <- data1[outliers1, "L3"]  # to which original obs (L values) do these outliers belong?
diag_values1 <- shift1[outliers1,] # what diag values do these outliers have
plot(l_values1, diag_values1$diagnostic)

# interesting: for threshold =0.05, expected those in [4,6] to have lowest EDP/support, 
#              but they're almost judged to have most
#              for threshold =0.25 also quite a few in [4,6] with high EDP & why this pattern in the plot?

l_values1 <- data1[outliers1, "L1"]
diag_values1 <- shift1[outliers1,]
plot(l_values1, diag_values1$diagnostic)  # clear trend for L1 tho: those with few support for IV=1 mostly have L1=1 

l_values1 <- data1[outliers1, "L2"]
diag_values1 <- shift1[outliers1,]
plot(l_values1, diag_values1$diagnostic)  # more scattered here, so few EDP mostly for L2=0

# also check for IV=2
shift2 <- res5[res5$shift == 2,]
outliers2 <- shift2$diagnostic < quantile(shift2$diagnostic, probs = 0.05)  # create indices for the "outliers"
l_values2 <- data1[outliers2, "L3"]  # original L3 values
diag_values2 <- shift1[outliers2,] # diag values
plot(l_values2, diag_values2$diagnostic)
# no clear pattern: obv less support for tails bc fewer obs there (did not have viol with P(A=0)~0)


# essence: seems like kbsd cannot PINPOINT that for IV=1 there is low support for L3 =[4,6], i.e. the "hole"
#          whereas PoRT could



# 5 Confounders ----

set.seed(28092025)
L3_1 <- rnorm(500, 3, 1)
L3_2 <- rnorm(500, 7, 1)
L3_1_2 <- c(L3_1, L3_2)
plot(density(L3_1_2))
sem1 <- DAG.empty() +
  node("L1", distr = "rbern", prob = 0.1) +
  node("L2", distr = "rbern", prob = 0.2) +
  node("L3", distr = "rconst", const = L3_1_2) +  # L3 follows a bimodal distribution,
  #        where vals around 5 have v low prob & vals to left & right have higher prob
  node("L4", distr = "rbern", prob = 0.1) +
  node("L5", distr = "rbern", prob = 0.2) +
  node("A", distr = "rbern", prob = plogis(0.2*L1 + 0.3*L2 + 0.1*L4 + 0.3*L5 - 0.5*L3*(L3 < 6)))
# all conf L_i=1 is v unlikely -> treatment not too likely, too but not extreme either
dag1 <- set.DAG(sem1)
plotDAG(dag1)
data1 <- sim(dag1, n = 1000)

plot(density(data1$L3), main = "bimodal L3 distribution from mixture") # again expect viol in [4,6]

# never any extreme P(A) for L1, L2, L4, L5 and their intersections, but for when L3 is involved:

data1 %>% filter(L3 >= 4 & L3 <= 6 & A==1) %>% nrow()/
  data1 %>% filter(L3 >= 4 & L3 <= 6) %>% nrow()  # viol #1: if L3 in [4,6], sample prop = 15.2% -> should find for g>=1, b=0.1

data1 %>% filter(L3 >= 4 & L3 <= 6 & L1==0 & A==1) %>% nrow()/
  data1 %>% filter(L3 >= 4 & L3 <= 6 & L1 == 0) %>% nrow()  # viol #2 (subviol of #1): if L3 in [4,6] AND L1=0, sample prop = 13.5% -> g>=2, b=0.1

data1 %>% filter(L3 >= 4 & L3 <= 6 & L2==0 & A==1) %>% nrow()/
  data1 %>% filter(L3 >= 4 & L3 <= 6 & L2 == 0) %>% nrow()  # viol #3 (subviol of #1): if L3 in [4,6] AND L2=0, sample prop = 11.9% -> g>=2, b=0.1

## PoRT: continuous var L3 uncategorised ----
source("data/port_utils.R")
lst5 <- list()
a_values <- c(0.01, 0.025, 0.05, 0.1)
b_values <- c(0.01, 5/(sqrt(nrow(dat))*log(nrow(data1))), 0.05, 0.1)
g_values <- 1:5
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst5[[paste0("gamma=", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A", cov.quanti = c("L3"), cov.quali = c("L1", "L2", "L4", "L5"),
             data = data1, alpha = a, beta = b, gamma = g)
    }
  }
}
lst5
# gamma = 1: detected , most informative for a=0.05/0.1 & b=0.1 (not too precise like for small a)
# gamma = 2-5: viol #1 detected, 2+3 only for a=0.05, b=gruber, else only covered


## PoRT: continuous var L3 categorised ----
data1_cat <- data1
source("data/port_utils.R")
data1_cat$L3 <- cut(data1_cat$L3, breaks = c(-Inf, 2, 4, 6, 8, Inf))
lst5_cat <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst5_cat[[paste0("gamma=", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A", cov.quanti = NULL, cov.quali = c("L1", "L2", "L3", "L4", "L5"),
             data = data1_cat, alpha = a, beta = b, gamma = g)
    }
  }
}
lst5_cat
# gamma = 1: viol #1 detected for b=0.1 as should
# gamma = 2: viol #1 again for b=0.1, and for a =0.01, new subgroup "L3=(4,6] & L4=1" but v small again
# gamma = 3-5: due to cat, L3 groups are larger, i.e. as contrary to uncategorised case,
#              where many small viol groups with L3 could be split, L3 is now large to that
#              viol subgroups are created by intersecting with e.g. 2 other vars! imp obs
#              viol #1 always det for b=0.1
# but so general trend is that for uncategorised case, small alpha means in case of 
# an existent viol subgroup, that subgroup is split into smaller strata to adhere to alpha
# -> for cat case, this is not poss anymore so that finding strata small enough that match the 
# small alpha is achieved by intersecting the categorised var with other vars (only poss for g>1 then tho)
data1_cat %>% filter(L3=="(4,6]" & L1==0 & L2==0 & L5==0 ) %>% nrow()


# KBSD ---
source("kbsd.R")
set.seed(28092025)
o5 <- data1[-1]
o5_1 <- o5
o5_1$A <- 1
o5_2 <- o5
o5_2$A <- 0
res5 <- kbsd(data = o5,
             int_data_list = list(o5_1, o5_2),
             disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3=sd(o5$L3), 
                            L4 = sd(o5$L4), L5 = sd(o5$L5), A=0.5*sd(o5$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
res5_plot <- kbsd(data = o5,
                  int_data_list = list(o5_1, o5_2),
                  disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3=sd(o5$L3),
                                 L4 = sd(o5$L4), L5 = sd(o5$L5), A=0.5*sd(o5$A)))  # use 1 SD for L_i, 0.5 SD for A
res5_plot  # again fewer support for IV=1
table(data1$A)  # which alr indicated here by fewer obs in A=1

# acc to viol, few support for IV=1 (A=1) if would estimate Y|A=1 further, for subgroup L3=(4,6]:
shift1 <- res5[res5$shift == 1,]
outliers1 <- shift1$diagnostic < quantile(shift1$diagnostic, probs = .25)  # create indices for the "outliers"
l_values1 <- data1[outliers1, "L3"]  # to which original obs (L values) do these outliers belong?
diag_values1 <- shift1[outliers1,] # what diag values do these outliers have
plot(l_values1, diag_values1$diagnostic)
# now clearer that L3 around 5 have fewer EDP than surrounding values! so problem of dimensionality or random?
# shows that if we intervened all on A=1, those with L3 around 5 have fewer support which reflects
# that there just don't exist many obs with L3~5 & A=1

# also check for IV=2 (A=0)
shift2 <- res5[res5$shift == 2,]
outliers2 <- shift2$diagnostic < quantile(shift2$diagnostic, probs = 0.05)  # create indices for the "outliers"
l_values2 <- data1[outliers2, "L3"]  # original L3 values
diag_values2 <- shift1[outliers2,] # diag values
plot(l_values2, diag_values2$diagnostic)
# indicates fewer support for L3 < 6, so check:
data1 %>% filter(L3 < 6 & A==0) %>% nrow()/data1 %>% filter(L3 < 6) %>% nrow()
# but actually many obs among A=0 with such values, also EDP threshold here v low
# so prob enough support overall (also did not have viol with P(A=0)~0)



# 10 Confounders ----

set.seed(29092025)
L3_1 <- rnorm(500, 3, 1)
L3_2 <- rnorm(500, 7, 1)
L3_1_2 <- c(L3_1, L3_2)
plot(density(L3_1_2))
sem1 <- DAG.empty() +
  node("L1", distr = "rbern", prob = 0.1) +
  node("L2", distr = "rbern", prob = 0.2) +
  node("L3", distr = "rconst", const = L3_1_2) +  # L3 follows a bimodal distribution,
  #        where vals around 5 have v low prob & vals to left & right have higher prob
  node("L4", distr = "rbern", prob = 0.1) +
  node("L5", distr = "rbern", prob = 0.2) +
  node("L6", distr = "rbern", prob = 0.1) +
  node("L7", distr = "rbern", prob = 0.2) +
  node("L8", distr = "rbern", prob = 0.1) +
  node("L9", distr = "rbern", prob = 0.2) +
  node("L10", distr = "rbern", prob = 0.1) +
  node("A", distr = "rbern", prob = plogis(0.2*L1 + 0.3*L2 + 0.1*L4 + 0.3*L5 - 0.5*L3*(L3 < 6)))
# all conf L_i=1 is v unlikely -> treatment not too likely, too but not extreme either
dag1 <- set.DAG(sem1)
plotDAG(dag1)
data1 <- sim(dag1, n = 1000)

data1 %>% filter(L3> 4 & L3 <6 & A==1) %>% nrow()/data1 %>% filter(L3 <6) %>% nrow()  # viol for g>=1, b>=0.05, sample prop =57.4%


## PoRT: continuous var L3 uncategorised ----
source("data/port_utils.R")
lst5 <- list()
a_values <- c(0.01, 0.025, 0.05, 0.1)
b_values <- c(0.01, 5/(sqrt(nrow(dat))*log(nrow(data1))), 0.05, 0.1)
g_values <- 1:10
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst5[[paste0("gamma=", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A", cov.quanti = c("L3"),
             cov.quali = c("L1", "L2", "L4", "L5", "L6", "L7", "L8", "L9", "L10"),
             data = data1, alpha = a, beta = b, gamma = g)
    }
  }
}
lst5
# gamma = 1:
# gamma = 2-10: 


## PoRT: continuous var L3 categorised ----
data1_cat <- data1
source("data/port_utils.R")
data1_cat$L3 <- cut(data1_cat$L3, breaks = c(-Inf, 2, 4, 6, 8, Inf))
lst5_cat <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst5_cat[[paste0("gamma=", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A", cov.quanti = NULL,
             cov.quali = c("L1", "L2", "L4", "L5", "L6", "L7", "L8", "L9", "L10"),
             data = data1, alpha = a, beta = b, gamma = g)
    }
  }
}
lst5_cat
# gamma = 1:
# gamma = 2: 
# gamma = 3-10:


# KBSD ---
source("kbsd.R")
set.seed(20092025)
o5 <- data1[-1]
o5_1 <- o5
o5_1$A <- 1
o5_2 <- o5
o5_2$A <- 0
res5 <- kbsd(data = o5,
             int_data_list = list(o5_1, o5_2),
             disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3=sd(o5$L3), 
                            L4 = sd(o5$L4), L5 = sd(o5$L5), L6 = sd(o5$L6),
                            L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9),
                            L10 = sd(o5$L10), A=0.5*sd(o5$A)), plot.out = F)
res5_plot <- kbsd(data = o5,
                  int_data_list = list(o5_1, o5_2),
                  disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3=sd(o5$L3),
                                 L4 = sd(o5$L4), L5 = sd(o5$L5), L6 = sd(o5$L6),
                                 L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9),
                                 L10 = sd(o5$L10), A=0.5*sd(o5$A)))  # use 1 SD for L_i, 0.5 SD for A
res5_plot  # again fewer support for IV=1 & NB: lower EDP overall bc more dims <=> more diff for obs to be close
table(data1$A)  # which alr indicated here by fewer obs in A=1


# acc to viol, few support for IV=1 (A=1) if would estimate Y|A=1 further, for subgroup L3=(4,6]:
shift1 <- res5[res5$shift == 1,]
outliers1 <- shift1$diagnostic < quantile(shift1$diagnostic, probs = .25)  # create indices for the "outliers"
l_values1 <- data1[outliers1, "L3"]  # to which original obs (L values) do these outliers belong?
diag_values1 <- shift1[outliers1,] # what diag values do these outliers have
plot(l_values1, diag_values1$diagnostic)
# now clearer that L3 around 5 have fewer EDP than surrounding values! so problem of dimensionality or random?
# shows that if we intervened all on A=1, those with L3 around 5 have fewer support which reflects
# that there just don't exist many obs with L3~5 & A=1

# also check for IV=2 (A=0)
shift2 <- res5[res5$shift == 2,]
outliers2 <- shift2$diagnostic < quantile(shift2$diagnostic, probs = 0.05)  # create indices for the "outliers"
l_values2 <- data1[outliers2, "L3"]  # original L3 values
diag_values2 <- shift1[outliers2,] # diag values
plot(l_values2, diag_values2$diagnostic)

### formula for EDP for high-dim covar set ----

# type = "minval" instead of default type = "Rfast"
res5_mv_plot <- kbsd(data = o5, int_data_list = list(o5_1, o5_2), type = "minval",
                     minval_vec = rep(0.5, 11),
                     disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3 = sd(o5$L3), L4 = sd(o5$L4), L5 = sd(o5$L5),
                                    L6 = sd(o5$L6), L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9), L10 = sd(o5$L10),
                                    A=0.5*sd(o5$A)), plot.out = T)
#ggsave("kbsd_10_mv.png")
res5_mv <- kbsd(data = o5, int_data_list = list(o5_1, o5_2), type = "minval",
                minval_vec = rep(0.5, 11),
                disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3 = sd(o5$L3), L4 = sd(o5$L4), L5 = sd(o5$L5),
                               L6 = sd(o5$L6), L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9), L10 = sd(o5$L10),
                               A=0.5*sd(o5$A)), plot.out = F)

# type = "harmonicmean" instead of default type = "Rfast"
res5_hm_plot <- kbsd(data = o5, int_data_list = list(o5_1, o5_2), type = "harmonicmean",
                     disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3 = sd(o5$L3), L4 = sd(o5$L4), L5 = sd(o5$L5),
                                    L6 = sd(o5$L6), L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9), L10 = sd(o5$L10),
                                    A=0.5*sd(o5$A)), plot.out = T)
#ggsave("kbsd_10_hm.png", width = 6, height = 3)
res5_hm <- kbsd(data = o5, int_data_list = list(o5_1, o5_2), type = "harmonicmean",
                disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3 = sd(o5$L3), L4 = sd(o5$L4), L5 = sd(o5$L5),
                               L6 = sd(o5$L6), L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9), L10 = sd(o5$L10),
                               A=0.5*sd(o5$A)), plot.out = F)

