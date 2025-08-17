
# Simulating causal settings with 1,2 and 3 confounders with sparse data for some strata
# and checking if different diagnostics detect these strata with positivity violations.
library(simcausal)
set.seed(15082025)

## 1) One Confounder ----
# Baseline covariate: Health score (L) ~ Normal(0, 1) -> neg = bad, pos = healthy
# Treatment: Vaccination (A) ~ Ber with logit depending on L
# Outcome: Hospitalisation (Y) ~ Ber depending on A and L

# simulate so that neg health score <=> high probs of treatment, pos health score <=> low probs of treatment
# expected viol: P(A=0|neg health)~0 or P(A=1|neg health)~1, and P(A=0|pos health)~1 or P(A=1|pos health)~0

sem1 <- DAG.empty() +
  node("L", distr = "rnorm", mean = 0, sd = 1) +
  node("A", dist = "rbern", prob = ifelse(L<0, 1-0.1*plogis(L), 0.1*plogis(L)))
dag1 <- set.DAG(sem1)
obs1 <- sim(dag1, rndseed = 05082025, n = 1000)

# categorise Health Score L to prevent greedy categorisation (Schomaker & Chatton, Danelian & Chatton):
min(obs1$L)
max(obs1$L)
obs1$L <- cut(obs1$L, breaks = c(min(obs1$L), 0, max(obs1$L)))  # binary categorisation
obs1$L <- cut(obs1$L, breaks = seq(-4,4,by = 1)) # finer categorisation


# PoRT ---
setwd("C:/Users/victo/OneDrive/Desktop/Uni/Statistik/BA_thesis")
source('data/port_utils.R')
lst1 <- list(port = NULL, port_risca = NULL)
a_values <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.1)
gruber <- 5/(sqrt(nrow(obs1))*log(nrow(obs1)))
b_values <- c(0.001, 0.01, gruber, 0.05, 0.1)
for (a in a_values) {
  for (b in b_values) {
    # without categorisation: L as cov.quanti
    #lst1$port[[paste0("alpha = ", a, ", beta = ", b)]] <- port(A = "A", cov.quanti = "L", cov.quali = NULL,
    #                                                           data = obs1, alpha = a, beta = b, gamma = 1)
    # with categorisation: L as cov.quali
    lst1$port[[paste0("alpha = ", a, ", beta = ", b)]] <- port(A = "A", cov.quali = "L", cov.quanti = NULL,
                                                               data = obs1, alpha = a, beta = b, gamma = 1) 
    #lst1$port_risca[[paste0("alpha value = ", a)]] <- RISCA::port(group = "A", cov.quanti = "L", cov.quali = NULL, alpha = a, beta = 0.05, data = obs1, gamma = 1)
  }
}
lst1
# without categorisation:
 # a = 0.01 & any b: v small subgroups esp for small b -> a value prob too small
 # a = 0.02/0.03: if b v small (until b = gruber), then subgroups still correct, from b = 0.05 extra 2 groups we did not want
 # a = 0.04: critical group with pos health only detected from b >= gruber
 # a = 0.05/0.1: " " only detected from b >= 0.5
# with binary categorisation into neg & pos:
 # identified violation, but only for beta >= gruber (if too small, P(A=0/1) not sufficiently extreme)
 # crticial subgroup of healthy and treated only identified for beta = 0.1 (probs not too extreme)
# with finer categorisation: similar as for binary

# essence: PoRT detected both subgroups that we wanted to be detected, but only when considering beta=0.1


# RISCA::port (only works with L as conv.quanti, i.e. not categorised) ---
# a = 0.05 (default), b = 0.05 (default), g = 1
RISCA::port(group = "A", cov.quanti = "L", cov.quali = NULL, data = obs1, gamma = 1)
# if health status pos, then P(A=1) v low
# if health status neg, then P(A=1) v high -> P(A=0) v unlikely (as for PoRT)

# a = 0.05 (default), b = gruber ~ 0.023, g = 1
gruber <- 5/(sqrt(nrow(obs1))*log(nrow(obs1)))
RISCA::port(group = "A", cov.quanti = "L", cov.quali = NULL, data = obs1, alpha = 0.05, beta = gruber, gamma = 1)
# err? first subgroup: 4.4% > b but should be beta < ß
# err? second subgroup: 4.7% < a but should be alpha >= a
# err? third: 3.6% > b
# err? fourth: 3.4% < a and 2.9% > b
# err? fifth: 2.3% < a

# !!!: do both cond have to be met? both alpha > a & beta < ß ?
#      but why first subgroup has alpha < a & beta > ß?

# a = 0.01, b = 0.05 (default), g = 1
RISCA::port(group = "A", cov.quanti = "L", cov.quali = NULL, data = obs1, alpha = 0.01, beta = 0.05, gamma = 1)
# covers subgroup form first case, just split into 2 "preciser" groups bc a = 0.01 allows for smaller groups
# but could've shown same large group as in first case or first case could've shown these 2 already

# essence: RISCA::port detected subgroup that we wanted to be detected, but also returned 2 other subgroups for beta = gruber


# KBSD ---
source("kbsd.R")

# need to retract categorisation:
obs1 <- sim(dag1, rndseed = 05082025, n = 1000)

o1 <- obs1[-1]
o1_1 <- o1
o1_1["A"] <- 1
o1_2 <- o1
o1_2["A"] <- 0
res1 <- kbsd(data = o1, int_data_list = list(o1_1, o1_2), disthalf_vec=c(L=1, A=0.5), plot.out = F)
res1_plot <- kbsd(data = o1,
                  int_data_list = list(o1_1, o1_2),
                  disthalf_vec=c(L=1, A=0.5))
res1_plot
# obs: patients treated have more support than patients not treated
sum(o1$L > 0)
sum(o1$L < 0)  # possible explanation: among treated, more have neg health -> more with neg health in general?

# but overall, most obs have EDP between 50 - 200
# question is: who are those among treated with low EDP (prob those with good health -> only few of such)
#              -> pos concern: problematic to assess what would happen to them if they do take vaccination
#              who are those among non-treated with low EDP (prob those with bad health and e.g. contraindication to A)
#              -> pos concern: problematic to assess what would happen if they were treated


# check what strata/confounder values those with low EDP in IV = 1 represent (prob those with pos health score)
subset1 <- res1[res1$shift == 1,]
subset1 <- subset1$diagnostic < 100#20 for outliers
obs1[subset1,]
plot(seq_len(nrow(obs1[subset1,])), obs1[subset1,"L"])   # more pos values as expected: healthy ppl in IV=1 have few EDP

# check what strata/confounder values those with low EDP in IV = 2 represent 
subset2 <- res1[res1$shift == 2,]
subset2 <- subset2$diagnostic < 100
obs1[subset2,]
plot(seq_len(nrow(obs1[subset2,])), obs1[subset2,"L"])
  # obs with low EDP among A=0 (not treated): neg health (rare that those with bad health not treated)

# essence: KBSD identified the few healthy&treated (subset1) and nonhealthy&non-treated (subset2) -> matches with port




## 2) Two Confounders ----
# Baseline confounders: Age (L1) ~ N(50, 10) and Fitness (L2) ~ Ber(0.05) -> most people in sample are not fit
# Treatment: BP medication (A) ~ Ber with logit(Age, Fitness)
# Cont. Outcome: Systolic BP at follow-up (Y) ~ (Age, Treatment)

# simulate obs data so that for stratum age>75 and fit==0 most treated -> few non-treated
sem2 <- DAG.empty() +
  node("age", distr = "rnorm", mean = 50, sd = 10) +
  #node("fit", distr = "rbern", prob = 0.1) +
  node("fit", distr = "rbern", prob = ifelse(age > 75, 0.15, 0.65)) +  # fitness depends on age
  # if old & unfit -> def treated, i.e. P(A=0|old & unfit) should be low
  node("A", distr = "rbern", prob = ifelse((age > 75 & fit == 0), 0.95, 1-0.95))
dag2 <- set.DAG(sem2)
plotDAG(dag2)
obs2 <- sim(dag2, rndseed = 30072025, n = 1000)

# check how many observations with age > 75
obs2[obs2$age > 75,]  # as wanted: only 1 old & unfit that was not treated
  # expected viol #1: P(A=0|age >75 & fit == 0) = 0 or P(A=1|age >75 & fit == 0) = 1
nrow(obs2[obs2$age <= 75 & obs2$A == 1,])  # young rarely treated
  # expected viol #2: P(A=0|age <= 75)~1 or P(A=1|age <= 75)=46/990 = 0.04
nrow(obs2[obs2$fit == 1 & obs2$A == 1, ])
  # expected viol #3: P(A=1|fit == 1) ~0 since only 35/613 = 0.057
nrow(obs2[obs2$age <= 75 & obs2$fit == 1 & obs2$A == 1,])  # young & fit rarely treated, only 35/612= 0.057
  # expected viol #4 (intersection of the 2 above)

table(cut(obs2$age, breaks = c(0,75,90)), obs2$fit)
table(cut(obs2$age, breaks = c(0,20,40,75,90)), obs2$fit)
  # most young (age <= 75) are fit, almost all old are unfit
table(cut(obs2$age, breaks = c(0,20,40,75,90)), obs2$A)
  # most young (age <= 75) not treated, almost all old treated
table(obs2$fit, obs2$A, cut(obs2$age, breaks = c(0,75,90)))
  # only always treated if old & unfit -> should be more treated if young & unfit: include  "| (age <= 75 & fit == 0)" above?

# categorise age
obs2$age <- cut(obs2$age, breaks = c(0,75,90))
obs2$age <- cut(obs2$age, breaks = c(0,20,40,75,90))


# PoRT ---
source('data/port_utils.R')
lst2 <- list(port = NULL, port_risca = NULL)
a_values <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.1)
gruber2 <- 5/(sqrt(nrow(obs2))*log(nrow(obs2)))
b_values <- c(0.001, 0.01, gruber2, 0.05, 0.1)
for (a in a_values) {
  for (b in b_values) {
    # uncategorised age<s
    lst2$port[[paste0("alpha = ", a, ", beta = ", b)]] <- port(A = "A", cov.quanti = "age", cov.quali = "fit",
                                                               data = obs2, alpha = a, beta = b, gamma = 2) 
    # categorised age
    #lst2$port[[paste0("alpha = ", a, ", beta = ", b)]] <- port(A = "A", cov.quanti =NULL, cov.quali = c("age", "fit"),
    #                                                           data = obs2, alpha = a, beta = b, gamma = 2)
    #lst2$port_risca[[paste0("alpha value = ", a)]] <- RISCA::port(group = "A", cov.quanti = "L", cov.quali = NULL, alpha = a, beta = 0.05, data = obs1, gamma = 1)
  }
}
lst2
# without categorisation: identical for gamma = 1 and gamma = 2 bc only split by age always
  # only subgroup of age < 75 ever identified as "too rarely treated", never group of age>75 & fit==0
  # BUT: when are the violating subgroups returned and when "The whole sample presents at least one exposure modality's ..."?
  #      bc maybe the violating subgroup with fit is hidden behind the statement
# first categorisation: too imprecise, only group splitted by age (never fit) identified & only if beta >= 0.05
# second categorisation: also only identified if beta >= 0.05

# essence: all problematic subgroups found, but NB: group of age>75 & fit==0 always treated
#          but such a small subgroup that never detected by port (even for a=0.01 -> sample prop even smaller with 0.009)
#          -> would be negligeable for causality estimation as too little impact when so small


# RISCA::port ---
# a = 0.05 (default), b = 0.05 (default), g = 2 (default)
RISCA::port(group = "A", cov.quanti = c("age", "fit"), cov.quali = NULL, data = obs2,
            alpha = 0.05, beta = 0.05, gamma = 2)
  # as above: P(A=1|age < 75) ~ 0 and P(A=1|fit = 1) ~0
### NB: says "greater than or equal to 0 on variable 'fit'" but actually is only >0 !! ###

# a = 0.05 (default), b = gruber, g = 2 (all default)
RISCA::port(group = "A", cov.quanti = c("age", "fit"), cov.quali = NULL, data = obs2,
            alpha = 0.05, beta = gruber2, gamma = 2)
  # err? first subgroup here should not be included: 4.6% > gruber

# a = 0.01, b = 0.05, g = 2 (default)
RISCA::port(group = "A", cov.quanti = c("age", "fit"), cov.quali = NULL, data = obs2, alpha = 0.01, beta = 0.05, gamma = 2)
  # more subgroups bc defined more precisely since lower a permits them to be smaller -> all disjoint
  # but all 3 could've been returned for a = 0.05 too? subgroup sizes all > 0.05

# essence: all subgroups disjoint & detected that younger ppl/fit ppl rarely treated


# KBSD ---
source("kbsd.R")
# revert categorisation
obs2 <- sim(dag2, rndseed = 30072025, n = 1000)

o2 <- obs2[-1]
o2_1 <- o2
o2_1$A <- 1
o2_2 <- o2
o2_2$A <- 0
res2 <- kbsd(data = o2, int_data_list = list(o2_1, o2_2), disthalf_vec=c(age=10, fit=0.5, A=0.2), plot.out = F)
res2_plot <- kbsd(data = o2,
                  int_data_list = list(o2_1, o2_2),
                  disthalf_vec=c(age=10, fit = 0.5, A=0.2))
res2_plot  # what strata have low EDP? same for outliers among those treated: what covar values do they have?
table(o2$A)  # already imbalanced here: way more untreated compared to treated! 

# inspect what strata are those with low EDP in IV rule 1 (A=1)
subset_3 <- res2[res2$diagnostic < 50 & res2$shift == 1, ]  # all 1,000 obs for A=1 since all have diag < 50
sum(obs2[subset_3$observation,][, "age"] > 75)
sum(obs2[subset_3$observation,][, "age"] <= 75) # confirms that subgroup with low support for A=1 has more younger people
table(obs2[subset_3$observation,][, "fit"])  # confirms that those with low support for A=1 are mostly fit ones
# more younger ppl makes sense -> treated far less often
# fit ppl makes sense -> treated less often

# check what strata are those with low EDP in IV rule 2 (A=0)
subset4 <- res2[res2$diagnostic < 100 & res2$shift == 2, ]
table(obs2[subset_4$observation,][, "fit"], cut(obs2[subset_4$observation,][, "age"], breaks = c(0, 75, 90)))
  # most with low EDP are unfit (first row), the vast majority is young (first col: prob due to class imbalance for age)
  # biggest intersecting subgroup is unfit & young -> rarely not treated
# class imbalance: more younger ppl is prob side effect bc so much more young ppl overall
table(cut(obs2$age, breaks = c(0,75,90)))
# more unfit ppl makes sense: rare that unfit & not treated <=> low EDP

# essence: KBSD identified 3 of 4 groups that risk to violate pos: young& A=1, fit& A=1, unfit& A=0, young& A=0



## 3) Three Confounders ----
sem3 <- DAG.empty() +
  node("L1", distr = "rbern", prob = 0.3) +
  node("L2", distr = "rbern", prob = 0.1) + 
  node("L3", distr = "rbern", prob = 0.6) +
  node("A", distr = "rbern", prob = 0.2*L1 + 0.3*L2 + 0.4*L3)  # if L1, L2, L3=1 -> prob A=1
dag3 <- set.DAG(sem3)
plotDAG(dag3)
obs3 <- sim(dag3, rndseed = 12082025, n = 1000)

# check through all strata
# gamma = 1
table(obs3$L1, obs3$A)/rowSums(table(obs3$L1, obs3$A))
table(obs3$L2, obs3$A)/rowSums(table(obs3$L2, obs3$A))
table(obs3$L3, obs3$A)/rowSums(table(obs3$L3, obs3$A))  # low P(A=1|L3=0) -> expected viol #1
# gamma =2
table(obs3[obs3$L1 == 0 & obs3$L2 == 0,"A"])
table(obs3[obs3$L1 == 0 & obs3$L3 == 0,"A"])  # v imbalanced -> expected viol #2
table(obs3[obs3$L2 == 0 & obs3$L3 == 0,"A"])  # v imbalanced -> expected viol #3

table(obs3[obs3$L1 == 0 & obs3$L2 == 1,"A"])
table(obs3[obs3$L1 == 1 & obs3$L2 == 0,"A"])

table(obs3[obs3$L1 == 0 & obs3$L3 == 1,"A"])
table(obs3[obs3$L1 == 1 & obs3$L3 == 0,"A"])

table(obs3[obs3$L2 == 0 & obs3$L3 == 1,"A"])
table(obs3[obs3$L2 == 1 & obs3$L3 == 0,"A"])

table(obs3[obs3$L1 == 1 & obs3$L2 == 1,"A"])
table(obs3[obs3$L1 == 1 & obs3$L3 == 1,"A"])
table(obs3[obs3$L2 == 1 & obs3$L3 == 1,"A"])


# PoRT ---
source('data/port_utils.R')
lst3 <- list(port = NULL, port_risca = NULL)
a_values <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.1)
gruber3 <- 5/(sqrt(nrow(obs2))*log(nrow(obs3)))
b_values <- c(0.001, 0.01, gruber2, 0.05, 0.1)
for (a in a_values) {
  for (b in b_values) {
    lst3$port[[paste0("alpha = ", a, ", beta = ", b)]] <- port(A = "A", cov.quanti = NULL, cov.quali = c("L1", "L2", "L3"),
                                                               data = obs3, alpha = a, beta = b, gamma = 2) 
    #lst3$port_risca[[paste0("alpha value = ", a)]] <- RISCA::port(group = "A", cov.quanti = "L", cov.quali = NULL, alpha = a, beta = 0.05, data = obs1, gamma = 1)
  }
}
lst3
# gamma = 1:
  # only for b = 0.1 the one wanted subgroup with L3=0 was detected
# gamma = 2:
  # b = gruber/0.05/0.1 best, no matter a -> b = 0.05 even better as gruber/0.1 only find one subgroup each

# essence: all violations found! but requires optimising/trying out multiple HP values


# KBSD ---
o3 <- obs3[-1]
o3_1 <- o3
o3_1$A <- 1
o3_2 <- o3
o3_2$A <- 0
res3 <- kbsd(data = o3,
             int_data_list = list(o3_1, o3_2),
             disthalf_vec=c(L1=0.5, L2 = 0.3, L3 = 0.5, A=0.5),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
res3_plot <- kbsd(data = o3,
                  int_data_list = list(o3_1, o3_2),
                  disthalf_vec=c(L1=0.5, L2 = 0.3, L3 = 0.5, A=0.5))  # use 1 SD for L_i, 0.5 SD for A
res3_plot

# what strata are those with low support for treated (IV = 1)
# according to port, the subgroups with low probs for treatment are: L3=0 & L1=0, L3=0 & L2=0, L3=0
subset_5 <- res3[res3$diagnostic < 100 & res3$shift == 1,]
obs3[subset_5$observation, ]
table(obs3[subset_5$observation, ][, c("L1", "L2")]) 
table(obs3[subset_5$observation, ][, c("L1", "L3")]) # confirmed that many from L3=0 & L1=0
table(obs3[subset_5$observation, ][, c("L2", "L3")]) # confirmed that many from L3=0 & L2=0
table(obs3[subset_5$observation, ][, "L3"])  # many from L3=0

# what strata are those with low support for untreated (IV = 2)
subset_6 <- res3[res3$diagnostic < 20 & res3$shift == 2,]
obs3[subset_6$observation, ]
table(obs3[subset_6$observation, ][, "L1"],
      obs3[subset_6$observation, ][, "L2"],
      obs3[subset_6$observation, ][, "L3"])  # split by L3: most are L3=0 & L1=1 & L2=0
table(o3$A)  # also: way more untreated compared to treated in general!
  # for diag < 20: all from group L2==1

# but generally, most EDP between 50 & 300 for both IV rules which is good compared to ex.2 

# essence: KBSD found the strata in which pos is violated (<=> low EDP), even one more among
#          untreated (A=0 / shift=2) where all have L2=1 


