
# Simulating causal settings with 1,2 and 3 confounders with sparse data for some strata
# and checking if different diagnostics detect these strata with positivity violations.
library(simcausal)
library(tidyverse)
theme_set(theme_minimal())
set.seed(15082025)

## 1) One Confounder ----
# Baseline covariate: Health score (L) ~ Normal(0, 1) -> neg = bad, pos = healthy
# Treatment: Vaccination (A) ~ Ber with logit depending on L
# Outcome: Hospitalisation (Y) ~ Ber depending on A and L

# simulate so that neg health score <=> high probs of treatment, pos health score <=> low probs of treatment
# expected viol #1: P(A=1|neg health)~1
# expected viol #2: P(A=1|pos health)~0

sem1 <- DAG.empty() +
  node("L", distr = "rnorm", mean = 0, sd = 1) +
  node("A", dist = "rbern", prob = ifelse(L<0, 1-0.1*plogis(L), 0.1*plogis(L)))
  # if neg health, then v probably treated; if pos health, treat with v low probs
dag1 <- set.DAG(sem1)
obs1 <- sim(dag1, rndseed = 05082025, n = 1000)
plot(obs1[-1])

# categorise L to prevent greedy categorisation (Schomaker & Chatton, Danelian & Chatton):
min(obs1$L)
max(obs1$L)
obs1$L <- cut(obs1$L, breaks = c(min(obs1$L), 0, max(obs1$L)))  # binary categorisation
obs1$L <- cut(obs1$L, breaks = seq(-4,4,by = 1)) # finer categorisation


# PoRT ---
setwd("C:/Users/victo/OneDrive/Desktop/Uni/Statistik/BA_thesis/positivity-violations-diag")
source("data/port_utils.R")
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
  }
}
lst1
# without categorisation:
 # subgroup of healthy only from beta = 0.05, as P(A=1|L>0)=0.04 not that extreme
 # except for v small alpha, but too small as subgroups too precise, e.g. alpha = 0.01/0.02
 # BEST COMBO: a=0.05, b=0.05 or a=0.05, b=0.1 (= a=0.1,b=0.1) to be stricter with probs for pos viol
 # a = 0.1, b = 0.05 already misses 0 < L < 1, because P(A|0<L<1) = 0.07 -> only detected with b=0.1
# with binary categorisation into neg & pos health:
 # both viol only identified for beta = 0.1 as for healthy, P(A=1|L>0)=0.04 not that extreme
 # but also viol #1 only from beta >= gruber bc exposure not too extreme
 # BEST COMBO: a = 0.05, b = 0.1
# with finer categorisation:
 # similar as for binary, although pos group partly from beta=0.05 already thanks to finer categories

# essence: PoRT found both subgroups that we wanted to be detected, but needed beta=0.1 to include all


# KBSD ---
source("kbsd.R")

# need to retract categorisation:
obs1 <- sim(dag1, rndseed = 05082025, n = 1000)

o1 <- obs1[-1]
o1_1 <- o1
o1_1["A"] <- 1
o1_2 <- o1
o1_2["A"] <- 0
res1 <- kbsd(data = o1, int_data_list = list(o1_1, o1_2), disthalf_vec=c(L=1, A=0.5*0.5), # HALF of SD for A
             plot.out = F)
res1_plot <- kbsd(data = o1,
                  int_data_list = list(o1_1, o1_2),
                  disthalf_vec=c(L=1, A=0.5))
res1_plot
# obs: patients treated have more support than patients not treated

# but overall, most obs have EDP between 200-400
# question is: who are those among treated with low EDP (prob those with good health -> only few of such)
#              -> pos concern: problematic to assess what would happen to them if they do take vaccination
#              who are those among non-treated with low EDP (prob those with bad health and e.g. contraindication to A)
#              -> pos concern: problematic to assess what would happen if they were treated


# check what strata/conf values those with low EDP in IV = 1 represent (prob those with pos health score)
shift1 <- res1[res1$shift == 1,]
outliers1 <- shift1$diagnostic < quantile(res1[res1$shift==1, "diagnostic"], probs = 0.05)  # outliers (below whiskers' ends)
l_values1 <- obs1[outliers1, "L"]  # to which original obs (&L values) do these outliers belong?
diag_values1 <- shift1[outliers1,] # what diag values do these outliers have
plot(l_values1, diag_values1$diagnostic)
# support visibly decreases for more extreme L -> L~0 have a lot of support so not included in the plot
# also mostly positive L values as expected -> means few healthy ppl received IV=1 (confirms viol #2)

# check what strata/confounder values those with low EDP in IV = 2 (everyone received A=0) represent
shift2 <- res1[res1$shift == 2,]
outliers2 <- shift2$diagnostic < quantile(res1[res1$shift==2, "diagnostic"], probs = 0.05)
l_values2 <- obs1[outliers2, "L"]  # to which original obs (&L values) do these outliers belong?
diag_values2 <- shift2[outliers2,] # what diag values do these outliers have
plot(l_values2, diag_values2$diagnostic)
# the smaller L, the less support -> obs with negative L are rare in A=0 (rare that those with bad health not treated) -> confirms viol #1

# essence: KBSD identified both critical strata
# but ex. IV = 2: when checking for EDP threshold =100, not all neg are identified
# i.e. the higher the EDP threshold, the more matches with port results, the more is "covered"
# but kbsd, esp plotting L values against diag values v useful for first identification of critical strata!


# best case disthalf vec: doubled
res1 <- kbsd(data = o1, int_data_list = list(o1_1, o1_2), disthalf_vec=c(L=2*1, A=2*0.5*0.5),
             plot.out = F)
res1_plot <- kbsd(data = o1,
                  int_data_list = list(o1_1, o1_2),
                  disthalf_vec= c(L=2*1, A=2*0.5*0.5))
res1_plot

# worst case disthalf vec: halved
res1 <- kbsd(data = o1, int_data_list = list(o1_1, o1_2), disthalf_vec=c(L=0.5*1, A=0.5*0.5*0.5),
             plot.out = F)
res1_plot <- kbsd(data = o1,
                  int_data_list = list(o1_1, o1_2),
                  disthalf_vec= c(L=0.5*1, A=0.5*0.5*0.5))
res1_plot

# essence: always same trend across IV levels (generally better support for IV=1)
# plots for strata also yield same results for underlying L values


### alternative metrics for disthalf vec ----
# 1) MAD within L, A instead of SD
res1 <- kbsd(data = o1, int_data_list = list(o1_1, o1_2), disthalf_vec=c(L=mad(o1$L), A=0.5*mad(o1$A)), # HALF of SD for A
             plot.out = F)
res1_plot <- kbsd(data = o1,
                  int_data_list = list(o1_1, o1_2),
                  disthalf_vec=c(L=mad(o1$L), A=0.5*mad(o1$A)))
res1_plot
# unsuitable: assume that the error is due to mad(A) = 0 ?

# 2) IQR for L, A instead of SD
res1 <- kbsd(data = o1, int_data_list = list(o1_1, o1_2), disthalf_vec=c(L=IQR(o1$L), A=0.5*IQR(o1$A)), # HALF of SD for A
             plot.out = F)
res1_plot <- kbsd(data = o1,
                  int_data_list = list(o1_1, o1_2),
                  disthalf_vec=c(L=IQR(o1$L), A=0.5*IQR(o1$A)))
res1_plot

# 3) average pairwise distance within L, A
res1 <- kbsd(data = o1, int_data_list = list(o1_1, o1_2),
             disthalf_vec=c(L=mean(dist(matrix(o1$L), method = "euclidean")),
                            A=0.5*mean(dist(matrix(o1$A), method = "euclidean"))),
             plot.out = F)
res1_plot <- kbsd(data = o1,
                  int_data_list = list(o1_1, o1_2),
                  disthalf_vec=c(L=mean(dist(matrix(o1$L), method = "euclidean")),
                                 A=0.5*mean(dist(matrix(o1$A), method = "euclidean"))))
res1_plot

# essence: all have same trend, but good for robustness as less sensitive to outliers compared to SD
# plots for strata also yield same results for underlying L values



## 2) Two Confounders ----
# Baseline confounders: Age (L1) ~ N(50, 10) and Fitness (L2) ~ Ber(0.05) -> most people in sample are not fit
# Treatment: BP medication (A) ~ Ber with logit(Age, Fitness)
# Cont. Outcome: Systolic BP at follow-up (Y) ~ (Age, Treatment)

# simulate obs data so that for stratum age>60 and fit==0 most treated -> few non-treated
sem2 <- DAG.empty() +
  node("age", distr = "rnorm", mean = 50, sd = 10) +
  node("fit", distr = "rbern", prob = ifelse(age > 60, 0.2, 0.65)) +  # fitness depends on age
  node("A", distr = "rbern", prob = ifelse((age > 60 & fit == 0), 0.95, 1-0.95)) # A1: beta less extreme
  #node("A", distr = "rbern", prob = ifelse((age > 60 & fit == 0), 0.97, 1-0.97)) # A2
  # if old & unfit -> def treated, i.e. P(A=1|old & unfit) should be high
dag2 <- set.DAG(sem2)
plotDAG(dag2)
obs2 <- sim(dag2, rndseed = 30072025, n = 1000)

# check how many observations with age > 60
nrow(obs2[obs2$age > 60 & obs2$fit == 0,])  # as wanted: old & unfit rarely not treated
  # expected viol #1: P(A=0|age >60 & fit == 0) = 0.046, i.e. should be found with b=0.05/0.1, sample size large enough with 130/1000= 0.13

nrow(obs2[obs2$age <= 60 & obs2$A == 1,])/nrow(obs2[obs2$age <= 60,])  # young rarely treated
  # expected viol #2: P(A=1|age <= 75)~0, sample size large enough with 835/1000

nrow(obs2[obs2$fit == 1 & obs2$A == 1, ])/nrow(obs2[obs2$fit == 1, ])  # fit rarely treated
  # expected viol #3: P(A=1|fit == 1) ~0, should be detected for b=0.1, sample size large enough with 556/1000

nrow(obs2[obs2$age <= 60 & obs2$fit == 1 & obs2$A == 1,])/nrow(obs2[obs2$age <= 60 & obs2$fit == 1,])  # young & fit rarely treated
  # expected viol #4 (intersection of the 2 above):
  # only detectable for b = 0.1, sample size sufficient with 521/1000 (NB: not detected if viol 2 or 3 already were)

# categorise age
obs2$age <- cut(obs2$age, breaks = c(0,60,90))
obs2$age <- cut(obs2$age, breaks = c(0,20,40,60,90))


# among A=1
table(obs2[obs2$A == 1, "age"], obs2[obs2$A == 1, "fit"])
# among A=0
table(obs2[obs2$A == 0, "age"], obs2[obs2$A == 0, "fit"])


# PoRT ---
source('data/port_utils.R')
lst2 <- list(port = NULL, port_risca = NULL)
a_values <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.1)
gruber2 <- 5/(sqrt(nrow(obs2))*log(nrow(obs2)))
b_values <- c(0.001, 0.01, gruber2, 0.05, 0.1)
for (a in a_values) {
  for (b in b_values) {
    # uncategorised age
    lst2$port[[paste0("alpha = ", a, ", beta = ", b)]] <- port(A = "A", cov.quanti = "age", cov.quali = "fit",
                                                               data = obs2, alpha = a, beta = b, gamma = 2) 
    # categorised age
    #lst2$port[[paste0("alpha = ", a, ", beta = ", b)]] <- port(A = "A", cov.quanti =NULL, cov.quali = c("age", "fit"),
    #                                                           data = obs2, alpha = a, beta = b, gamma = 2)
    #lst2$port_risca[[paste0("alpha value = ", a)]] <- RISCA::port(group = "A", cov.quanti = "L", cov.quali = NULL, alpha = a, beta = 0.05, data = obs1, gamma = 1)
  }
}
lst2

# for A1:
# without categorisation: identical for gamma = 1 and gamma = 2
  # only split by 1 var, as proba.exposure small enough already -> viol 2+3 identified
  # subgroup fit only from beta = 0.1 as expected
  # viol #4 not detected (would've only been detected for b=0.1 but alr included in viol 2+3)
  # viol #1 covered for a = 0.01, b = 0.1: covered by age -> split by 1 var was sufficient
  # -> but don't understand why for beta = 0.05/0.1 and other alphas (as a=0.13) never reported "age > 60" then, only for a = 0.01 subgroup "age > 60 & age < 65"
# first categorisation: identical for gamma = 1,2 as only ever split by 1 var
  # viol #1 neither detected nor covered
  # too imprecise, identified viol #2 for beta = 0.05, both viol 2+3 for beta = 0.1
# second categorisation:
  # gamma = 1: same as for first categorisation
  # gamma = 2: viol #1 neither detected nor covered
  #            for alpha=0.01-0.05, beta=gruber: new unexpected subgroup! split with 2 vars
  #            intersection fit=0 & age=(20,40], as allows for smaller proba.exposure that only intersection of groups fulfills
  #            else split by 1 var, as above: viol #2 for beta=0.05, viol 2+3 for beta=0.1

# essence: all problematic subgroups detected if not categorised (viol #1 "detected" in terms of covered)
#          if categorised, all found except viol #1: group of age>60 & fit==0 always treated (sample size prop = 0.13) not even covered
#          and 1 new subgroup "fit=0 & age=(20,40]": unfit ppl between 20&40 rarely treated

# for A2:
# without cat: same for gamma = 1,2
#              viol #1 covered (split by 1 var sufficient, that's why no further intersection as in viol #1),
#              but only for alpha = 0.02/0.03 for b = 0.05/0.1, or for alpha = 0.01 for any b
#              but again, why only for small alpha when split into tiny groups, not for alpha = 0.05 as one big group?
#              viol 2,3 also detected
# first cat: viol #1 (old & unfit) neither detected nor covered
#            viol 2+3 for beta = 0.05/0.1
#            new subgroups young & unfit (subset of viol #2) for beta = gruber, and old & fit (subset of viol #3) for alpha <= 0.03
# second cat:
  # gamma = 1: viol #1 neither detected nor covered
  #            viol 2+3 for beta = 0.05/0.1
  # gamma = 2: same as first categorisation

# essence: viol #1 only "detected" in terms of covered, & only if not categorised, no matter how extreme beta (both in A1 & A2)


# KBSD ---
source("kbsd.R")
# revert categorisation
obs2 <- sim(dag2, rndseed = 30072025, n = 1000)

o2 <- obs2[-1]
o2_1 <- o2
o2_1$A <- 1
o2_2 <- o2
o2_2$A <- 0
res2 <- kbsd(data = o2, int_data_list = list(o2_1, o2_2), disthalf_vec=c(age=10, fit=0.5, A=0.5*0.4),  # half of the SD for A!
             plot.out = F)
res2_plot <- kbsd(data = o2,
                  int_data_list = list(o2_1, o2_2),
                  disthalf_vec=c(age=10, fit = 0.5, A=0.2))
res2_plot
# IV=1 encompasses all points that originally had A=0 & artificially got A=1 AND those that originally had A=1
# can tell with table below that prob so few support for IV=1 due to many obs from A=0
#  for which we intervened on AND THAT CANNOT BE WELL SUPPORTED BY THE FEW OBS ALR IN A=1
table(o2$A)  # i.e. imbalance does have impact!!!!!!!
table(cut(obs2$age, breaks = c(0,60,90)))  # also more younger ppl

# inspect what strata are those with low EDP in IV rule 1 (A=1)
subset_3 <- res2[res2$diagnostic < 50 & res2$shift == 1, ]  # most treated obs have low support (771)
sum(obs2[subset_3$observation,][, "age"] > 60)
sum(obs2[subset_3$observation,][, "age"] <= 60) # confirms that subgroup with low support for A=1 has more younger people
table(obs2[subset_3$observation,][, "fit"])  # confirms that those with low support for A=1 are mostly fit ones
obs2[subset_3$observation,] %>% filter(age<=60 & fit == 1) %>% nrow()  # most are young & fit
 # more younger ppl makes sense -> treated far less often (confirms viol #2)
 # fit ppl makes sense -> treated less often (confirms viol #3)

# inspect what strata are those with low EDP in IV rule 2 (A=0)
subset_4 <- res2[res2$diagnostic < 100 & res2$shift == 2, ]
table(obs2[subset_4$observation,][, "fit"], cut(obs2[subset_4$observation,][, "age"], breaks = c(0, 60, 90)))
  # biggest intersecting subgroup is unfit & old -> rarely not treated (confirms viol #1)
# check again how many among unfit & old are treated/untreated
table(obs2[obs2$age > 60 & obs2$fit == 0, "A"])/nrow(obs2[obs2$age > 60 & obs2$fit == 0,])
  # indeed critical as proba.exposure = 0.95 (most are treated) -> should've been detected by port for beta = 0.05/0.1

# essence: better than port, KBSD identified all critical groups (viol 1,2,3)

# compare prop of obs where freq for each treatment level < 0.1:
table(obs2 %>% filter(fit == 1 & age == "(0,60]") %>% select(A))  # pos viol for 0.1
table(obs2 %>% filter(fit == 1 & age == "(60,90]") %>% select(A)) # viol for 0.05
table(obs2 %>% filter(fit == 0 & age == "(0,60]") %>% select(A)) # viol for 0.05
table(obs2 %>% filter(fit == 0 & age == "(60,90]") %>% select(A)) # viol for 0.05
# how to use this as metric, kbsd uses same data so will be same & after constructing
#    hypothetical scenario in kbsd, does not realy make sense

### alternative metrics for disthalf_vec ----




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
table(obs3$L3, obs3$A)/rowSums(table(obs3$L3, obs3$A))  # low P(A=1|L3=0) -> expected viol #1 (sample large enough with 386/1000)
# gamma =2
table(obs3[obs3$L1 == 0 & obs3$L2 == 0,"A"])
table(obs3[obs3$L1 == 0 & obs3$L3 == 0,"A"])  # v imbalanced -> expected viol #2 (sample large enough with 278/1000)
table(obs3[obs3$L2 == 0 & obs3$L3 == 0,"A"])  # v imbalanced -> expected viol #3 (sample large enough with 346/1000)

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
gruber3 <- 5/(sqrt(nrow(obs3))*log(nrow(obs3)))
b_values <- c(0.001, 0.01, gruber3, 0.05, 0.1)
for (a in a_values) {
  for (b in b_values) {
    lst3$port[[paste0("alpha = ", a, ", beta = ", b)]] <- port(A = "A", cov.quanti = NULL, cov.quali = c("L3", "L2", "L1"),
                                                               data = obs3, alpha = a, beta = b, gamma = 3) 
    #lst3$port_risca[[paste0("alpha value = ", a)]] <- RISCA::port(group = "A", cov.quanti = "L", cov.quali = NULL, alpha = a, beta = 0.05, data = obs1, gamma = 1)
  }
}
lst3
# gamma = 1:
  # only for b = 0.1 the one wanted subgroup with L3=0 was detected; makes sense bc proba.exposure=0.073 which is only <0.1, not <0.05 etc.
# gamma = 2:
  # never v small proba.exposure (< gruber)
  # BEST COMBO: any alpha, b = 0.05 to find all 3 viol (beta=gruber/0.1 only find 1 critical subgroup each)
# gamma = 3:
  # same as gamma = 2, but now for smaller beta also found subgroups with L1&L2&L3 intersection
  # reason: proba.exposure not extreme enough for split by 1 var/by 2 vars, so continued to 
  #         build tree with splits by 3 vars where a viol was eventually found

# essence: all violations found, even one more (L3=0 & L1=0 & L2=0)
#          but only if aggregate results from diff HP values, never all 3 viol at once


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
  # detected by kbsd as detected by port!
table(obs3[subset_5$observation, ][, "L1"],
      obs3[subset_5$observation, ][, "L2"],
      obs3[subset_5$observation, ][, "L3"])  # kbsd shows that they most small EDP from there

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




# 4.1) 10 Confounders ----
sem4 <- DAG.empty() +
  node("L1", distr = "rbern", prob = 0.5) +
  node("L2", distr = "rbern", prob = 0.5) +
  # if all =0, then probs=0.5, if all=1, then probs=0.94 (i.e. many that have L1=1, L2=1, L3=1 & few with L1=1, L2=1, L3=0)
  # if L1=0 & L2=1, then probs=0.73, if L1=1 & L2=0, then probs=0.57
  node("L3", distr = "rbern", prob = plogis(0.3*L1 + L2 + L1*L2)) +
  node("L4", distr = "rbern", prob = 0.3) +
  node("L5", distr = "rbern", prob = 0.1) +
  node("L6", distr = "rbern", prob = 0.4) +
  node("L7", distr = "rbern", prob = 0.5) +
  node("L8", distr = "rbern", prob = 0.3) +
  node("L9", distr = "rbern", prob = 0.2) +
  node("L10", distr = "rbern", prob = 0.7) +
  node("A", dist = "rbern", prob = plogis(L1 + 2*L2 + 3*L3 + 0.4*L4 + 0.5*L5*L6))
dag4 <- set.DAG(sem4)
obs4 <- sim(dag4, rndseed = 16092025, n = 1000)
boxplot(obs4[-1])
table(obs4$A)  # v imbalanced alr


# 4.2) 10 Confounders ----
set.seed(22092025)
# L1-L10 (mix of cont & binary)
DAG <- DAG.empty() +
  node("L1", distr = "rnorm", mean = 1, sd = 1) +
  node("L2", distr = "rnorm", mean = 2, sd = 1) +
  node("L3", distr = "rnorm", mean = 3, sd = 1) +
  node("L4", distr = "rnorm", mean = 4, sd = 1) +
  node("L5", distr = "rnorm", mean = 5, sd = 1) +
  node("L6", distr = "rbern", prob = 0.5) +
  node("L7", distr = "rbern", prob = 0.5) +
  node("L8", distr = "rbern", prob = 0.5) +
  node("L9", distr = "rbern", prob = 0.5) +
  node("L10", distr = "rbern", prob = 0.5)
  
#  do not use for loop anymore: returned weird obs -> rnorm obs were all similar 
# for (i in 1:5) {
#   DAG <- DAG + node(paste0("L", i), distr = "rnorm", mean = i, sd = 1)
# }
# for (i in 6:10) {
#   DAG <- DAG + node(paste0("L", i), distr = "rbern", prob = 0.5)
# }

# A dep on a nonlinear function of bernoulli L -> ensures some conf combos lead to low/high treatment probs
DAG <- DAG + node("A", distr = "rbern",
                  prob = plogis(2*L6*L7*L8))  # A only depends on these conf; should help characterising problematic strata
DAG <- set.DAG(DAG)
dat <- sim(DAG, n = 1000)

# table to check for extreme treatment probs only for binary conf: results in table of dimension 2^5
make_strata_table <- function(dat, A = "A", binary_vars){
  dat %>%
    group_by(across(all_of(binary_vars))) %>%
    summarise(
      n      = n(),
      n_treat = sum(.data[[A]] == 1),
      n_control = sum(.data[[A]] == 0),
      proba_exp= mean(.data[[A]] == 1), .groups = "drop") %>%
    arrange(proba_exp)
}
binary_vars <- paste0("L", 6:10) # again, L6...L10 are binary
print(make_strata_table(dat, A = "A", binary_vars = binary_vars), n = 32)
# as wanted: all violating subgroups have L8=1 & L6=1 & L7=1 (NB: overall only have extremely HIGH proba_exp, not extr LOW)
# sample size of L8=1 & L6=1 & L7=1: 0.13, proba_exp = 106/117 = 0.906 i.e. to detect with any alpha & beta = 0.1

### continuous vars uncategorised ----
lst5 <- list(port = NULL)
a_values <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.1)
b_values <- c(0.001, 0.01, 5/(sqrt(nrow(dat))*log(nrow(dat))), 0.05, 0.1)
for (a in a_values) {
  for (b in b_values) {
    lst5$port[[paste0("alpha = ", a, ", beta = ", b)]] <- port(A = "A", cov.quanti = c("L1", "L2", "L3", "L4", "L5"),
                                                               cov.quali = c("L6", "L7", "L8", "L9", "L10"),
                                                               data = dat, alpha = a, beta = b, gamma = 10) 
  }
}
lst5
# gamma = 1,2: predefined critical stratum not yet poss to cover bc intersection of 3
#  -> so other viol among cont conf, but only for alpha=0.01/0.02
#     (where small groups poss -> reason: greedy categorisation esp for small alphas)
#  -> the smaller you alow the subgroups to be, the extremer the pos viol can be defined WITH CONT CONF
# gamma = 3-10: from now on strata by 3 conf, so should have L6=1 & L7=1 & L8=1
#  -> included! BUT never for a=0.02/0.03/0.04 & beta=0.1!


port("A", cov.quali = c("L10", "L9", "L8", "L7", "L6"),
     cov.quanti = c("L3", "L2"), data = dat, alpha = 0.05, beta = 0.1, gamma = 3)
# important finding: order of specifying CONTINUOUS covars as argument matters -> returns diff subgroups!
# e.g. gamma = 3:
#   diff subgroups for c("L1", "L2", "L3", "L4", "L5"), c("L3", "L1", "L2", "L4", "L5"), c("L2", "L3", "L1", "L4", "L5")
#   same for L2, L3 (L8=1 & L6=1 & L7=1 is undetected) and L3, L2 (detected)
# e.g. gamma = 4:
#   undetected if order in continuous vars is L2, L1 -> detected if order is L1, L2
#   same for L2, L3 (undetected) and L3, L2 (detected)
# how is this possible?

### continuous vars categorised ----
# check if problem of specifying order of cov.quanti persists with categorisation
dat_cat <- dat
for (i in names(dat_cat)[-1]) {
  if (typeof(dat_cat[[i]]) == "double") {
    dat_cat[[i]] <- cut(dat_cat[[i]], breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8))
  }
}
lst5_cat <- list(port=NULL)
for (a in a_values) {
  for (b in b_values) {
    lst5_cat$port[[paste0("alpha = ", a, ", beta = ", b)]] <- port(A = "A", cov.quanti = NULL,
                                                               cov.quali = c("L1", "L7", "L8", "L9", "L2", "L3", "L4", "L5",
                                                                             "L6", "L10"),
                                                               data = dat_cat, alpha = a, beta = b, gamma = 3) 
  }
}
lst5_cat
# gamma = 1: no critical subgroup
# gamma = 2: one critical subgroup for alpha = 0.01, but v small
# gamma = 3-8: now detected for all alpha, beta = 0.1 as wanted! shows importance of categorising cont conf!! also computationally faster
# gamma = 9,10  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# problem of changing order in cov.quanti irrelevant here bc now all vars in cov.quali due to categorisation

# KBSD ---
source("kbsd.R")
o5 <- dat[-1]
o5_1 <- o5
o5_1$A <- 1
o5_2 <- o5
o5_2$A <- 0
res5 <- kbsd(data = o5,
             int_data_list = list(o5_1, o5_2),
             disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3 = sd(o5$L3), L4 = sd(o5$L4), L5 = sd(o5$L5),
                            L6 = sd(o5$L6), L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9), L10 = sd(o5$L10),
                            A=0.5*sd(o5$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
res5_plot <- kbsd(data = o5,
                  int_data_list = list(o5_1, o5_2),
                  disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3 = sd(o5$L3), L4 = sd(o5$L4), L5 = sd(o5$L5),
                                 L6 = sd(o5$L6), L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9), L10 = sd(o5$L10),
                                 A=0.5*sd(o5$A)))
res5_plot
# overall v few EDP as high-dim adjustment set

# what strata are those with low support for treated (IV = 1)
subset_5_1 <- res5[res5$diagnostic < median(res5[res5$shift == 1, "diagnostic"]) & res5$shift == 1,]
o5[subset_5_1$observation, ]
# did not built in viol here on purpose so would have to search in detail which underlying covar values

# what strata are those with low support for treated (IV = 0): should be L8=1 & L6=1 & L7=1
subset_5_2 <- res5[res5$diagnostic < median(res5[res5$shift == 2, "diagnostic"]) & res5$shift == 2,]
o5[subset_5_2$observation, ]
table(o5[subset_5_2$observation, ][, c("L6", "L7", "L8")]) # highest count among untreated is L8=1 & L6=1 & L7=1 -> as expected


# disthalf_vec: best case values, i.e. double values <=> kernel less steep
res5 <- kbsd(data = o5,
             int_data_list = list(o5_1, o5_2),
             disthalf_vec=c(L1=2*sd(o5$L1), L2 = 2*sd(o5$L2), L3 = 2*sd(o5$L3), L4 = 2*sd(o5$L4), L5 = 2*sd(o5$L5),
                            L6 = 2*sd(o5$L6), L7 = 2*sd(o5$L7), L8 = 2*sd(o5$L8), L9 = 2*sd(o5$L9), L10 = 2*sd(o5$L10),
                            A=2*0.5*sd(o5$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
res5_plot <- kbsd(data = o5,
                  int_data_list = list(o5_1, o5_2),
                  disthalf_vec=c(L1=2*sd(o5$L1), L2 = 2*sd(o5$L2), L3 = 2*sd(o5$L3), L4 = 2*sd(o5$L4), L5 = 2*sd(o5$L5),
                                 L6 = 2*sd(o5$L6), L7 = 2*sd(o5$L7), L8 = 2*sd(o5$L8), L9 = 2*sd(o5$L9), L10 = 2*sd(o5$L10),
                                 A=2*0.5*sd(o5$A)))
res5_plot
# table for tracing back critical strata (IV=0) more balanced now -> L8=1 & L6=1 & L7=1 not particularly problematic

# " " : worst case values, i.e. halve values <=> kernel steeper
res5 <- kbsd(data = o5,
             int_data_list = list(o5_1, o5_2),
             disthalf_vec=c(L1=0.5*sd(o5$L1), L2 = 0.5*sd(o5$L2), L3 = 0.5*sd(o5$L3), L4 = 0.5*sd(o5$L4), L5 = 0.5*sd(o5$L5),
                            L6 = 0.5*sd(o5$L6), L7 = 0.5*sd(o5$L7), L8 = 0.5*sd(o5$L8), L9 = 0.5*sd(o5$L9), L10 = 0.5*sd(o5$L10),
                            A=0.5*0.5*sd(o5$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
res5_plot <- kbsd(data = o5,
                  int_data_list = list(o5_1, o5_2),
                  disthalf_vec=c(L1=0.5*sd(o5$L1), L2 = 0.5*sd(o5$L2), L3 = 0.5*sd(o5$L3), L4 = 0.5*sd(o5$L4), L5 = 0.5*sd(o5$L5),
                                 L6 = 0.5*sd(o5$L6), L7 = 0.5*sd(o5$L7), L8 = 0.5*sd(o5$L8), L9 = 0.5*sd(o5$L9), L10 = 0.5*sd(o5$L10),
                                 A=0.5*0.5*sd(o5$A)))
res5_plot
# table for tracing back critical strata (IV=0) most imbalanced now -> prominence of L8=1 & L6=1 & L7=1 v striking

# over all disthalf_vec specifications: same general trend with better support for IV=1, but EDP values quite diff

### alternative metrics for disthalf vec ----




# no pattern: as soon as >2, stratum is not detected, but it is again for c("L1", "L3", "L2", "L4")

# seems like detected -> need to take closer look and compare
make_strata_table(dat, A = "A", binary_vars = binary_vars) %>%
  filter((L7==1 & L8==0 & L6==1) |
         (L7==1 & L6==1 & L9==1) |
         (L7==1 & L8==0 & L9==1) | L10 == 1) %>%
  filter(proba_exp <= 0.1 | proba_exp >= 0.9)  # indeed all critical treatment probs; most covered by L10==1


sem5 <- DAG.empty() +
  node("L1", distr = "rbern", prob = 0.3) +
  node("L2", distr = "rbern", prob = 0.1) + 
  node("L3", distr = "rbern", prob = 0.6) +
  node("L4", distr = "rbern", prob = 0.4) +
  node("L5", distr = "rbern", prob = 0.5) +
  node("A", distr = "rbern", prob = 0.1*L1 + 0.1*L2 + 0.1*L3 + 0.1*L4 + 0.6*L5)  # if L1, L2, L3=1 -> prob A=1
dag5 <- set.DAG(sem5)
plotDAG(dag5)
obs5 <- sim(dag5, rndseed = 12082025, n = 1000)
port("A", cov.quali = c("L1", "L2", "L3", "L4", "L5"), cov.quanti = NULL,
                        data = obs5, alpha = 0.05, beta = 0.1, gamma = 5)


# 5) 20 Confounders ----
DAG2 <- DAG.empty()
# L1-L10 (mix of continuous & binary for realism)
for (i in 1:10) {
  DAG2 <- DAG2 + node(paste0("L", i), distr = "rnorm", mean = i, sd = 1)
}
for (i in 11:20) {
  DAG2 <- DAG2 + node(paste0("L", i), distr = "rbern", prob = 0.5)
}
DAG2 <- set.DAG(DAG2)
dat2 <- sim(DAG2, rndseed = 12082025, n = 1000)[-1]
port("A", cov.quanti = c("L1","L2","L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10"),
     cov.quali = c("L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19", "L20"),  # why is this a problem??
     data = dat2, alpha = 0.05, beta = 0.1, gamma = 100)


# 6) 50 Confounders ----
DAG2 <- DAG.empty()
# L1-L10 (mix of continuous & binary for realism)
for (i in 1:25) {
  DAG2 <- DAG2 + node(paste0("L", i), distr = "rnorm", mean = i, sd = 1)
}
for (i in 26:50) {
  DAG2 <- DAG2 + node(paste0("L", i), distr = "rbern", prob = 0.5)
}
DAG2 <- set.DAG(DAG2)
dat2 <- sim(DAG2, rndseed = 12082025, n = 1000)[-1]
port("A", cov.quanti = c("L1","L2","L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10",
                         "L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19", "L20",
                         "L21", "L22", "L23", "L24" ,"L25"),
     cov.quali = c("L26", "L27", "L28", "L29", "L30", "L31", "L32", "L33", "L34", "L35", "L36",
                   "L37", "L38", "L39", "L40", "L41", "L42", "L43", "L44", "L45", "L46", "L47",
                   "L48", "L49", "L50"),  # why is this a problem??
     data = dat2, alpha = 0.05, beta = 0.1, gamma = 100)
