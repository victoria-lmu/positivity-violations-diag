
# Simulating causal settings with 1,2 and 3 confounders with sparse data for some strata
# and checking if different diagnostics detect these strata with positivity violations.
source("setup.R")
set.seed(15082025)

## 1) One Confounder ----
# Baseline covariate: Health score (L) ~ Normal(0, 1) -> neg = bad, pos = healthy
# Treatment: Vaccination (A) ~ Ber with logit depending on L
# Outcome: Hospitalisation (Y) ~ Ber depending on A and L

# simulate so that neg health score <=> high probs of treatment, pos health score <=> low probs of treatment
# expected viol #1: P(A=1|neg health)~1 -> among all with neg L, i.e. if seen as one category as lateron: P(A=1) = 0.97, subgroup size = 51.9% -> from beta = gruber
# expected viol #2: P(A=1|pos health)~0 -> among all with pos L i seen as one cat: P(A=0) = 0.93, subgroup size = 48.1 % -> to be detected for beta = 0.1

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
b_values <- c(0.01, gruber, 0.05, 0.1)
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
 # both always covered, though for small alpha by dividing into small groups, for big alpha in broader groups
 # BEST COMBO: a=0.05, b=0.05 or a=0.05, b=0.1 (= a=0.1,b=0.1) to be stricter with probs for pos viol
 # a = 0.1, b = 0.05 already misses 0 < L < 1, because P(A|0<L<1) = 0.07 -> only detected with b=0.1
# with binary categorisation into neg & pos health:
 # both detected from b=0.1 as they should
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
# stratum with few support if intervened on A=1: L>0 bc few obs there with L>0 & A=1
#    "     "        "             "  on A=0: L<0 bc few obs there with L<0 & A=0 that could give intervened-on obs some support
# but ex. IV = 2: when checking for EDP threshold =100, not all neg are identified
# i.e. the higher the EDP threshold, the more matches with port results, the more is "covered"
# but kbsd, esp plotting L values against diag values v useful for first identification of critical strata!


# best case disthalf vec: doubled
res1_plot <- kbsd(data = o1,
                  int_data_list = list(o1_1, o1_2),
                  disthalf_vec= c(L=2*1, A=2*0.5*0.5))
res1_plot

# worst case disthalf vec: halved
res1_plot <- kbsd(data = o1,
                  int_data_list = list(o1_1, o1_2),
                  disthalf_vec= c(L=0.5*1, A=0.5*0.5*0.5))
res1_plot

# essence: always same trend across IV levels (generally better support for IV=1)
# plots for strata also yield same results for underlying L values


### alternative metrics for disthalf vec ----
# 1) MAD within L, A instead of SD
res1_plot <- kbsd(data = o1,
                  int_data_list = list(o1_1, o1_2),
                  disthalf_vec=c(L=mad(o1$L), A=0.5*mad(o1$A)))
res1_plot
# unsuitable: assume that the error is due to mad(A) = 0 that leads to kernel calc
# crushing down (bc mad(A)² in denominator)? but there is protection mechanism actually

# 2) IQR for L, A instead of SD
res1_plot <- kbsd(data = o1,
                  int_data_list = list(o1_1, o1_2),
                  disthalf_vec=c(L=IQR(o1$L), A=0.5*IQR(o1$A)))
res1_plot
# slightly better support overall (more optimistic) than for regular calc

# 3) average pairwise distance within L, A
res1_plot <- kbsd(data = o1,
                  int_data_list = list(o1_1, o1_2),
                  disthalf_vec=c(L=mean(dist(matrix(o1$L), method = "euclidean")),
                                 A=0.5*mean(dist(matrix(o1$A), method = "euclidean"))))
res1_plot
# similar to regular plot

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
  #node("A", distr = "rbern", prob = ifelse((age > 60 & fit == 0), 0.95, 1-0.95)) # A1: beta less extreme
  node("A", distr = "rbern", prob = ifelse((age > 60 & fit == 0), 0.97, 1-0.97)) # A2
  # if old & unfit -> def treated, i.e. P(A=1|old & unfit) should be high
dag2 <- set.DAG(sem2)
plotDAG(dag2)
obs2 <- sim(dag2, rndseed = 30072025, n = 1000)
table(obs2$A)  # v imbalanced! 

# check how many observations with age > 60
nrow(obs2[obs2$age > 60 & obs2$fit == 0,])  # as wanted: old & unfit rarely not treated
  # expected viol #1: P(A=0|age >60 & fit == 0) = 0.046 (0.023), i.e. should be found with b=0.05/0.1, sample size large enough with 130/1000= 0.13

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
a_values <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.1)
gruber2 <- 5/(sqrt(nrow(obs2))*log(nrow(obs2)))
b_values <- c(0.01, gruber2, 0.05, 0.1)
g_values <- 1:2
lst2 <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      # uncategorised age
      #lst2[[paste0("gamma = ",g, ", alpha = ", a, ", beta = ", b)]] <-
      #  port(A = "A", cov.quanti = "age", cov.quali = "fit", data = obs2, alpha = a, beta = b, gamma = g) 
      # categorised age
      lst2[[paste0("gamma = ",g, ", alpha = ", a, ", beta = ", b)]] <-
      port(A = "A", cov.quanti =NULL, cov.quali = c("age", "fit"), data = obs2, alpha = a, beta = b, gamma = g)
    }
  }
}
lst2

# for A1:
# without categorisation: identical for gamma = 1 and gamma = 2
  # only split by 1 var, as proba.exposure small enough already -> viol 2+3 identified
  # subgroup fit only from beta = 0.1 as expected
  # viol #4 covered by viol 2+3
  # viol #1 partly covered for a = 0.01, b = 0.1, by age -> split by 1 var was sufficient
  # -> but for other alphas (as a=0.13) with beta = 0.05/0.1, never reported "age > 60", only for a = 0.01 subgroup "age > 60 & age < 65"
  # i.e. viol #1 undetected for a=0.02-0.1
# first categorisation: identical for gamma = 1,2 as only ever split by 1 var
  # viol #1 neither detected nor covered for any alpha & beta =0.1, viol 2+3 for beta = 0.1
# second categorisation:
  # gamma = 1: same as for first categorisation
  # gamma = 2: viol #1 neither detected nor covered
  #            for alpha=0.01-0.05, beta=gruber: new unexpected subgroup! split with 2 vars
  #            fit=0 & age=(20,40] rarely treated -> confirmed as P(A=1) = 0.016
  #            as above: det viol #2 for beta=0.05, viol 2+3 for beta=0.1 as should be

# essence: all viol except viol #1 "detected"/covered both for uncat & cat (i.e. effect of cat
#          not that visible yet with few conf) -> viol #1 only partly covered if uncategorised
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
table(o2$A)
# this table shows that prob so few support for IV=1 due to many obs from A=0
#  for which we intervened on AND THAT CANNOT BE WELL SUPPORTED BY THE FEW OBS ALR IN A=1
#  -> i.e. low EDP mainly caused by difference in A for many obs with A=0 getting IV=1, additionally stems from age/fit
#  -> if obs in one treatment level v scarce alr (A=1), then intervening obs from other on it results in few support for them (in IV=1)
# check if few support mainly stems from age or fit
table(cut(obs2$age, breaks = c(0,60,90)))  # distr of A goes hand in hand with
#  distr of age: those many ppl having A=0 are all mostly young ppl

# inspect what strata are those with low EDP in IV = 1 (A=1)
# take c=50 as EDP threshold here, but depends on scenario (maybe should take same threshold for both BPs?)
subset_3 <- res2[res2$diagnostic < 50 & res2$shift == 1, ]
sum(obs2[subset_3$observation,][, "age"] > 60)
sum(obs2[subset_3$observation,][, "age"] <= 60) # confirms that stratum with low support for IV=1 has more younger people
table(obs2[subset_3$observation,][, "fit"])  # confirms that those with low support for A=1 are mostly fit ones
obs2[subset_3$observation,] %>% filter(age<=60 & fit == 1) %>% nrow()  # most are young & fit: 521/776
 # more younger ppl makes sense -> treated less often so few supporting obs (confirms viol #2)
 # fit ppl makes sense -> also treated less often (confirms viol #3)

# inspect what strata are those with low EDP in IV rule 2 (A=0)
subset_4 <- res2[res2$diagnostic < 100 & res2$shift == 2, ]
table(obs2[subset_4$observation,][, "fit"],
      cut(obs2[subset_4$observation,][, "age"], breaks = c(0, 60, 90)))
# most with low support in IV=0 are unfit + old ppl -> few such obs have received A=0
#  (makes sense bc should be treated), so few support (confirms viol #1)
# check again how many among unfit & old are treated/untreated
table(obs2[obs2$age > 60 & obs2$fit == 0, "A"])/nrow(obs2[obs2$age > 60 & obs2$fit == 0,])
  # indeed most are treated, even critical as proba.exposure = 0.95 -> should've been detected by port for beta = 0.05/0.1

# essence: better than port, KBSD identified all critical groups (viol 1,2,3)

# compare prop of obs where freq for each treatment level < 0.1: requires categorisation
table(obs2 %>% filter(fit == 1 & age == "(0,60]") %>% select(A))  # pos viol for 0.1
table(obs2 %>% filter(fit == 1 & age == "(60,90]") %>% select(A)) # viol for 0.05
table(obs2 %>% filter(fit == 0 & age == "(0,60]") %>% select(A)) # viol for 0.05
table(obs2 %>% filter(fit == 0 & age == "(60,90]") %>% select(A)) # viol for 0.05
# how to use this as metric? kbsd uses same data so will be same & after constructing
#    hypothetical scenario in kbsd, does not really make sense

### alternative metrics for disthalf_vec ----

# best case disthalf vec: doubled
res2_plot <- kbsd(data = o2,
                  int_data_list = list(o2_1, o2_2),
                  disthalf_vec=c(age=2*10, fit = 2*0.5, A=2*0.5*sd(o2$A)))
res2_plot # makes sense that overall higher EDP bc kernel info higher on avg for every point now

# worst case disthalf vec: halved
res2_plot <- kbsd(data = o2,
                  int_data_list = list(o2_1, o2_2),
                  disthalf_vec=c(age=0.5*10, fit = 0.5*0.5, A=0.5*0.5*sd(o2$A)))
res2_plot

# 1) MAD within L, A instead of SD
# again unsuitable: mad(A) = 0

# 2) IQR for L, A instead of SD
res2_plot <- kbsd(data = o2,
                  int_data_list = list(o2_1, o2_2),
                  disthalf_vec=c(age=IQR(o2$age), fit = IQR(o2$fit), A=0.5*IQR(o2$A)))
res2_plot
# not possible bc due to imbalanced treatment distr, most have A=0, i.e. IQR(A)=0

# 3) average pairwise distance within L, A
res2_plot <- kbsd(data = o2,
                  int_data_list = list(o2_1, o2_2),
                  disthalf_vec=c(age=mean(dist(matrix(o2$age), method = "euclidean")),
                                 fit=mean(dist(matrix(o2$fit), method = "euclidean")),
                                 A=0.5*mean(dist(matrix(o2$A), method = "euclidean"))))
# NB: method=euclid <=> manhattan here, bc only 1D -> sqrt((age[1]-age[2])²)=age[1]-age[2]
res2_plot



## 2) Two Confounders But Balanced ----
# Baseline confounders: Age (L1) ~ N(50, 10) and Fitness (L2) ~ Ber(0.05) -> most people in sample are not fit
# Treatment: BP medication (A) ~ Ber with logit(Age, Fitness)
# Cont. Outcome: Systolic BP at follow-up (Y) ~ (Age, Treatment)
set.seed(15082025)
sem2 <- DAG.empty() +
  node("age", distr = "rnorm", mean = 50, sd = 10) +
  node("fit", distr = "rbern", prob = ifelse(age > 60, 0.2, 0.7)) +  # fitness depends on age
  node("A", distr = "rbern", prob = plogis(4*(age > 60) + 2.5*(1-fit) - 3*fit))
   # if unfit -> P(A=1)=0.92, if old -> P(A=1) = 0.98, if unfit + old -> P(A=1)~1, if fit -> P(A=1)=0.047
# if old & unfit -> def treated, i.e. P(A=1|old & unfit) should be high
dag2 <- set.DAG(sem2)
plotDAG(dag2)
obs2 <- sim(dag2, rndseed = 30072025, n = 1000)
table(obs2$A)  # more balance now

obs2 %>% filter(age > 60 & A==1) %>% nrow()/obs2 %>% filter(age > 60) %>% nrow()  # viol #1: for beta=0.1

obs2 %>% filter(fit == 0 & A==1) %>% nrow()/obs2 %>% filter(fit ==0) %>% nrow()  # viol #2: for beta=0.05,0.1

obs2 %>% filter(age > 60 & fit == 0 & A==1) %>% nrow()/obs2 %>% filter(age > 60 & fit ==0) %>% nrow()  # viol #3: combo of the 2 above with beta = 0.01!

obs2 %>% filter(fit == 1 & A==1) %>% nrow()/obs2 %>% filter(fit ==1) %>% nrow()  # viol #3: for beta=0.1



# PoRT ---
source('data/port_utils.R')
a_values <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.1)
gruber2 <- 5/(sqrt(nrow(obs2))*log(nrow(obs2)))
b_values <- c(0.01, gruber2, 0.05, 0.1)
g_values <- 1:2
lst2 <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      # uncategorised age
      #lst2[[paste0("gamma = ",g, ", alpha = ", a, ", beta = ", b)]] <-
      #  port(A = "A", cov.quanti = "age", cov.quali = "fit", data = obs2, alpha = a, beta = b, gamma = g) 
      # categorised age
      lst2[[paste0("gamma = ",g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A", cov.quanti =NULL, cov.quali = c("age", "fit"), data = obs2, alpha = a, beta = b, gamma = g)
    }
  }
}
lst2
# uncategorised:
# gamma = 1, 2:
#   a=0.01/0.02 & b<=gruber: only viol #1, bc age as only cont var that can be split into smaller subgroups
#   from a=0.03, all detected (or say viol #3 covered by detection of 1,2) for b=0.1 (viol #2 from 0.05) as they should be!
#   viol #3 (intersecting group) poss to detect for gamma =2, but only covered bc of greedy cat for age

# categorised:
obs2$age <- cut(obs2$age, breaks = c(0,20,40,60,90))
# gamma = 1: same as above, always detected from beta = 0.1 (viol #2 from 0.05) as should
# gamma = 2: now also intersecting group detected! bc no greedy cat by age, rest always det for b=0.1 as wanted


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




## 3) Three Confounders ----
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

# (OLD!! bc only intersections with all 3, not with 2 as well) table to check for binary conf where P(A) extreme
make_strata_table <- function(dat, A = "A", binary_vars){
  dat %>%
    group_by(across(all_of(binary_vars))) %>%
    summarise(
      n      = n(),
      n_treat = sum(.data[[A]] == 1),
      n_control = sum(.data[[A]] == 0),
      proba_exp= mean(.data[[A]] == 1), .groups = "drop",
      sample_prop = n()/nrow(dat)) %>%
    arrange(proba_exp)
}
binary_vars <- paste0("L", 1:3)
print(make_strata_table(obs3, A = "A", binary_vars = binary_vars), n = 8)  # all =0 or all =1, both sample prop large enough 



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

