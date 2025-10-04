
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

# check manually for cont conf
obs2 %>% filter(L1 > 60 & A==1) %>% nrow()/obs2 %>% filter(L1 > 60) %>% nrow()  # viol #1: P(A)~1 for b=0.1, any a
obs2 %>% filter(L1 < 60 & A==1) %>% nrow()/obs2 %>% filter(L1 < 60) %>% nrow()

obs2 %>% filter(L2 == 0 & A==1) %>% nrow()/obs2 %>% filter(L2 ==0) %>% nrow()  # viol #2: P(A)~1 for beta=0.05/0.1, any a
obs2 %>% filter(L1 > 60 & L2 == 0 & A==1) %>% nrow()/obs2 %>% filter(L1 > 60 & L2 ==0) %>% nrow()  # viol #3: combo of the 2 above with beta = 0.01, any a!

obs2 %>% filter(L2 == 1 & A==1) %>% nrow()/obs2 %>% filter(L2 ==1) %>% nrow()  # viol #4: P(A)~0 for beta=0.1, any a

# check via table for all combos of binary conf if there are any with extreme P(A): function defined in setup.R
tab <- make_strata_table(obs2, A = "A", binary_vars = "L2")
tab %>% filter((proba_exp <= 0.1 | proba_exp >= 0.9) & sample_prop >= 0.01)
# confirms viol 2+4 & sample prop large enough

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
# always det viol #1 for any a & b, bc L1 as only cont var that can be split into smaller subgroups
# always det viol #2 for any a & b = 0.05/0.1 as should; always det viol #4 for any a & b =0.1 as should
# viol #3 (intersecting group) poss to detect for gamma =2, but only covered bc of greedy cat for L1
# i.e. best choice to detected all: a=0.05/0.1 & b=0.1 (bc for smaller a other "irrelevant" strata returned due to greedy cat of L1)

# categorised:
obs2$L1 <- cut(obs2$L1, breaks = c(0,20,40,60,90))
# gamma = 1: same as above, always detected from b = 0.1 (viol #2 from 0.05) as should
# gamma = 2: now also intersecting group detected! bc no greedy cat by L1, rest always det for b=0.1 as wanted
# g=1,2: only viol 1,2,4 all detected for any a, b=0.1! do not have "irrelevant" small strata for small a anymore
# bc L1 categorised now

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
# visibly more support for IV=2 (A=0), overall range between 0-400
table(obs2$A)

# expected viol (confirmed by PoRT): P(A=1|age > 60)~1, P(A=1|fit = 0)~1, P(A=1|fit = 1)~0

# what strata are those with low EDP in IV = 1 (A=1) -> expect fit=1
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

# what strata are those with low EDP in IV = 2 (A=0) -> expect age>60 & fit=0
shift2 <- res2[res2$shift == 2,]
outliers2 <- shift2$diagnostic < quantile(res2[res2$shift==2, "diagnostic"], probs = 0.25)
l_values2 <- obs2[outliers2, "L1"] 
diag_values2 <- shift2[outliers2,]
plot(l_values2, diag_values2$diagnostic)
# most with few EDP are v old pp (but also some v young ones!) -> sensible bc few old ppl
# that got A=0 so that few support for such obs if IV=2 (A=0)
l_values2 <- obs2[outliers2, "L2"]
diag_values2 <- shift2[outliers2,]
plot(l_values2, diag_values2$diagnostic)
# most with few support for A=0 are those with L2=0, i.e. unfit ppl bc unfit ppl are
# usually treated so that if intervene on A=0 for them, few unfit obs in neighbourhood

# essence: det viol 1,2,4 & even hints at obs with age <50 having few support for A=1
# (although not extreme enough to be considered a viol by PoRT, still important insight)


# 3 Confounders ----
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

# table to check for all combos of binary conf if there are any with extreme P(A): def in setup.R
binary_vars <- paste0("L", 1:3)
tab <- make_strata_table(obs3, A = "A", binary_vars = binary_vars)
tab %>% filter((proba_exp <= 0.1 | proba_exp >= 0.9) & sample_prop >= 0.01)
# g=1: low P(A=1|L3=0) -> viol #1 for b=0.1 (sample large enough with 38.6%)
# g=2: L1=0 & L3=0 -> viol #2: P(A=1)~0 for b=gruber (sample large enough with 27.8%)
#      L2=0 & L3=0 -> viol #3: P(A=1)~0 for b=0.1 (34.6%) -> will prob be covered by viol #1
# g=3: L1=1 & L2=1 & L3=1 -> viol #4: P(A=1)~1 for a=0.01 & b=0.1
#      L1=0 & L2=0 & L3=0 -> viol #5: P(A=1)~0 for any a & b


# PoRT ---
source('data/port_utils.R')
a_values <- c(0.01, 0.025, 0.05, 0.1)
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
#sink("port_3_binary.txt")
# gamma = 1:
# only for b = 0.1 the one wanted subgroup with L3=0 was detected; makes sense bc proba.exposure=0.073 which is only <0.1, not <0.05 etc.
# gamma = 2:
# all 3 viol detected for appropriate beta values
# BEST COMBO: any alpha, b = 0.05 to find all 3 viol (beta=gruber/0.1 only find 1 critical subgroup each)
# gamma = 3:
# for a=0.01, b=0.1, viol #4 with L1==1 & L2==1 & L3==1 not found

# essence: viol 1,2,5 always det, #3 always covered by #1, only #4 is missing


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
# overall EDP range between 0-300

# what strata are those with low support for treated (IV = 1): should be viol 1,2,3,5
subset_5 <- res3[res3$diagnostic < quantile(res3[res3$shift==1, "diagnostic"], probs = 0.5)
                 & res3$shift == 1,]
obs3[subset_5$observation, ]
table(obs3[subset_5$observation, ][, c("L1", "L2")]) # many L1=0 & L2=0 among treated
obs3 %>% filter(L1==0 & L2==0 & A==1) %>% nrow()/obs3 %>% filter(L1==0 & L2==0) %>% nrow() 
   # while few support (prob tied to viol #5), P(A=1|L1=0 & L2=0) not extreme
table(obs3[subset_5$observation, ][, c("L1", "L3")]) # confirmed that many from L3=0 & L1=0 (#2), ACTUALLY from L3=0 (#1)
table(obs3[subset_5$observation, ][, c("L2", "L3")]) # confirmed that many from L3=0 & L2=0 (#3)

# check for viol #5:
table(obs3[subset_5$observation, ][, "L1"],
      obs3[subset_5$observation, ][, "L2"],
      obs3[subset_5$observation, ][, "L3"])
# kbsd shows that most with small EDP from L3=0 & L1=0 & L2=0 -> det viol #5


# what strata are those with low support for untreated (IV = 2): should be viol #4
subset_6 <- res3[res3$diagnostic < quantile(res3[res3$shift==2, "diagnostic"], probs = 0.25)
                 & res3$shift == 2,]
obs3[subset_6$observation, ] %>% filter(L1==1 & L2==1 & L3==1)
table(obs3[subset_6$observation, ][, "L1"],
      obs3[subset_6$observation, ][, "L2"],
      obs3[subset_6$observation, ][, "L3"]) # split by L3: only 18 obs of viol #4 (not well det)
# instead flagged L3=1 & L1=0 & L2=1 as most obs with low EDP and L3=0 & L1=0 & L2=1
# check among L3=1 & L1=0 & L2=1:
table(obs3[obs3$L1 == 0 & obs3$L2 == 1 & obs3$L3 == 1, "A"])  # not extreme enough: proba.exposure = 0.78
# check among L3=0 & L1=0 & L2=1:
table(obs3[obs3$L1 == 0 & obs3$L2 == 1 & obs3$L3 == 0, "A"])  # not extreme enough: proba.exposure = 0.16
table(o3$A)  # also: way more untreated compared to treated in general!

# essence: KBSD found all 4 viol among A=1
#          but for A=0, just like PoRT failed to det viol #4 (L1==1 & L2==1 & L3==1),
#          instead flagged 2 new strata among A=0 that do not present problem actually



# 3 Confounders With Middle-Gap ----

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
  data1 %>% filter(L3 >= 4 & L3 <= 6) %>% nrow()  # viol #1: P(A=1)~0 if L3 in [4,6], sample prop = 14.1% -> should find for g>=1, b=0.1

data1 %>% filter(L3 >= 4 & L3 <= 6 & L1==0 & A==1) %>% nrow()/
  data1 %>% filter(L3 >= 4 & L3 <= 6 & L1 == 0) %>% nrow()  # viol #2 (subviol of #1): if L3 in [4,6] AND L1=0, sample prop = 12.9% -> g>=2, b=0.1

data1 %>% filter(L3 >= 4 & L3 <= 6 & L2==0 & A==1) %>% nrow()/
  data1 %>% filter(L3 >= 4 & L3 <= 6 & L2 == 0) %>% nrow()  # viol #3 (subviol of #1): if L3 in [4,6] AND L2=0, sample prop = 11.3% -> g>=2, b=0.1


# table to check for all combos of binary conf if there are any with extreme P(A): def in setup.R
binary_vars <- paste0("L", c(1,2))
tab <- make_strata_table(data1, A = "A", binary_vars = binary_vars)
tab
# no extreme proba only with L1, L2


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
# gamma =2: viol #3 with combo of L2=0 and L3< 5.626 & L3>=5.199 only for a=0.025 & b=0.01/gruber
#           bc smaller stratum for L3 so that proba.exp extremer there
#           viol #2 always covered by viol #1 & viol #3 also covered by #1 for larger a=0.05/0.1 & b=0.1
# gamma =3: always identified where appropriate
# but wrong numbers for several subgroups? ex. g=1, a=0.01, b=0.1: ------
# for viol 1,2,5,6,7 the sample prop are wrong & for viol 2,5,6 the proba are wrong


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


# essence: overall viol #1 found in both cat & uncat, but when uncat,
#          can return preciser subgroup that includes L2=0, whereas after cat
#          the broader subgroup only is returned -> so emphasises importance of cat thresholds
#          appropriate/meaningful for the context as Chatton et al also emphasised -> to limit info loss
#          however, problem with wrong proba.exposure computed by PoRT in uncat case.. 


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

# expect subgroup L3=(4,6] to have few support for IV=1 (A=1) if would estimate Y|A=1 further:
shift1 <- res5[res5$shift == 1,]
outliers1 <- shift1$diagnostic < quantile(shift1$diagnostic, probs = 0.25)  # create indices for the "outliers"
l_values1 <- data1[outliers1, "L3"]  # to which original obs (L values) do these outliers belong?
diag_values1 <- shift1[outliers1,] # what diag values do these outliers have
plot(l_values1, diag_values1$diagnostic)

# interesting: for threshold =0.05, those in [4,6] do not have lowest EDP/support, 
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
# no clear pattern except less support for tails bc fewer obs there (is ok, did not have viol with P(A=0)~0)


# essence: kbsd cannot PINPOINT that for IV=1 there is low support for L3 =[4,6], i.e. the "hole" whereas PoRT could
