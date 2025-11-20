source("setup.R")

# 5 Confounders ----

set.seed(28092025)
sem1 <- DAG.empty() +
  node("L1", distr = "rbern", prob = 0.1) +
  node("L2", distr = "rbern", prob = 0.2) +
  node("L3", distr = "rnorm", mean = 6, sd = 1) +
  node("L4", distr = "rbern", prob = 0.1) +
  node("L5", distr = "rbern", prob = 0.2) +
  node("A", distr = "rbern", prob = plogis(0.2*L1 + 0.3*L2 + 0.1*L4 + 0.3*L5 - 2*L3*(L3 < 6)))
# all conf L_i=1 is v unlikely -> treatment not too likely, too but not extreme either
dag1 <- set.DAG(sem1)
plotDAG(dag1)
data1 <- sim(dag1, n = 1000, rndseed = 28092025)

# check for extreme probs with cont conf L3:
data1 %>% filter(L3<6 & A==1) %>% nrow()/data1 %>% filter(L3<6) %>% nrow()  # viol: for g >=1, any a & any b as sample prop = 52%

# check via table for extreme probs for all combos of binary conf: def in setup.R
binary_vars <- paste0("L", c(1,2,4,5))
tab <- make_strata_table(data1, A = "A", binary_vars = binary_vars)
tab %>% filter((proba_exp <= 0.1 | proba_exp >= 0.9) & sample_prop >= 0.01)
# 2 more viol subgroups for a = 0.01 expected

## PoRT: continuous var L3 uncategorised ----
source("data/port_utils.R")
lst5 <- list()
a_values <- c(0.01, 0.025, 0.05, 0.1)
b_values <- c(0.01, 5/(sqrt(nrow(data1))*log(nrow(data1))), 0.05, 0.1)
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
#sink("output_port/port_5_uncat.txt")
# gamma = 1-5: always det when should


## PoRT: continuous var L3 categorised ----
source("data/port_utils.R")
data1_cat <- data1
data1_cat$L3 <- cut(data1_cat$L3, breaks = c(-Inf, 1, 2, 3, 4, 5, 6, 7, 8, Inf))
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
#sink("output_port/port_5_cat.txt")
# g = 1-5: also always det

## kbsd ----
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
res5_plot  # EDP range between 0-200, fewer support for IV=1
table(data1$A)  # which alr indicated here by fewer obs in A=1

# identify obs with low EDP: subgroup L3<6 with P(A=1)~0 should have few support for IV=1 (A=1) if would estimate Y|A=1 further
indices_outliers <- unique(res5[res5$shift == 1 & res5$diagnostic < quantile(res5$diagnostic, 0.05), "observation"])
o5[indices_outliers, ] %>% nrow()
o5[indices_outliers, ] %>% filter(L3 < 6) %>% nrow()
# q=0.05: 77/93 from critical stratum -> majority of obs reported as critical are from viol stratum!
# q=0.01: 19/20 from critical stratum!
o5[indices_outliers, ] %>% filter(!(L3 < 6))
# q=0.05: quite a lot from (L1== 0 & L2==1 & L4==1 & L5==1) -> indeed flagged initially by make_strata_table

# essence: both PoRT & kbsd det viol


## KBSD cat ----
table(data1_cat$L3)
data1_cat <- data1_cat %>% 
  mutate(L3 = case_when(L3 == "(-Inf,1]" ~ 1, L3 == "(1,2]" ~ 2, L3 == "(2,3]" ~3, 
                        L3 == "(3,4]" ~4, L3 =="(4,5]" ~ 5, L3 == "(5,6]" ~6,
                        L3 == "(6,7]" ~7, L3 == "(7,8]"~8, L3 == "(8, Inf]"~9))
make_strata_table(data1_cat, binary_vars = c(paste0("L", 1:5))) %>%
  filter((proba_exp <= 0.1 | proba_exp >= 0.9) & sample_prop >= 0.01) %>% print(n=110)
source("kbsd.R")
set.seed(15082025)
o5 <- data1_cat[-1]
o5_1 <- o5
o5_1$A <- 1
o5_2 <- o5
o5_2$A <- 0
res5 <- kbsd(data = o5,
             int_data_list = list(o5_1, o5_2),
             disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3=sd(o5$L3), 
                            L4 = sd(o5$L4), L5 = sd(o5$L5), A=0.5*sd(o5$A)),
             plot.out = F)
res5_plot <- kbsd(data = o5,
                  int_data_list = list(o5_1, o5_2),
                  disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3=sd(o5$L3),
                                 L4 = sd(o5$L4), L5 = sd(o5$L5), A=0.5*sd(o5$A)))
res5_plot
# check what strata have low support for IV=1
shift1 <- res5[res5$shift == 1,]
outliers1 <- shift1$diagnostic < quantile(shift1$diagnostic, probs = 0.25)  # create indices for the "outliers"
l_values1 <- data1_cat[outliers1, "L3"]   # original cat L3 values
mfv(l_values1)  # most with few support for IV=1 are from "(5,6]" -> poss bc viol is for L3<6

# check where few support for IV=2 (A=0)
shift2 <- res5[res5$shift == 2,]
outliers2 <- shift2$diagnostic < quantile(shift2$diagnostic, probs = 0.25)
l_values2 <- data1_cat[outliers2, "L3"]
mfv(l_values2)  # least support for IV=1 is for L3="(7,8]" (not expected, not planned)


# identify obs with low EDP
indices_outliers <- unique(res5[res5$shift == 1 & res5$diagnostic < quantile(res5$diagnostic, 0.05), "observation"])
o5[indices_outliers, ] %>% nrow()
o5[indices_outliers, ] %>% filter(L3 < 7) %>% nrow()
# q=0.05: 75/87 from critical stratum -> majority of obs reported as critical are from viol stratum!
# q=0.01: 18/19 from critical stratum!
o5[indices_outliers, ] %>% filter(!(L3 < 7))
# q=0.05: quite a few viol which prob alr flagged initially by make_strata_table



# 5 Confounders Correlated ----

set.seed(28092025)
sem1 <- DAG.empty() +
  node("L1", distr = "rbern", prob = 0.1) +
  node("L2", distr = "rbern", prob = 0.2) +
  node("L3", distr = "rnorm", mean = 6, sd = 1) +
  node("L4", distr = "rbern", prob = 0.1) +
  node("L5", distr = "rbern", prob = 0.2) +
  node("A", distr = "rbern", prob = plogis(0.2*L1 + 0.3*L2 + 0.1*L4 + 0.3*L5 - 2*L3*(L3 < 6)))
# all conf L_i=1 is v unlikely -> treatment not too likely, too but not extreme either
dag1 <- set.DAG(sem1)
plotDAG(dag1)
data1 <- sim(dag1, n = 1000, rndseed = 28092025)

data1$L2 <- data1$L1  # perfect corr between these 2

# check for extreme probs with cont conf L3:
data1 %>% filter(L3<6 & A==1) %>% nrow()/data1 %>% filter(L3<6) %>% nrow()  # viol: P(A=1)~0 for g >=1, any a & any b as sample prop = 52%

# check via table for extreme probs for all combos of binary conf: def in setup.R
binary_vars <- paste0("L", c(1,2,4,5))
tab <- make_strata_table(data1, A = "A", binary_vars = binary_vars)
tab %>% filter((proba_exp <= 0.1 | proba_exp >= 0.9) & sample_prop >= 0.01)
# no other viol among the binary vars expected

## PoRT: continuous var L3 uncategorised ----
source("data/port_utils.R")
lst5 <- list()
a_values <- c(0.01, 0.025, 0.05, 0.1)
b_values <- c(0.01, 5/(sqrt(nrow(data1))*log(nrow(data1))), 0.05, 0.1)
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
#sink("output_port/port_5_corr_uncat.txt")
# gamma = 1-5: always det


## PoRT: continuous var L3 categorised ----
source("data/port_utils.R")
data1_cat <- data1
data1_cat$L3 <- cut(data1_cat$L3, breaks = c(-Inf, 1, 2, 3, 4, 5, 6, 7, 8, Inf))
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
#sink("output_port/port_5_corr_cat.txt")
# g = 1-5: always det when should


## kbsd ----
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
res5_plot  # EDP range between 0-350, fewer support for IV=1
table(data1$A)  # which alr indicated here by fewer obs in A=1

# identify obs with low EDP: subgroup L3<6 with P(A=1)~0 should have few support for IV=1 (A=1)
indices_outliers <- unique(res5[res5$shift == 1 & res5$diagnostic < quantile(res5$diagnostic, 0.05), "observation"])
o5[indices_outliers, ] %>% nrow()
o5[indices_outliers, ] %>% filter(L3 <6) %>% nrow()
# q=0.01 & 0.05: majority of obs reported as critical are from viol stratum

# essence: both PoRT & kbsd det viol & for PoRT, the induced correlation here between the other vars 
# not involved in viol did not impact detection (as it did in 10 Conf scenario, but then again
# in 10 Conf scenario the problem with nondetection alr existed for uncorr setting so rather dimensionality problem?)



## KBSD cat ----
table(data1_cat$L3)
data1_cat <- data1_cat %>% 
  mutate(L3 = case_when(L3 == "(-Inf,1]" ~ 1, L3 == "(1,2]" ~ 2, L3 == "(2,3]" ~3, 
                        L3 == "(3,4]" ~4, L3 =="(4,5]" ~ 5, L3 == "(5,6]" ~6,
                        L3 == "(6,7]" ~7, L3 == "(7,8]"~8, L3 == "(8, Inf]"~9))
source("kbsd.R")
set.seed(15082025)
o5 <- data1_cat[-1]
o5_1 <- o5
o5_1$A <- 1
o5_2 <- o5
o5_2$A <- 0
res5 <- kbsd(data = o5,
             int_data_list = list(o5_1, o5_2),
             disthalf_vec=c(L1=0.5*sd(o5$L1), L2 = 0.5*sd(o5$L2), L3=0.5*sd(o5$L3), 
                            L4 = 0.5*sd(o5$L4), L5 = 0.5*sd(o5$L5), A=0.5*0.5*sd(o5$A)),
             plot.out = F)
res5_plot <- kbsd(data = o5,
                  int_data_list = list(o5_1, o5_2),
                  disthalf_vec=c(L1=0.5*sd(o5$L1), L2 = 0.5*sd(o5$L2), L3=0.5*sd(o5$L3),
                                 L4 = 0.5*sd(o5$L4), L5 = 0.5*sd(o5$L5), A=0.5*0.5*sd(o5$A)))
res5_plot

# identify obs with low EDP: should be for L3<7
indices_outliers <- unique(res5[res5$shift == 1 & res5$diagnostic < quantile(res5$diagnostic, 0.05), "observation"])
o5[indices_outliers, ] %>% nrow()
o5[indices_outliers, ] %>% filter(L3 <7) %>% nrow()
# q=0.01 & 0.05: majority of obs reported as critical are from viol stratum
o5[indices_outliers, ] %>% filter(!L3 <7)

# essence: majority of those with few EDP are indeed obs from the viol stratum!



# 5 Confounders With Internal Violation ----

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
  node("A", distr = "rbern", prob = plogis(0.2*L1 + 0.3*L2 + 0.1*L4 + 0.3*L5 - 0.5*L3*(L3 > 4 & L3 < 6)))
# all conf L_i=1 is v unlikely -> treatment not too likely, too but not extreme either
dag1 <- set.DAG(sem1)
plotDAG(dag1)
data1 <- sim(dag1, n = 1000, rndseed = 28092025)

plot(density(data1$L3), main = "bimodal L3 distribution from mixture") # again expect viol in [4,6]

# check for extreme probs with cont conf L3:
data1 %>% filter(L3 > 4 & L3 < 6 & A==1) %>% nrow()/
  data1 %>% filter(L3 > 4 & L3 < 6) %>% nrow()  # viol #1: P(A=1)~0, sample prop = 15.2% -> should find for g>=1, b=0.1, any a
# often in combo with L1=0/L2=0 (<- can be considered subviolations of viol #1)

# check via table for extreme probs for all combos of binary conf: def in setup.R
binary_vars <- paste0("L", c(1,2,4,5))
tab <- make_strata_table(data1, A = "A", binary_vars = binary_vars)
tab %>% print(n=80)
tab %>% filter((proba_exp <= 0.1 | proba_exp >= 0.9) & sample_prop >= 0.01)
# no extreme proba with only binary L_i expected (strata with p=0/1 are too small)

## PoRT: continuous var L3 uncategorised ----
source("data/port_utils.R")
lst5 <- list()
a_values <- c(0.01, 0.025, 0.05, 0.1)
b_values <- c(0.01, 5/(sqrt(nrow(data1))*log(nrow(data1))), 0.05, 0.1)
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
#sink("port_5_bimodal_uncat.txt")
# gamma = 1: detected , most informative for a=0.05/0.1 & b=0.1 (not too precise like for small a)
# gamma = 2-5: viol #1 detected, the other two with L1=2 or L2=0 only for a=0.05, b=gruber but subviol of #1 so does not matter
# viol #1 actually shows as "L3>=3.837 & L3< 5.996"


## PoRT: continuous var L3 categorised ----
source("data/port_utils.R")
data1_cat <- data1
data1_cat$L3 <- cut(data1_cat$L3, breaks = c(-Inf, 1, 2, 3, 4, 5, 6, 7, 8, Inf))
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
#sink("port_5_bimodal_cat.txt")
# gamma = 1,2: viol #1 detected for b=0.1 as should
# gamma = 2-3: viol #1 again for b=0.1, and for a =0.05, preciser subgroup "(5,6]"
# gamma = 4-5: due to cat, L3 groups are larger, i.e. as contrary to uncategorised case,
#              where many small viol groups with L3 could be split, L3 is now large so that
#              viol subgroups are created by intersecting with e.g. 2 other vars! imp obs 
#              (i.e. in cat, higher gamma values are "exploited"); viol #1 always det for b=0.1
# but so general trend is that for uncategorised case, small alpha means in case of 
# an existent viol subgroup, that subgroup is split into smaller strata to adhere to alpha
# -> for cat case, this is not poss anymore so that finding strata small enough that match the 
# small alpha is achieved by intersecting the categorised var with other vars (only poss for g>1 then tho)

# essence: viol #1 always detected when should (for b=0.1), both when uncat & cat, 
# but uncat even better bc larger subgroup found i.e. exceeds [45], [5,6]: lower bound is 3.2 actually


## kbsd ----
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
                            L4 = sd(o5$L4), L5 = sd(o5$L5), A=0.5*sd(o5$A)), 
             plot.out = F)
res5_plot <- kbsd(data = o5,
                  int_data_list = list(o5_1, o5_2),
                  disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3=sd(o5$L3),
                                 L4 = sd(o5$L4), L5 = sd(o5$L5), A=0.5*sd(o5$A))) 
res5_plot  # EDP range between 0-250, again fewer support for IV=1
table(data1$A)  # which alr indicated here by fewer obs in A=1

# identify obs with low EDP: subgroup L3=(4,6] with P(A=1)~0 should have few support for IV=1 (A=1)
indices_outliers <- unique(res5[res5$shift == 1 & res5$diagnostic < quantile(res5$diagnostic, 0.05), "observation"])
o5[indices_outliers, ] %>% nrow()
o5[indices_outliers, ] %>% filter(L3 >4 & L3 <6) %>% nrow()
# not even for q=0.01, majority of obs reported as critical are from viol stratum!
o5[indices_outliers, ] %>% filter(!(L3> 4 & L3 <6))

# essence: PoRT found the viol of P(A=1|L3 in [4,6])~0, kbsd didn't have it as majority among obs with low EDP



## KBSD cat ----
table(data1_cat$L3)
data1_cat <- data1_cat %>% 
  mutate(L3 = case_when(L3 == "(-Inf,1]" ~ 1, L3 == "(1,2]" ~ 2, L3 == "(2,3]" ~3, 
                        L3 == "(3,4]" ~4, L3 =="(4,5]" ~ 5, L3 == "(5,6]" ~6,
                        L3 == "(6,7]" ~7, L3 == "(7,8]"~8, L3 == "(8, Inf]"~9))
source("kbsd.R")
set.seed(15082025)
o5 <- data1_cat[-1]
o5_1 <- o5
o5_1$A <- 1
o5_2 <- o5
o5_2$A <- 0
res5 <- kbsd(data = o5,
             int_data_list = list(o5_1, o5_2),
             disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3=sd(o5$L3), 
                            L4 = sd(o5$L4), L5 = sd(o5$L5), A=0.5*sd(o5$A)),
             plot.out = F)
res5_plot <- kbsd(data = o5,
                  int_data_list = list(o5_1, o5_2),
                  disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3=sd(o5$L3),
                                 L4 = sd(o5$L4), L5 = sd(o5$L5), A=0.5*sd(o5$A)))
res5_plot

# identify obs with low EDP
indices_outliers <- unique(res5[res5$shift == 1 & res5$diagnostic < quantile(res5$diagnostic, 0.05), "observation"])
o5[indices_outliers, ] %>% nrow()
o5[indices_outliers, ] %>% filter(L3 > 4 & L3 < 6) %>% nrow()
# q=0.01: only 1 out of 15 obs reported as critical is from viol stratum!?
o5[indices_outliers, ] %>% filter(!(L3 > 4 & L3 <6)) # no real coherence & too small strata if standalone

# essence: viol strata were not the majority of obs with few EDP!






# 5 Confounders With External Violation ----

set.seed(28092025)
L3_1 <- rnorm(100, 3, 1)
L3_2 <- rnorm(900, 7, 1)
L3_1_2 <- c(L3_1, L3_2)
plot(density(L3_1_2))
sem1 <- DAG.empty() +
  node("L1", distr = "rbern", prob = 0.1) +
  node("L2", distr = "rbern", prob = 0.2) +
  node("L3", distr = "rconst", const = L3_1_2) +  # L3 follows a bimodal distribution,
  #        where vals around 5 have v low prob & vals to left & right have higher prob
  node("L4", distr = "rbern", prob = 0.1) +
  node("L5", distr = "rbern", prob = 0.2) +
  node("A", distr = "rbern", prob = plogis(0.2*L1 + 0.3*L2 + 0.1*L4 + 0.3*L5 - 2*L3*(L3 < 4)))
# all conf L_i=1 is v unlikely -> treatment not too likely, too but not extreme either
dag1 <- set.DAG(sem1)
plotDAG(dag1)
data1 <- sim(dag1, n = 1000, rndseed = 28092025)

plot(density(data1$L3), main = "bimodal L3 distribution from mixture") # again expect viol in [4,6]

# check for extreme probs with cont conf L3:
data1 %>% filter(L3<4 & A==1) %>% nrow()/data1 %>% filter(L3<4) %>% nrow()  # viol: for g >=1, all b, a<=0.05 as sample prop = 8.9%

# check via table for extreme probs for all combos of binary conf: def in setup.R
binary_vars <- paste0("L", c(1,2,4,5))
tab <- make_strata_table(data1, A = "A", binary_vars = binary_vars)
tab %>% print(n=80)
tab %>% filter((proba_exp <= 0.1 | proba_exp >= 0.9) & sample_prop >= 0.01)
# no extreme proba with only binary L_i expected (strata with p=0/1 are too small)

## PoRT: continuous var L3 uncategorised ----
source("data/port_utils.R")
lst5 <- list()
a_values <- c(0.01, 0.025, 0.05, 0.1)
b_values <- c(0.01, 5/(sqrt(nrow(data1))*log(nrow(data1))), 0.05, 0.1)
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
#sink("output_port/port_5_leftgap_uncat.txt")
# gamma = 1-5: always det when should, i.e. for beta >= gruber &  a <= 0.05


## PoRT: continuous var L3 categorised ----
source("data/port_utils.R")
data1_cat <- data1
data1_cat$L3 <- cut(data1_cat$L3, breaks = c(-Inf, 1, 2, 3, 4, 5, 6, 7, 8, Inf))
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
#sink("output_port/port_5_leftgap_cat.txt")
# g = 1-5: interestingly only returned [1-4] covered as interval as viol -> not [-Inf,4]!
#          the former makes sense for the small b=0.01, but the latter was should've been det for b>= gruber


## kbsd ----
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
                            L4 = sd(o5$L4), L5 = sd(o5$L5), A=0.5*sd(o5$A)),  # use 1 sd for L_i, 0.5 sd for A
             plot.out = F)
res5_plot <- kbsd(data = o5,
                  int_data_list = list(o5_1, o5_2),
                  disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3=sd(o5$L3),
                                 L4 = sd(o5$L4), L5 = sd(o5$L5), A=0.5*sd(o5$A)))  # use 1 sd for L_i, 0.5 sd for A
res5_plot  # EDP range between 0-200, fewer support for IV=2

# identify obs with low EDP: subgroup L3<4 with P(A=1)~0 should have few support for IV=1 (A=1)
indices_outliers <- unique(res5[res5$shift == 1 & res5$diagnostic < quantile(res5$diagnostic, 0.05), "observation"])
o5[indices_outliers, ] %>% nrow()
o5[indices_outliers, ] %>% filter(L3 <4) %>% nrow()

o5[indices_outliers, ] %>% filter(!(L3 <4)) # interesting that many with L4=1
data1 %>% filter(L4==1 & A==1) %>% nrow()/data1 %>% filter(L4==1) %>% nrow()  # but no viol

# essence: both PoRT & kbsd det viol that was a left side gap



## KBSD cat ----
table(data1_cat$L3)
data1_cat <- data1_cat %>% 
  mutate(L3 = case_when(L3 == "(-Inf,1]" ~ 1, L3 == "(1,2]" ~ 2, L3 == "(2,3]" ~3, 
                        L3 == "(3,4]" ~4, L3 =="(4,5]" ~ 5, L3 == "(5,6]" ~6,
                        L3 == "(6,7]" ~7, L3 == "(7,8]"~8, L3 == "(8, Inf]"~9))
source("kbsd.R")
set.seed(15082025)
o5 <- data1_cat[-1]
o5_1 <- o5
o5_1$A <- 1
o5_2 <- o5
o5_2$A <- 0
res5 <- kbsd(data = o5,
             int_data_list = list(o5_1, o5_2),
             disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3=sd(o5$L3), 
                            L4 = sd(o5$L4), L5 = sd(o5$L5), A=0.5*sd(o5$A)),
             plot.out = F)
res5_plot <- kbsd(data = o5,
                  int_data_list = list(o5_1, o5_2),
                  disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3=sd(o5$L3),
                                 L4 = sd(o5$L4), L5 = sd(o5$L5), A=0.5*sd(o5$A)))
res5_plot

# identify obs with low EDP
indices_outliers <- unique(res5[res5$shift == 1 & res5$diagnostic < quantile(res5$diagnostic, 0.05), "observation"])
o5[indices_outliers, ] %>% nrow()
o5[indices_outliers, ] %>% filter(L3<5) %>% nrow()
# only for q=0.01: majority of obs reported as critical are from viol stratum
o5[indices_outliers, ] %>% filter(!(L3<5)) # no viol bc strata just consist of 1 obs each

# essence: with checking for majority of flagged obs, viol stratum not well identified!



# 5 Confounders Binary ----

set.seed(28092025)
L3_1 <- rnorm(500, 3, 1)
L3_2 <- rnorm(500, 7, 1)
L3_1_2 <- c(L3_1, L3_2)
plot(density(L3_1_2))
sem1 <- DAG.empty() +
  node("L1", distr = "rbern", prob = 0.1) +
  node("L2", distr = "rbern", prob = 0.2) +
  node("L3", distr = "rbern", prob = 0.8) +
  node("L4", distr = "rbern", prob = 0.1) +
  node("L5", distr = "rbern", prob = 0.2) +
  node("A", distr = "rbern", prob = plogis(0.2*L1 + 0.3*L2 + 0.1*L4 + 0.3*L5 - 3*L3))
dag1 <- set.DAG(sem1)
plotDAG(dag1)
data1 <- sim(dag1, n = 1000, rndseed = 28092025)

# check for extreme probs with cont conf L3:
data1 %>% filter(L3==1 & A==1) %>% nrow()/data1 %>% filter(L3==1) %>% nrow()  # viol: P(A=1)~0 for g >=1, b>= 0.05, any a as sample prop = 80.7%

# check via table for extreme probs for all combos of binary conf: def in setup.R
binary_vars <- paste0("L", 1:5)
tab <- make_strata_table(data1, A = "A", binary_vars = binary_vars)
tab %>% filter((proba_exp <= 0.1 | proba_exp >= 0.9) & sample_prop >= 0.01) %>% print(n = 70)
# almost all involve L3=1 as wanted
# 8 others do not involve L3 at all!

## PoRT ----
source("data/port_utils.R")
a_values <- c(0.01, 0.025, 0.05, 0.1)
b_values <- c(0.01, 5/(sqrt(nrow(data1))*log(nrow(data1))), 0.05, 0.1)
g_values <- 1:5
lst5 <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst5[[paste0("gamma=", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A", cov.quanti = NULL, cov.quali = c("L1", "L2", "L3", "L4", "L5"),
             data = data1, alpha = a, beta = b, gamma = g)
    }
  }
}
lst5
#sink("output_port/port_5_binary.txt")
# g = 1-5: always det for any a & b=0.05/0.1 as should
# the higher gamma, the more intersections / smaller subgroups (usually with smaller b)


## kbsd ----
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
                            L4 = sd(o5$L4), L5 = sd(o5$L5), A=0.5*sd(o5$A)),  # use 1 sd for L_i, 0.5 sd for A
             plot.out = F)
res5_plot <- kbsd(data = o5,
                  int_data_list = list(o5_1, o5_2),
                  disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3=sd(o5$L3),
                                 L4 = sd(o5$L4), L5 = sd(o5$L5), A=0.5*sd(o5$A)))  # use 1 sd for L_i, 0.5 sd for A
res5_plot  # EDP range between 0-200, fewer support for IV=1
table(data1$A)  # which alr indicated here by fewer obs in A=1


# identify obs with low EDP: subgroup L3<4 with P(A=1)~0 should have few support for IV=1 (A=1)
indices_outliers <- unique(res5[res5$shift == 1 & res5$diagnostic < quantile(res5$diagnostic, 0.01), "observation"])
o5[indices_outliers, ] %>% nrow()
o5[indices_outliers, ] %>% filter(L3 ==1) %>% nrow()
# q=0.01 & 0.05: majority of obs reported as critical are from viol stratum
o5[indices_outliers, ] %>% filter(!(L3 ==1)) # no viol bc strata just consist of 1 obs each

# essence: both PoRT & kbsd det viol


