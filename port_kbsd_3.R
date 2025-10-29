source("setup.R")


# 3 Confounders ----

set.seed(15082025)
sem3 <- DAG.empty() +
  node("L1", distr = "rbern", prob = 0.3) +
  node("L2", distr = "rbern",prob = 0.5) + 
  node("L3", distr = "rnorm", mean = 0, sd = 1) +
  node("A", distr = "rbern", prob = 0.2*L1 + 0.3*L2 + 0.4*(L3<0))  # if L1=1,L2=1,L3<0 -> prob A=1 (viol #1) 
#                                                                    if L1=0, L2=0, L3 >0 -> prob A=0 (viol #2)
dag3 <- set.DAG(sem3)
plotDAG(dag3)
obs3 <- sim(dag3, rndseed = 12082025, n = 1000)
table(obs3$A)

obs3 %>% filter(L1==1 & L2==1 & L3 <0 & A==1) %>% nrow()/
  obs3 %>% filter(L1==1 & L2==1 & L3 <0) %>% nrow()  # viol #1: P(A=0|L)~0 <=> P(A=1|L)~1 for g=3, a<=0.05, b = 0.1

obs3 %>% filter(L1==0 & L2==0 & L3 >0 & A==1) %>% nrow()/
  obs3 %>% filter(L1==0 & L2==0 & L3 >0) %>% nrow()  # viol #2: P(A=1|L)~0 for g=3, any a bc sample prop = 18.7% & any b


## PoRT uncategorised ----
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
        port(A = "A", cov.quanti = "L3", cov.quali = c("L2", "L1"),
             data = obs3, alpha = a, beta = b, gamma = g) 
    }
  }
}
lst3
# g=3: viol #1 only det for a = 0.05 & b=0.1 but should've been for smaller a, too! 
#      for smaller a only viol with cont var L3. 
#      viol #2 for a=0.025 & b<=0.05, a=0.05 & b<=0.05, a = 0.1 & b<=0.05 -> covered for b=0.1
#      but undet for a=0.01 is not good! again suspicion that for a=0.01 only focus on small strata esp here with cont conf?

## PoRT: continuous var L3 categorised ----
obs3_cat <- obs3
obs3_cat$L3 <- cut(obs3_cat$L3, breaks = c(-Inf, -4, -3, -2, -1, 0, 1, 2, 3, 4, Inf))
lst3_cat <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst3_cat[[paste0("gamma=", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A", cov.quanti = NULL, cov.quali = c("L1", "L2", "L3"),
             data = obs3_cat, alpha = a, beta = b, gamma = g)
    }
  }
}
lst3_cat
# g=3: viol #1 and #2 now always det when should!

## KBSD ----
source("kbsd.R")
set.seed(15082025)
o3 <- obs3[-1]
o3_1 <- o3
o3_1$A <- 1
o3_2 <- o3
o3_2$A <- 0
res5 <- kbsd(data = o3,
             int_data_list = list(o3_1, o3_2),
             disthalf_vec=c(L1=0.5*sd(o3$L1), L2 = 0.5*sd(o3$L2), L3=0.5*sd(o3$L3), A=0.5*0.5*sd(o3$A)),  # use 1 0.5*sd for L_i, 0.5 0.5*sd for A
             plot.out = F)
res5_plot <- kbsd(data = o3,
                  int_data_list = list(o3_1, o3_2),
                  disthalf_vec=c(L1=0.5*sd(o3$L1), L2 = 0.5*sd(o3$L2), L3=0.5*sd(o3$L3), A=0.5*0.5*sd(o3$A)))  # use 1 0.5*sd for L_i, 0.5 0.5*sd for A
res5_plot
table(o3$A)  # see that alr fewer support among IV=1 bc fewer obs in A=1

# expect subgroup L1==0 & L2==0 & L3>0 to have few support for IV=1 (A=1) if would estimate Y|A=1 further:
shift1 <- res5[res5$shift == 1,]
outliers1 <- shift1$diagnostic < quantile(shift1$diagnostic, probs = 0.25)  # create indices for the "outliers"
l_values1 <- obs3[outliers1, "L3"] 
diag_values1 <- shift1[outliers1,]
plot(l_values1, diag_values1$diagnostic)  # most with few support for A=1 are indeed those with L3>0
subset_3 <- obs3[outliers1,]
table(subset_3[, c("L1", "L2")]) # most are indeed from L1=0 & L2=0!

# expect L1==1 & L2==1 & L3 <0 to have few support for IV=2 (A=0)
shift2 <- res5[res5$shift == 2,]
outliers2 <- shift2$diagnostic < quantile(shift2$diagnostic, probs = 0.25)
l_values2 <- obs3[outliers2, "L3"]  # original L3 values
diag_values2 <- shift1[outliers2,] # diag values
plot(l_values2, diag_values2$diagnostic)  # most with few support for A=0 are those with neg L3 values
subset_3 <- obs3[outliers2,]
table(subset_3[, c("L1", "L2")]) # most are indeed from L1=1 & L2=1

# essence: both viol det by PoRT & kbsd, but for PoRT uncat the critical strata were
#          not flagged for smaller a, so again suspicion with smaller a meaning PoRT only focuses on small strata
#          ! only consider viol #2 tho bc only that one also poss in correlated scenario!

# identify what strata the obs with low EDP belong to
indices_outliers <- res5[res5$shift == 1 & res5$diagnostic < quantile(res5$diagnostic, 0.05), "observation"]
indices_outliers <- unique(indices_outliers)
o3[indices_outliers,] %>% nrow()
o3[indices_outliers,] %>% filter(L1==1 & L2==1 & L3 <0) %>% nrow()
o3[indices_outliers,] %>% filter(L1==0 & L2==0 & L3 >0) %>% nrow()
# !!!: the smaller the quantile, the better can the critical strata be characterised! best for q=0.01
#      but better to know what to look for, esp for higher dim cases bc can verify viol preciser then
# check what strata the other obs are part of then: maybe detect new critical groups!
o3[indices_outliers,] %>% filter(!((L1==1 & L2==1 & L3 <0) | (L1==0 & L2==0 & L3 >0)))
# q=0.01: just one that is close to viol #2
# q=0.1: quite a few that have one of L1/L2 with other value
# this part above tells me how many from the critical group were correctly flagged by kbsd ###

# the part below searches for EDP value from which we say we have good support <-> when can construct CI ###
res5[res5$diagnostic < quantile(res5$diagnostic, 0.1),] %>% summary()
# being in lower q=0.01 corresponds to having max 16 obs close by/a diag value of 3.8 at most
# q=0.1: >30 obs for both viol strata necessitate >24.6 EDP -> for such obs could construct CI bc >30 other obs
# q=0.08: alr insufficient EDP in one of the treatment levels -> so for CI need EDP > EDP-8-percentile, i.e. diag >21.6

# exploratively: most with A=1 have L2=1 & L3 < -1 -> will be the obs that have few support for A=0 (IV=2)
# most with A=0 have L2=0 & L3 >1.2 -> will prob have few support for A=1 (IV=1)
# check P(A=a|L) for these strata now
o3 %>% filter(L2==1 & L3 < -1 & A==1) %>% nrow() / o3 %>% filter(L2 == 1 & L3 < -1) %>% nrow()  # almost all have A=1 -> viol for P(A=0|L)
o3 %>% filter(L2==0 & L3 > 1.2 & A==1) %>% nrow() / o3 %>% filter(L2==0 & L3 > 1.2) %>% nrow()  # almost all have A=0 -> viol for P(A=1|L)
# better to know what viol to look for: can verify viol 1+2 preciser then!
# L1==1 & L2==1 & L3 <0 and L1==0 & L2==0 & L3 >0

# mid-essence: above part shows prop of correctly flagged higher if choose lower quantile
#   whereas lower part shows we need a bit higher quantile to have sufficient obs to be able to construct CI
#   so use both bc interested in both! first checks how well kbsd performs, second gives researcher guidance
#   on threshold from when should be careful with trusting estimation of estimand 
#   -> to know how to proceed (e.g. obs under what threshold should be discarded) -> BUT NOT THAT EASY



## KBSD cat ----
table(obs3_cat$L3)
obs3_cat <- obs3_cat %>% 
  mutate(L3 = case_when(L3 == "(-Inf,-4]" ~ 1, L3 == "(-4,-3]" ~ 2, L3 == "(-3,-2]" ~3, 
                        L3 == "(-2,-1]" ~4, L3 =="(-1,0]" ~ 5, L3 == "(0,1]" ~6,
                        L3 == "(1,2]" ~7, L3 == "(2,3]"~8, L3 == "(3,4]"~9, L3 == "(4, Inf]"~10))
source("kbsd.R")
set.seed(15082025)
o3 <- obs3_cat[-1]
o3_1 <- o3
o3_1$A <- 1
o3_2 <- o3
o3_2$A <- 0
res5 <- kbsd(data = o3,
             int_data_list = list(o3_1, o3_2),
             disthalf_vec=c(L1=sd(o3$L1), L2 = sd(o3$L2), L3=sd(o3$L3), A=0.5*sd(o3$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
res5_plot <- kbsd(data = o3,
                  int_data_list = list(o3_1, o3_2),
                  disthalf_vec=c(L1=sd(o3$L1), L2 = sd(o3$L2), L3=sd(o3$L3), A=0.5*sd(o3$A)))  # use 1 SD for L_i, 0.5 SD for A
res5_plot
table(o3$A)

# check what strata have low support for IV=1
shift1 <- res5[res5$shift == 1,]
outliers1 <- shift1$diagnostic < quantile(shift1$diagnostic, probs = 0.25)  # create indices for the "outliers"
l_values1 <- obs3_cat[outliers1, "L3"]   # original cat L3 values
mfv(l_values1)  # most with few support for IV=1 are from "(1,2]" -> makes sense bc pos (even preciser than uncat, bc know more exactly where few obs!)
subset_3 <- obs3_cat[outliers1,]
table(subset_3[, c("L1", "L2")])  # among "(1,2]", most are from L1=1/0 & L2=0

# expect L1==1 & L2==1 & L3 <0 to have few support for IV=2 (A=0)
shift2 <- res5[res5$shift == 2,]
outliers2 <- shift2$diagnostic < quantile(shift2$diagnostic, probs = 0.25)
l_values2 <- obs3_cat[outliers2, "L3"]
mfv(l_values2)  # least support for IV=1 is for L3="(-1,0]" as expected 
subset_3 <- obs3_cat[outliers2,]
table(subset_3[, c("L1", "L2")]) # looking within L3="(-1,0]", most are indeed from L1=1 & L2=1

# subset strata with low support
indices_outliers <- unique(res5[res5$shift == 1 & res5$diagnostic < quantile(res5$diagnostic, 0.05), "observation"])
o3[indices_outliers,] %>% nrow()
o3[indices_outliers,] %>% filter(L1==1 & L2==1 & L3 %in% 1:5) %>% nrow()
o3[indices_outliers,] %>% filter(L1==0 & L2==0 & L3 %in% 6:10) %>% nrow()
# majority in q=0.05 is from viol stratum -> i.e. majority correctly identified
o3[indices_outliers,] %>% filter(!(L1==0 & L2==0 & L3 %in% 6:10)) # check how other obs with few EDP characterised



# 3 Confounders Correlated ----
set.seed(15082025)
sem3 <- DAG.empty() +
  node("L1", distr = "rbern", prob = 0.3) +
  node("L2", distr = "rbern", prob = L1) +   # perfect corr between L1 & L2
  node("L3", distr = "rnorm", mean = 0, sd = 1) +
  node("A", distr = "rbern", prob = 0.2*L1 + 0.3*L2 + 0.4*(L3<0))  # if L1=0, L2=0, L3 >0 (alr suffices to say L1=0 bc L1=L2 anyway) -> prob A=0
dag3 <- set.DAG(sem3)
plotDAG(dag3)
obs3 <- sim(dag3, rndseed = 12082025, n = 1000)
table(obs3$A)

obs3 %>% filter(L1==0 & L2==0 & L3 >0 & A==1) %>% nrow()/
  obs3 %>% filter(L1==0 & L2==0 & L3 >0) %>% nrow()  # expect viol P(A=1|L)~0 for g=3, any a bc sample prop = 36% & any b

## PoRT uncategorised ----
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
        port(A = "A", cov.quanti = "L3", cov.quali = c("L2", "L1"),
             data = obs3, alpha = a, beta = b, gamma = g) 
    }
  }
}
lst3
# g=3: det for all except a=0.01 & any b, a=0.025 & b=0.05, a=0.025/0.05 & b=0.1

## PoRT: continuous var L3 categorised ----
source("data/port_utils.R")
obs3_cat <- obs3
obs3_cat$L3 <- cut(obs3_cat$L3, breaks = c(-Inf, -4, -3, -2, -1, 0, 1, 2, 3, 4, Inf))
lst3_cat <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst3_cat[[paste0("gamma=", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A", cov.quanti = NULL, cov.quali = c("L1", "L2", "L3"),
             data = obs3_cat, alpha = a, beta = b, gamma = g)
    }
  }
}
lst3_cat
# g=3: always det

## KBSD ----
source("kbsd.R")
set.seed(15082025)
o3 <- obs3[-1]
o3_1 <- o3
o3_1$A <- 1
o3_2 <- o3
o3_2$A <- 0
res5 <- kbsd(data = o3,
             int_data_list = list(o3_1, o3_2),
             disthalf_vec=c(L1=sd(o3$L1), L2 = sd(o3$L2), L3=sd(o3$L3), A=0.5*sd(o3$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
res5_plot <- kbsd(data = o3,
                  int_data_list = list(o3_1, o3_2),
                  disthalf_vec=c(L1=sd(o3$L1), L2 = sd(o3$L2), L3=sd(o3$L3), A=0.5*sd(o3$A)))  # use 1 SD for L_i, 0.5 SD for A
res5_plot
table(o3$A)  # see that alr fewer support among IV=1 bc fewer obs in A=1

# expect subgroup L1==0 & L2==0 & L3>0 to have few support for IV=1 (A=1) if would estimate Y|A=1 further:
shift1 <- res5[res5$shift == 1,]
outliers1 <- shift1$diagnostic < quantile(shift1$diagnostic, probs = 0.25)  # create indices for the "outliers"
l_values1 <- obs3[outliers1, "L3"] 
diag_values1 <- shift1[outliers1,]
plot(l_values1, diag_values1$diagnostic)  # most with few support for A=1 are indeed those with L3>0
subset_3 <- obs3[outliers1,]
table(subset_3[, c("L1", "L2")]) # most are indeed from L1=0 & L2=0!

# no viol expected for IV=2 (A=0)
shift2 <- res5[res5$shift == 2,]
outliers2 <- shift2$diagnostic < quantile(shift2$diagnostic, probs = 0.25)
l_values2 <- obs3[outliers2, "L3"] 
diag_values2 <- shift1[outliers2,]
plot(l_values2, diag_values2$diagnostic)  # i.e. makes sense that no real pattern here
subset_3 <- obs3[outliers2,]
table(subset_3[, c("L1", "L2")]) # but most with few support for A=0 are from L1=1 & L2=1

# essence: PoRT & kbsd det the viol, but PoRT missed it again for uncat case -> perfect after cat though

# identify obs with low EDP
indices_outliers <- unique(res5[res5$shift == 1 & res5$diagnostic < quantile(res5$diagnostic, 0.05), "observation"])
o3[indices_outliers, ] %>% nrow()
o3[indices_outliers, ] %>% filter(L1==0 & L2==0 & L3 >0) %>% nrow()
# q=0.05: 69/94=73% of obs with few EDP are from critical stratum
# q=0.01: all are from critical stratum! -> no other potential viol strata
o3[indices_outliers, ] %>% filter(!(L1==0 & L2==0 & L3 >0))
# q=0.05. most of rest is L1==1 & L2==1 & L3 <0 -> check if really a new viol stratum!
o3 %>% filter(L1==1 & L2==1 & L3 <0 & A==1) %>% nrow()/o3 %>% filter(L1==1 & L2==1 & L3 <0) %>% nrow() # not a viol

# from q=0.02: >30 obs of the critical stratum are detected -> see to how many EDP that equates to
res5[res5$diagnostic < quantile(res5$diagnostic, 0.02),] %>% summary()



## KBSD cat ----
table(obs3_cat$L3)
obs3_cat <- obs3_cat %>% 
  mutate(L3 = case_when(L3 == "(-Inf,-4]" ~ 1, L3 == "(-4,-3]" ~ 2, L3 == "(-3,-2]" ~3, 
                        L3 == "(-2,-1]" ~4, L3 =="(-1,0]" ~ 5, L3 == "(0,1]" ~6,
                        L3 == "(1,2]" ~7, L3 == "(2,3]"~8, L3 == "(3,4]"~9, L3 == "(4, Inf]"~10))
source("kbsd.R")
set.seed(15082025)
o3 <- obs3_cat[-1]
o3_1 <- o3
o3_1$A <- 1
o3_2 <- o3
o3_2$A <- 0
res5 <- kbsd(data = o3,
             int_data_list = list(o3_1, o3_2),
             disthalf_vec=c(L1=sd(o3$L1), L2 = sd(o3$L2), L3=sd(o3$L3), A=0.5*sd(o3$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
res5_plot <- kbsd(data = o3,
                  int_data_list = list(o3_1, o3_2),
                  disthalf_vec=c(L1=sd(o3$L1), L2 = sd(o3$L2), L3=sd(o3$L3), A=0.5*sd(o3$A)))  # use 1 SD for L_i, 0.5 SD for A
res5_plot
table(o3$A)
# check what strata have low support for IV=1
shift1 <- res5[res5$shift == 1,]
outliers1 <- shift1$diagnostic < quantile(shift1$diagnostic, probs = 0.25)  # create indices for the "outliers"
l_values1 <- obs3_cat[outliers1, "L3"]   # original cat L3 values
mfv(l_values1)  # most with few support for IV=1 are from "(1,2]" -> makes sense bc pos (even preciser than uncat, bc know more exactly where few obs!)
subset_3 <- obs3_cat[outliers1,]
table(subset_3[, c("L1", "L2")])  # among "(1,2]", most are indeed from L1=0 & L2=0

# expect L1==1 & L2==1 & L3 <0 to have few support for IV=2 (A=0)
shift2 <- res5[res5$shift == 2,]
outliers2 <- shift2$diagnostic < quantile(shift2$diagnostic, probs = 0.25)
l_values2 <- obs3_cat[outliers2, "L3"]
mfv(l_values2)  # least support for IV=1 is for L3="(-1,0]" as expected 
subset_3 <- obs3_cat[outliers2,]
table(subset_3[, c("L1", "L2")]) # looking within L3="(-1,0]", most are indeed from L1=1 & L2=1

# identify obs with low EDP
indices_outliers <- unique(res5[res5$shift == 1 & res5$diagnostic < quantile(res5$diagnostic, 0.05), "observation"])
o3[indices_outliers, ] %>% nrow()
o3[indices_outliers, ] %>% filter(L1==0 & L2==0 & L3 %in% 6:10) %>% nrow()
# q=0.05: 12/19=63% of obs with few EDP are from critical stratum -> majority
o3[indices_outliers, ] %>% filter(!(L1==0 & L2==0 & L3 %in% 6:10))
# q=0.05. most of rest is L1==1 & L2==1 & L3=(-3,-2] -> check if really a new viol stratum!
o3 %>% filter(L1==1 & L2==1 & L3 ==3 & A==1) %>% nrow()/o3 %>% filter(L1==1 & L2==1 & L3 ==3) %>% nrow() # not a viol




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
  node("A", distr = "rbern", prob = plogis(L1 + L2 - 0.5*L3*(L3 > 4 & L3 < 6)))
dag1 <- set.DAG(sem1)
plotDAG(dag1)
data1 <- sim(dag1, n = 1000)

# visualise distribution of L3
plot(density(data1$L3), main = "bimodal L3 distribution from mixture")
# i.e. few values of L3 that are close to 5 -> will have viol there

data1 %>% filter(L3 > 4 & L3 < 6 & A==1) %>% nrow()/
  data1 %>% filter(L3 > 4 & L3 < 6) %>% nrow()  # viol #1: P(A=1)~0 if L3 in [4,6], sample prop = 14.1% -> should find for g>=1, b=0.1

# table to check for all combos of binary conf if there are any with extreme P(A): def in setup.R
binary_vars <- paste0("L", c(1,2))
tab <- make_strata_table(data1, A = "A", binary_vars = binary_vars)
tab  # no extreme proba only with L1, L2

# distr of A within L1 and L2 is not too extreme
table(data1$L1, data1$A)
table(data1$L2, data1$A)


## PoRT: continuous var L3 uncategorised ----
source("data/port_utils.R")
lst5 <- list()
a_values <- c(0.01, 0.025, 0.05, 0.1)
b_values <- c(0.01, 5/(sqrt(nrow(data1))*log(nrow(data1))), 0.05, 0.1)
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
# gamma =1-3: viol #1 always identified where appropriate, viol 2+3 often covered
# but wrong numbers for several subgroups? ex. g=1, a=0.01, b=0.1:
# sometimes proba_exp wrong, sometimes sample prop -> a bit imprecision due to rounding


## PoRT: continuous var L3 categorised ----
source("data/port_utils.R")
data1_cat <- data1
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
# gamma = 1, 2, 3: det viol #1 always and also only! with cat being ad broad, other viol not poss to cover anymore

# essence: overall viol #1 found in both cat & uncat, but when uncat,
#          can return preciser subgroup that includes L2=0, whereas after cat
#          the broader subgroup only is returned -> so emphasises importance of cat thresholds
#          appropriate/meaningful for the context as Chatton et al also emphasised -> to limit info loss


## KBSD ----
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

# plot with subplot "values of var with viol" & subplot "EDP boxplots"
internal_3 <- ggplot(data1, aes(x=factor(A), y=L3)) + 
  geom_point(size=3, alpha=0.1) +
  theme(axis.title.x = element_blank()) + 
  scale_x_discrete(labels = c("A = 0", "A = 1"))
internal_3 / res5_plot
#ggsave("internal_3_edp.png", width = 4, height = 9)

# expect subgroup L3=(4,6] to have few support for IV=1 (A=1) if would estimate Y|A=1 further:
shift1 <- res5[res5$shift == 1,]
outliers1 <- shift1$diagnostic < quantile(shift1$diagnostic, probs = 0.25)
l_values1 <- data1[outliers1, "L3"] 
diag_values1 <- shift1[outliers1,]
plot(l_values1, diag_values1$diagnostic)
# interesting: both for threshold =0.25 & 0.5 not clear that those in [4,6] have lowest EDP/support, 
#              i.e. no clear pattern except fewer support for tail values & why this pattern in the plot?

# also check for IV=2
shift2 <- res5[res5$shift == 2,]
outliers2 <- shift2$diagnostic < quantile(shift2$diagnostic, probs = 0.25)  # create indices for the "outliers"
l_values2 <- data1[outliers2, "L3"]  # original L3 values
diag_values2 <- shift1[outliers2,] # diag values
plot(l_values2, diag_values2$diagnostic)
# no clear pattern except less support for tails bc fewer obs there (is ok, did not have viol with P(A=0)~0)

# essence: kbsd cannot PINPOINT that for IV=1 there is low support for L3 =[4,6], i.e. the "hole" whereas PoRT could

# identify obs with low EDP
indices_outliers <- unique(res5[res5$shift == 1 & res5$diagnostic < quantile(res5$diagnostic, 0.05), "observation"])
o5[indices_outliers, ] %>% nrow()
o5[indices_outliers, ] %>% filter(L3 > 4 & L3 < 6) %>% nrow()
# q=0.1: 12/118=10% from critical stratum only.. not well det
# q=0.05: only 2 from critical stratum
# q=0.01: none from critical stratum!
o5[indices_outliers, ] %>% filter(!(L3 > 4 & L3 < 6))



## KBSD cat ----
table(data1_cat$L3)
data1_cat <- data1_cat %>% 
  mutate(L3 = case_when(L3 == "(-Inf,2]" ~ 1, L3 == "(2,4]" ~ 2, L3 == "(4,6]" ~3,
                        L3 == "(6,8]"~4, L3 == "(8, Inf]"~5))
# check if any other viol now that categorised
make_strata_table(data1_cat, A = "A", binary_vars = c(binary_vars, "L3")) %>%
  filter((proba_exp <=0.1 | proba_exp >= 0.9) & sample_prop >= 0.01) %>% arrange(L3) %>% print(n=52) # no

source("kbsd.R")
set.seed(15082025)
o3 <- data1_cat[-1]
o3_1 <- o3
o3_1$A <- 1
o3_2 <- o3
o3_2$A <- 0
res5 <- kbsd(data = o3,
             int_data_list = list(o3_1, o3_2),
             disthalf_vec=c(L1=sd(o3$L1), L2 = sd(o3$L2), L3=sd(o3$L3), A=0.5*sd(o3$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
res5_plot <- kbsd(data = o3,
                  int_data_list = list(o3_1, o3_2),
                  disthalf_vec=c(L1=sd(o3$L1), L2 = sd(o3$L2), L3=sd(o3$L3), A=0.5*sd(o3$A)))  # use 1 SD for L_i, 0.5 SD for A
res5_plot

# check what strata have low support for IV=1
shift1 <- res5[res5$shift == 1,]
outliers1 <- shift1$diagnostic < quantile(shift1$diagnostic, probs = 0.25)  # create indices for the "outliers"
l_values1 <- data1_cat[outliers1, "L3"]   # original cat L3 values
mfv(l_values1)  # most with few support for IV=1 are from "(6,8]" -> also imprecise! should be (4,6]

# check where few support for IV=2 (A=0)
shift2 <- res5[res5$shift == 2,]
outliers2 <- shift2$diagnostic < quantile(shift2$diagnostic, probs = 0.25)
l_values2 <- data1_cat[outliers2, "L3"]
mfv(l_values2)  # least support for IV=1 is for L3="(2,4]" (not expected, not planned)

# identify obs with low EDP
indices_outliers <- unique(res5[res5$shift == 1 & res5$diagnostic < quantile(res5$diagnostic, 0.05), "observation"])
o3[indices_outliers, ] %>% nrow()
o3[indices_outliers, ] %>% filter(L3 > 4 & L3 < 6) %>% nrow()
# q=0.1: 12/118=16% from critical stratum only.. not well det
# q=0.05: only 4 from critical stratum
# q=0.01: none from critical stratum!
o3[indices_outliers, ] %>% filter(!(L3 > 4 & L3 < 6))
# q=0.05: quite a lot of other strata flagged -> alr checked if viol with make_strata_table and weren't




# 3 Confounders With Left-Side-Gap ----

set.seed(27092025)
L3_1 <- rnorm(100, 3, 1)
L3_2 <- rnorm(900, 7, 1)
L3_1_2 <- c(L3_1, L3_2)
plot(density(L3_1_2))
sem1 <- DAG.empty() +
  node("L1", distr = "rbern", prob = 0.1) +
  node("L2", distr = "rbern", prob = 0.2) +
  node("L3", distr = "rconst", const = L3_1_2) +  # L3 follows a bimodal distribution,
  #        where vals around 5 have v low prob & vals to left & right have higher prob
  node("A", distr = "rbern", prob = plogis(L1 + L2 - 2*L3*(L3 < 4)))
dag1 <- set.DAG(sem1)
plotDAG(dag1)
data1 <- sim(dag1, n = 1000, rndseed = 12082025)

data1 %>% filter(L3 < 4 & A==1) %>% nrow()/data1 %>% filter(L3 < 4) %>% nrow() # viol for g>=1, b>=0.05, a <=0.05 (sample prop = 8.1%)

# table to check for all combos of binary conf if there are any with extreme P(A): def in setup.R
binary_vars <- paste0("L", c(1,2))
tab <- make_strata_table(data1, A = "A", binary_vars = binary_vars)
tab  # no viol solely with L1 & L2

## PoRT: continuous var L3 uncategorised ----
source("data/port_utils.R")
lst5 <- list()
a_values <- c(0.01, 0.025, 0.05, 0.1)
b_values <- c(0.01, 5/(sqrt(nrow(data1))*log(nrow(data1))), 0.05, 0.1)
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
# g=1-3: always found when should be


## PoRT: continuous var L3 categorised ----
source("data/port_utils.R")
data1_cat <- data1
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
# g=1-3: also always found when should be

## KBSD ----
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
table(data1$A)  # similar support makes sense bc as much support for treated as for untreated -> but of interest, if not so within particulat strata!

# plot with subplots "viol var" & "EDP"
external_3 <- ggplot(data1, aes(x=factor(A), y=L3)) +
  geom_point(size=3, alpha=0.1) +
  theme(axis.title.x = element_blank()) +
  scale_x_discrete(labels= c("A = 0", "A = 1"))
external_3/res5_plot
#ggsave("external_3_edp.png", width =4, height = 9)

# expect subgroup L3<4 to have few support for IV=1 (A=1) if would estimate Y|A=1 further:
shift1 <- res5[res5$shift == 1,]
outliers1 <- shift1$diagnostic < quantile(shift1$diagnostic, probs = 0.25)
l_values1 <- data1[outliers1, "L3"]
diag_values1 <- shift1[outliers1,]
plot(l_values1, diag_values1$diagnostic)  # indeed esp few EDP for L3<4 as expected! most with L3<4 have A=0 so that few support if intervene them on A=1 and check for neighbours

# also check for IV=2
shift2 <- res5[res5$shift == 2,]
outliers2 <- shift2$diagnostic < quantile(shift2$diagnostic, probs = 0.25)
l_values2 <- data1[outliers2, "L3"]
diag_values2 <- shift2[outliers2,]
plot(l_values2, diag_values2$diagnostic)  # no clear pattern, also no viol expected for A=0

# essence: viol det by both

# identify obs with low EDP
indices_outliers <- unique(res5[res5$shift == 1 & res5$diagnostic < quantile(res5$diagnostic, 0.05), "observation"])
o5[indices_outliers, ] %>% nrow()
o5[indices_outliers, ] %>% filter(L3<4) %>% nrow()
# q=0.05: all obs with few EDP are from critical stratum!
# so q=0.05 really suffices to capture viol group -> researcher's attention raised
o5[indices_outliers, ] %>% filter(!L3<4)


## KBSD cat ----
table(data1_cat$L3)
data1_cat <- data1_cat %>% 
  mutate(L3 = case_when(L3 == "(-Inf,2]" ~ 1, L3 == "(2,4]" ~ 2, L3 == "(4,6]" ~3,
                        L3 == "(6,8]"~4, L3 == "(8, Inf]"~5))
source("kbsd.R")
set.seed(15082025)
o3 <- data1_cat[-1]
o3_1 <- o3
o3_1$A <- 1
o3_2 <- o3
o3_2$A <- 0
res5 <- kbsd(data = o3,
             int_data_list = list(o3_1, o3_2),
             disthalf_vec=c(L1=sd(o3$L1), L2 =sd(o3$L2), L3=sd(o3$L3), A=0.5*sd(o3$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
res5_plot <- kbsd(data = o3,
                  int_data_list = list(o3_1, o3_2),
                  disthalf_vec=c(L1=sd(o3$L1), L2 =sd(o3$L2), L3=sd(o3$L3), A=0.5*sd(o3$A)))  # use 1 SD for L_i, 0.5 SD for A
res5_plot
table(o3$A)

# check what strata have low support for IV=1
shift1 <- res5[res5$shift == 1,]
outliers1 <- shift1$diagnostic < quantile(shift1$diagnostic, probs = 0.25)  # create indices for the "outliers"
l_values1 <- data1_cat[outliers1, "L3"]   # original cat L3 values
mfv(l_values1)  # most with few support for IV=1 are from "(2,4]" -> makes sense bc viol is for L3<4

# check where few support for IV=2 (A=0)
shift2 <- res5[res5$shift == 2,]
outliers2 <- shift2$diagnostic < quantile(shift2$diagnostic, probs = 0.25)
l_values2 <- data1_cat[outliers2, "L3"]
mfv(l_values2)  # least support for IV=1 is for L3="(6,8]" (not expected, not planned)

# identify obs with low EDP
indices_outliers <- unique(res5[res5$shift == 1 & res5$diagnostic < quantile(res5$diagnostic, 0.05), "observation"])
o3[indices_outliers, ] %>% nrow()
o3[indices_outliers, ] %>% filter(L3 <3) %>% nrow()
# q=0.05: perfectly captured the viol as above




# 3 Confounders Binary ----
set.seed(15082025)
sem3 <- DAG.empty() +
  node("L1", distr = "rbern", prob = 0.3) +
  node("L2", distr = "rbern", prob = 0.4) + 
  node("L3", distr = "rbern", prob = 0.2) +
  node("A", distr = "rbern", prob = plogis(2*L1+L2))
dag3 <- set.DAG(sem3)
plotDAG(dag3)
obs3 <- sim(dag3, rndseed = 12082025, n = 1000)
table(obs3$A)

obs3 %>% filter(L1==1 & L2==1 & A==1) %>% nrow()/
  obs3 %>% filter(L1==1 & L2==1) %>% nrow()  # viol for g>=2, all a, b=0.1 -> for A=0!
  

# table to check for all combos of binary conf if there are any with extreme P(A): def in setup.R
binary_vars <- paste0("L", 1:3)
tab <- make_strata_table(obs3, A = "A", binary_vars = binary_vars)
tab %>% filter((proba_exp <= 0.1 | proba_exp >= 0.9) & sample_prop >= 0.01)
# 3 are subviol of L1=1 & L2=1, last one not expected but not focus, too

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
# always found when should!


## KBSD ----
source("kbsd.R")
o3 <- obs3[-1]
o3_1 <- o3
o3_1$A <- 1
o3_2 <- o3
o3_2$A <- 0
res3 <- kbsd(data = o3,
             int_data_list = list(o3_1, o3_2),
             disthalf_vec=c(L1=0.5*sd(o3$L1), L2 = 0.5*sd(o3$L2), L3 = 0.5*sd(o3$L3), A=0.5*0.5*sd(o3$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
res3_plot <- kbsd(data = o3,
                  int_data_list = list(o3_1, o3_2),
                  disthalf_vec=c(L1=sd(o3$L1), L2 = sd(o3$L2), L3 = sd(o3$L3), A=0.5*sd(o3$A)))  # use 1 SD for L_i, 0.5 SD for A
res3_plot  # for A=0, majority has low support -> reflects that viol affects A=0


# identify obs with low EDP
# for A=0: should be L1==1 & L2==1
indices_outliers_0 <- unique(res3[res3$shift == 2 & res3$diagnostic < quantile(res3[res3$shift == 2, "diagnostic"], 0.05), "observation"])
o3[indices_outliers_0, ] %>% nrow()
o3[indices_outliers, ] %>% filter(L1==1 & L2==1) %>% nrow() # all from critical stratum!
# for A=1: no violation expected
indices_outliers <- unique(res3[res3$shift == 1 & res3$diagnostic < quantile(res3[res3$shift == 1, "diagnostic"], 0.05), "observation"])
o3[indices_outliers, ] %>% nrow()
o3[indices_outliers, ] %>% filter(L1==1 & L2==1) %>% nrow()  # interestingly also all from L1=1 & L2=1 although that's a viol for A=0?
res3 %>% arrange(diagnostic) %>% View()


