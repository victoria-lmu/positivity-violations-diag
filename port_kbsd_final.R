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

# expect L1==1 & L2==1 & L3 <0 to have few support for IV=2 (A=0)
shift2 <- res5[res5$shift == 2,]
outliers2 <- shift2$diagnostic < quantile(shift2$diagnostic, probs = 0.25)
l_values2 <- obs3[outliers2, "L3"]  # original L3 values
diag_values2 <- shift1[outliers2,] # diag values
plot(l_values2, diag_values2$diagnostic)  # most with few support for A=0 are those with neg L3 values
subset_3 <- obs3[outliers2,]
table(subset_3[, c("L1", "L2")]) # most are indeed from L1=1 & L2=1

# essence: both viol det by PoRT & kbsd, but for PoRT uncat the critical strata were
#          not flagged for smaller a, so again suspicion with smaller a meaning PoRT only focuese on small strata
#          ! only consider viol #2 tho bc only that one also poss in correlated scenario!


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
obs3_cat <- obs3
source("data/port_utils.R")
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



# 3 Confounders Binary ----
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


## KBSD ----
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
table(obs3[subset_5$observation, ][, c("L3")]) # most with few support among A=1 are from L3=0 as should (#1)
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

# essence: KBSD found the 4 viol among A=1,
#          but for A=0, just like PoRT failed to det viol #4 (L1==1 & L2==1 & L3==1),
#          instead flagged 2 new strata among A=0 that do not present problem actually



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
data1 <- sim(dag1, n = 1000)

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
l_values1 <- data1_cat[outliers1, "L3"]   # original cat L3 values
mfv(l_values1)  # most with few support for IV=1 are from "(2,4]" -> makes sense bc viol is for L3<4


# check where few support for IV=2 (A=0)
shift2 <- res5[res5$shift == 2,]
outliers2 <- shift2$diagnostic < quantile(shift2$diagnostic, probs = 0.25)
l_values2 <- data1_cat[outliers2, "L3"]
mfv(l_values2)  # least support for IV=1 is for L3="(6,8]" (not expected, not planned)



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

