source("setup.R")

# 10 Confounders ----

set.seed(22092025)
#  do not use for loop anymore: returned weird obs -> rnorm obs were all similar 

# L1-L10 (mix of cont & binary)
DAG <- DAG.empty() +
  node("L1", distr = "rnorm", mean = 1, sd = 1) +  # sd=10 did not change det of viol above, just shifted other viol (i.e. gamma=1,2)
  node("L2", distr = "rnorm", mean = 2, sd = 1) +
  node("L3", distr = "rnorm", mean = 3, sd = 1) +
  node("L4", distr = "rnorm", mean = 4, sd = 1) +
  node("L5", distr = "rnorm", mean = 5, sd = 1) +
  node("L6", distr = "rbern", prob = 0.5) +
  node("L7", distr = "rbern", prob = 0.5) +
  node("L8", distr = "rbern", prob = 0.5) +
  node("L9", distr = "rbern", prob = 0.5) +
  node("L10", distr = "rbern", prob = 0.5) +
  node("A", distr = "rbern", prob = plogis(3*L6*L7*L8))  # A only depends on these conf;
# if one of them =0, then P(A)=0.5; if all =1, then P(A)=0.95
DAG <- DAG 
DAG <- set.DAG(DAG)
dat <- sim(DAG, n = 1000, rndseed = 23092025)
table(dat$A)  # balanced

dat %>% filter(L6==1 & L7==1 & L8==1 & A==1) %>% nrow()/
  dat %>% filter(L6==1 & L7==1 & L8==1) %>% nrow()  # to be det for g>=3, b=0.1, any a

# table to check for all combos of binary conf if there are any with extreme P(A): fun defined in setup.R
binary_vars <- paste0("L", 6:10) # again, L6...L10 are binary
tab <- make_strata_table(dat, A = "A", binary_vars = binary_vars)
tab %>% filter((proba_exp <=0.1 | proba_exp >= 0.9) & sample_prop >= 0.01)
# as wanted: only viol stratum is L8=1 & L6=1 & L7=1, optionally L9 & L10
# sample prop of L8=1 & L6=1 & L7=1: 13.6%, proba_exp = 130/136 = 0.956 i.e. to detect with any a & b = 0.05/0.1



## PoRT: continuous vars uncategorised ----
source("data/port_utils.R")
lst5 <- list()
a_values <- c(0.01, 0.025, 0.05, 0.1)
b_values <- c(0.01, 5/(sqrt(nrow(dat))*log(nrow(dat))), 0.05, 0.1)
g_values <- 1:5
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst5[[paste0("gamma=", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A", cov.quanti = c("L1", "L2", "L3", "L4", "L5"),
             cov.quali = c("L6", "L7", "L8", "L9", "L10"), data = dat, alpha = a, beta = b, gamma = g)
    }
  }
}
lst5
# sink("port_10_uncat.txt")
# lst5
# sink()
# gamma = 1,2: irrelevant for defined critical stratum bc not yet poss to cover as intersection of 3
#  -> so other viol among cont conf, but only for a=0.01 (g=1) and a=0.01/0.02 (g=2) -> v small, prob meaningless viol
#  -> the smaller you allow the subgroups to be, the extremer the pos viol can be defined WITH CONT CONF
# gamma = 3-5: from now on strata by 3 conf, so should have L6=1 & L7=1 & L8=1
# -> included for a = 0.01/0.05/0.1 & b=0.1, not det for a = 0.025 & b=0.1
# -> included for a = 0.05/0.1 & b=0.05, not det for a=0.01/0.025 & b=0.05


## PoRT: continuous vars categorised ----
# check if problem of specifying order of cov.quanti persists with categorisation
dat_cat <- dat
for (i in names(dat_cat)[-1]) {
  if (typeof(dat_cat[[i]]) == "double") {
    dat_cat[[i]] <- cut(dat_cat[[i]], breaks = c(-Inf, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, Inf))
  }
}
lst5_cat <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst5_cat[[paste0("gamma = ", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A", cov.quanti = NULL, cov.quali = c("L1", "L2", "L3", "L4", "L5", 
                                                       "L6", "L7", "L8", "L9", "L10"),
             data = dat_cat, alpha = a, beta = b, gamma = g) 
    }
  }
}
lst5_cat
# sink("port_10_cat.txt")
# lst5_cat
# sink()
# gamma = 1: no critical subgroup
# gamma = 2: v small subgroups for a=0.01 only
# g= 3-5: det for all a & b=0.05
#         undet for a = 0.01, b = 0.1  -> weird that not det for a = 0.01
#    maybe bc if cat, then a=0.01 (too small a) lets focus on small strata only??
#    but was also problem in 20_cof_UNCAT setting so not sure if tied to categorisation -> maybe randomness?
#         det for all a & b=0.05 tho -> so weird, bc b=0.1 allows for more, hypothesis here:
# it's as if PoRT only recognises as viol if the violating gorups proba_exposure is close to b/1-b
# alr more detected which shows importance of cat cont conf!! also computationally faster


## KBSD ----
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
table(o5$A)  # a few less obs for A=0 could've indicated that there'll be less support for IV=2 (A=0)

# what strata are those with low support for treated (IV = 1)
subset_5_1 <- res5[res5$diagnostic < median(res5[res5$shift == 1, "diagnostic"]) & res5$shift == 1,]
table(o5[subset_5_1$observation, c("L6", "L7", "L8")])
# those with few support in IV=1 are L6=0 & L7=0 & L8=0, NOT L6=1 & L7=1 & L8=1 -> new viol found!!!!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
table(o5[subset_5_1$observation, c("L9", "L10")])
# many from L9=1 & L10=0 that have low support in IV=1

# what strata are those with low support for treated (IV = 0): should be L8=1 & L6=1 & L7=1
subset_5_2 <- res5[res5$diagnostic < median(res5[res5$shift == 2, "diagnostic"]) & res5$shift == 2,]
table(o5[subset_5_2$observation, ][, c("L6", "L7", "L8")]) # highest count among IV=0 is L8=1 & L6=1 & L7=1 -> as expected

# identify what points have low support, with cat data!

# essence: PoRT directly returns strata (can check if sensible in context -> don't need to know viol in advance),
#          for kbsd you should know where to look for viol, else explorative by tracing back the strata (good for overview tho)


### alternative values/metrics for disthalf vec ----

# disthalf_vec: best case values, i.e. double values <=> kernel less steep
res5_plot <- kbsd(data = o5,
                  int_data_list = list(o5_1, o5_2),
                  disthalf_vec=c(L1=2*sd(o5$L1), L2 = 2*sd(o5$L2), L3 = 2*sd(o5$L3), L4 = 2*sd(o5$L4), L5 = 2*sd(o5$L5),
                                 L6 = 2*sd(o5$L6), L7 = 2*sd(o5$L7), L8 = 2*sd(o5$L8), L9 = 2*sd(o5$L9), L10 = 2*sd(o5$L10),
                                 A=2*0.5*sd(o5$A)))
res5_plot
# table for tracing back critical strata (IV=0) more balanced now -> L8=1 & L6=1 & L7=1 not particularly problematic

# " " : worst case values, i.e. halve values <=> kernel steeper
res5_plot <- kbsd(data = o5,
                  int_data_list = list(o5_1, o5_2),
                  disthalf_vec=c(L1=0.5*sd(o5$L1), L2 = 0.5*sd(o5$L2), L3 = 0.5*sd(o5$L3), L4 = 0.5*sd(o5$L4), L5 = 0.5*sd(o5$L5),
                                 L6 = 0.5*sd(o5$L6), L7 = 0.5*sd(o5$L7), L8 = 0.5*sd(o5$L8), L9 = 0.5*sd(o5$L9), L10 = 0.5*sd(o5$L10),
                                 A=0.5*0.5*sd(o5$A)))
res5_plot
# table for tracing back critical strata (IV=0) most imbalanced now -> prominence of L8=1 & L6=1 & L7=1 v striking

# over all disthalf_vec specifications: same general trend with better support for IV=1, but EDP values quite diff


# prob similar results as for other metrics/plots
mad(o5$A)  # mad unsuitable
IQR(o5$A)  # IQR ok here, treatment distr not as imbalanced as in 2 Conf setting
summary(o5$A)


### formula for EDP for high-dim covar set ----

# type = "minval" instead of default type = "Rfast"
res5_mv_plot <- kbsd(data = o5, int_data_list = list(o5_1, o5_2), type = "minval",
             minval_vec = rep(0.5, 11),
             disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3 = sd(o5$L3), L4 = sd(o5$L4), L5 = sd(o5$L5),
                            L6 = sd(o5$L6), L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9), L10 = sd(o5$L10),
                            A=0.5*sd(o5$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = T)
ggsave("kbsd_10_mv.png")
res5_mv <- kbsd(data = o5, int_data_list = list(o5_1, o5_2), type = "minval",
                minval_vec = rep(0.5, 11),
                disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3 = sd(o5$L3), L4 = sd(o5$L4), L5 = sd(o5$L5),
                               L6 = sd(o5$L6), L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9), L10 = sd(o5$L10),
                               A=0.5*sd(o5$A)),  # use 1 SD for L_i, 0.5 SD for A
                plot.out = F)

# sink("kbsd_10_mv.txt")
# print(res5_mv)
# sink()
# did not work with minval?


# type = "harmonicmean" instead of default type = "Rfast"
res5_hm_plot <- kbsd(data = o5, int_data_list = list(o5_1, o5_2), type = "harmonicmean",
                     disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3 = sd(o5$L3), L4 = sd(o5$L4), L5 = sd(o5$L5),
                     L6 = sd(o5$L6), L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9), L10 = sd(o5$L10),
                     A=0.5*sd(o5$A)),  # use 1 SD for L_i, 0.5 SD for A
                     plot.out = T)
#ggsave("kbsd_10_hm.png", width = 6, height = 3)

#### essence: same trend, but much higher discrepancy between EDP ranges of the two IV types: IV=1 is 20 EDPs higher than IV=2 ----


res5_hm <- kbsd(data = o5, int_data_list = list(o5_1, o5_2), type = "harmonicmean",
                     disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3 = sd(o5$L3), L4 = sd(o5$L4), L5 = sd(o5$L5),
                                    L6 = sd(o5$L6), L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9), L10 = sd(o5$L10),
                                    A=0.5*sd(o5$A)),  # use 1 SD for L_i, 0.5 SD for A
                     plot.out = F)
#sink("kbsd_10_hm.txt")


## KBSD cat ----
table(dat_cat$L1)
table(dat_cat$L2)  # all have same categorisation -> give numerical repr to compute kbsd values
dat_cat <- dat_cat %>%
  mutate(across(all_of(c("L1", "L2", "L3", "L4", "L5")), 
                ~ case_when(
                  . == "(-Inf,-2]" ~ 1, . == "(-2,-1]"    ~ 2, . == "(-1,0]"    ~ 3,
                  . == "(0,1]"     ~ 4, . == "(1,2]"     ~ 5, . == "(2,3]"     ~ 6,
                  . == "(3,4]"     ~ 7, . == "(4,5]"     ~ 8, . == "(5,6]"     ~ 9,
                  . == "(6,7]"     ~ 10, . == "(7,8]"    ~ 11, . == "(8, Inf]"  ~ 12)))
table(dat_cat$L1)
table(dat_cat$L2)
source("kbsd.R")
o5 <- dat_cat[-1]
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
table(o5$A)  # a few less obs for A=0 could've indicated that there'll be less support for IV=2 (A=0)

# what strata are those with low support for treated (IV = 1)
subset_5_1 <- res5[res5$diagnostic < median(res5[res5$shift == 1, "diagnostic"]) & res5$shift == 1,]
table(o5[subset_5_1$observation, c("L1", "L2", "L3")])
table(o5[subset_5_1$observation, c("L3", "L4", "L5")])
table(o5[subset_5_1$observation, c("L6", "L7", "L8")])
table(o5[subset_5_1$observation, c("L9", "L10")]) # no viol expected; most from L9=1 & L10=0 tho that have low support in IV=1

# what strata are those with low support for treated (IV = 0): should be L8=1 & L6=1 & L7=1
subset_5_2 <- res5[res5$diagnostic < median(res5[res5$shift == 2, "diagnostic"]) & res5$shift == 2,]
table(o5[subset_5_2$observation, ][, c("L6", "L7", "L8")]) # highest count among IV=0 is L8=1 & L6=1 & L7=1 -> as expected

# identify obs with few support
indices_outliers <- unique(res5[res5$shift == 1 & 
                                  res5$diagnostic < quantile(res5[res5$shift ==1, "diagnostic"], 0.05), "observation"])
o5[indices_outliers, ] # exploratively -> v diff now to do manually 

indices_outliers0 <- unique(res5[res5$shift == 2 & 
                                   res5$diagnostic < quantile(res5[res5$shift ==2, "diagnostic"], 0.05), "observation"])
o5[indices_outliers0, ] %>% nrow()
o5[indices_outliers0, ] %>% filter(L6==1 & L7==1 & L8==1) %>% nrow()  # majority det as from the critical stratum
o5[indices_outliers0, ] %>% filter(!(L6==1 & L7==1 & L8==1))  # 21 remaining obs -> also diff to recognise strata manually


# essence: both for uncat & for cat data, kbsd correctly flagged the critical stratum as with low EDP



# 10 Confounders Binary ----

setwd("C:/Users/victo/OneDrive/Desktop/Uni/Statistik/BA_thesis/positivity-violations-diag")
set.seed(03102025)
DAG <- DAG.empty() +
  node("L1", distr = "rbern", prob = 0.5) +
  node("L2", distr = "rbern", prob = 0.5) +
  node("L3", distr = "rbern", prob = 0.5) +
  node("L4", distr = "rbern", prob = 0.5) +
  node("L5", distr = "rbern", prob = 0.5) +
  node("L6", distr = "rbern", prob = 0.5) +
  node("L7", distr = "rbern", prob = 0.5) +
  node("L8", distr = "rbern", prob = 0.5) +
  node("L9", distr = "rbern", prob = 0.5) +
  node("L10", distr = "rbern", prob = 0.5) +
  node("A", distr = "rbern", prob = plogis(3*L6*L7*L8))  # A only depends on these conf;
# if one of them =0, then P(A)=0.5; if all =1, then P(A=1)=0.95 -> viol for A=0

DAG <- DAG 
DAG <- set.DAG(DAG)
dat <- sim(DAG, n = 1000, rndseed = 03102025)
table(dat$A)  # balanced

dat %>% filter(L6==1 & L7==1 & L8==1 & A==1) %>% nrow()/
  dat %>% filter(L6==1 & L7==1 & L8==1) %>% nrow()  # to be det for g>=3, b=0.1, any a

# table to check for all combos of binary conf if there are any with extreme P(A)
binary_vars <- paste0("L", 1:10)
tab <- make_strata_table(dat, A = "A", binary_vars = binary_vars)
tab %>% 
  filter((proba_exp <=0.1 | proba_exp >= 0.9) & sample_prop >= 0.01) %>%
  arrange(desc(sample_prop)) %>% print(n=830)
# many viol with >5 covariates for stratum -> only interested in <= 5
tab %>% 
  filter((proba_exp <=0.1 | proba_exp >= 0.9) & sample_prop >= 0.01) %>%
  arrange(desc(sample_prop)) %>% print(n=830) %>% head(30) %>% print(n=50)
# as wanted: L6=1 & L7=1 & L8=1 with P(A=1)= 0.921, sample prop = 14% -> should be det for b=0.1 & any a

## PoRT ----
source("data/port_utils.R")
lst5 <- list()
a_values <- c(0.01, 0.025, 0.05, 0.1)
b_values <- c(0.01, 5/(sqrt(nrow(dat))*log(nrow(dat))), 0.05, 0.1)
g_values <- 1:5
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst5[[paste0("gamma=", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A", cov.quanti = NULL,
             cov.quali = c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10"), data = dat, alpha = a, beta = b, gamma = g)
    }
  }
}
lst5
# sink("port_10_binary.txt")
# lst5
# sink()
# g=1,2: no viol bc not yet inteserction of 3 poss
# g=3: always det for any a & b=0.1
# g=4: always det for any a & b=0.1, with additional L_i for any a & b=0.05
# g=5: always det for any a & b=0.1, with additional L_i for small a <= 0.025 & b=0.05


## KBSD ----
source("kbsd.R")
set.seed(03102025)
o5 <- dat[-1]
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
res5_plot  # again fewer support for IV=1 & NB: lower EDP range with 0-100 bc more dims <=> more diff for obs to be close
table(dat$A)  # which alr indicated here by fewer obs in A=1


# what strata are those with low support for treated (IV = 1)
subset_5_1 <- res5[res5$diagnostic < median(res5[res5$shift == 1, "diagnostic"]) & res5$shift == 1,]
table(o5[subset_5_1$observation, c("L6", "L7", "L8")])
# again, those with few support in IV=1 are L6=0 & L7=0 & L8=0 -> new viol found?
# check prevalence in data:
dat %>% filter(L6==0 &L7==0 &L8==0 & A==1) %>% nrow()/
  dat %>% filter(L6==0 &L7==0 &L8==0) %>% nrow()  # not extreme tho so actually no viol
table(o5[subset_5_1$observation, c("L9", "L10")])
# many from L9=1 & L10=0 that have low support in IV=1

# what strata are those with low support for treated (IV = 0): should be L8=1 & L6=1 & L7=1
subset_5_2 <- res5[res5$diagnostic < median(res5[res5$shift == 2, "diagnostic"]) & res5$shift == 2,]
table(o5[subset_5_2$observation, ][, c("L6", "L7", "L8")]) # highest count among IV=0 is L8=1 & L6=1 & L7=1 -> as expected


# identify obs with few support
indices_outliers <- unique(res5[res5$shift == 1 & 
                                  res5$diagnostic < quantile(res5[res5$shift ==1, "diagnostic"], 0.05), "observation"])
o5[indices_outliers, ] # exploratively -> v diff now to do manually 

indices_outliers0 <- unique(res5[res5$shift == 2 & 
                                   res5$diagnostic < quantile(res5[res5$shift ==2, "diagnostic"], 0.05), "observation"])
o5[indices_outliers0, ] %>% nrow()
o5[indices_outliers0, ] %>% filter(L6==1 & L7==1 & L8==1) %>% nrow()  # majority det as from the critical stratum
o5[indices_outliers0, ] %>% filter(!(L6==1 & L7==1 & L8==1))  # 5 remaining obs -> also diff to recognise strata manually


# essence: PoRT & kbsd detected viol -> binary case seems v agréable for both diags



# 10 Confounders With Left-Side-Gap ----

set.seed(04102025)
L3_1 <- rnorm(100, 3, 1)
L3_2 <- rnorm(900, 7, 1)
L3_1_2 <- c(L3_1, L3_2)
plot(density(L3_1_2))
sem1 <- DAG.empty() +
  node("L1", distr = "rbern", prob = 0.1) +
  node("L2", distr = "rbern", prob = 0.2) +
  node("L3", distr = "rconst", const = L3_1_2) +  # L3 follows a bimodal distribution with
  #        gap in the middle (i.e. values around 5 are rare & will receive A=1 v rarely),
  #        where vals around 5 have v low prob & vals to left & right have higher prob
  node("L4", distr = "rbern", prob = 0.1) +
  node("L5", distr = "rbern", prob = 0.2) +
  node("L6", distr = "rbern", prob = 0.1) +
  node("L7", distr = "rbern", prob = 0.2) +
  node("L8", distr = "rbern", prob = 0.1) +
  node("L9", distr = "rbern", prob = 0.2) +
  node("L10", distr = "rbern", prob = 0.1) +
  node("A", distr = "rbern", prob = plogis(0.2*L1 + 0.3*L2 + 0.1*L4 + 0.3*L5 - 2*L3*(L3 < 4)))
# treatment for L3~5 is highly unlikely
# all conf L_i=1 is v unlikely -> treatment not too likely, too but not extreme either
dag1 <- set.DAG(sem1)
plotDAG(dag1)
data1 <- sim(dag1, n = 1000, rndseed = 06102025)
table(data1$A)  # balanced

data1 %>% filter(L3 < 4 & A==1) %>% nrow()/
  data1 %>% filter(L3 < 4) %>% nrow()  # viol P(A=1)~0 for g>=1, b>=GRUBER (modified 20.10.!), a<=0.05 as sample prop =8.9%, maybe esp if all other L_i=0!

# table to check for all combos of binary conf if there are any with extreme P(A): fun defined in setup.R
binary_vars <- paste0("L", c(1,2,4:10))
tab <- make_strata_table(data1, A = "A", binary_vars = binary_vars)
tab %>% filter((proba_exp <= 0.1 | proba_exp >= 0.9) & sample_prop >= 0.01) %>% print(n=64)
# no viol with binary vars bc not large enough sample prop that could be relevant

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
        port(A = "A", cov.quanti = c("L3"),
             cov.quali = c("L1", "L2", "L4", "L5", "L6", "L7", "L8", "L9", "L10"),
             data = data1, alpha = a, beta = b, gamma = g)
    }
  }
}
lst5
# sink("port_10_leftgap_uncat.txt")
# lst5
# sink()
# same results for all g=1-5 as only viol for L3 expected: always det when should 
# tho most precise for a=0.05 & b=0.05
# over all g, if small a <= 0.025, many small strata of L3 bc cont confounder


## PoRT: continuous var L3 categorised ----
data1_cat <- data1
source("data/port_utils.R")
data1_cat$L3 <- cut(data1_cat$L3, breaks = c(-Inf, 2, 3, 4, 5, 6, 7, 8, Inf))
lst5_cat <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst5_cat[[paste0("gamma=", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A", cov.quanti = NULL,
             cov.quali = c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10"),
             data = data1_cat, alpha = a, beta = b, gamma = g)
    }
  }
}
lst5_cat
# sink("port_10_leftgap_cat.txt")
# lst5_cat
# sink()
# also same over all g=1-5, no viol ofr other binary conf as expected


## KBSD ----
source("kbsd.R")
set.seed(04102025)
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
res5_plot  # similar support, EDP range is 0-100
table(data1$A)  # which alr by treatment level distribution

# what strata are those with low support for treated (IV = 1)
shift1 <- res5[res5$shift == 1,]
outliers1 <- shift1$diagnostic < quantile(shift1$diagnostic, probs = .05)  # create indices for the "outliers"
l_values1 <- data1[outliers1, "L3"]
diag_values1 <- shift1[outliers1,]
plot(l_values1, diag_values1$diagnostic,
     xlab = "Values of Confounder L3", ylab = "EDP")
# plots with L_1 value vs EDP for other binary vars: no pattern expect that most
#   have L_i=0 (makes sense bc simulated with low P(A=1))
# not rly clear that there is few support for IV=1 (A=1) from L3<4, also have many with few support for [6,8] 
#   but at least visible for quantile=0.05 that smallest EDP are on the left side, so they fall into L3<4
# -> wanted to see that due to few obs from L3<4 with A=1, few support if we intervene these obs on IV=1 (A=1!)

# what strata are those with low support for untreated (IV = 0)
shift2 <- res5[res5$shift == 2,]
outliers2 <- shift2$diagnostic < quantile(shift2$diagnostic, probs = .25)
for (i in names(o5)[-11]) {
  l_values2 <- data1[outliers2, i]
  diag_values2 <- shift2[outliers2,]
  plot(l_values2, diag_values2$diagnostic, xlab = paste0("Values of Confounder ", i), ylab = "EDP")
}
# no real pattern: makes sense that support for A=0 is fine among L3<4 bc most in this stratum got A=0!
# worse for quantile = 0.5 than for 0.25 bc just overall fewer obs with L3<4 (see sim with low density on left)

# identify obs with few support
indices_outliers <- unique(res5[res5$shift == 1 & 
                                  res5$diagnostic < quantile(res5[res5$shift ==1, "diagnostic"], 0.05), "observation"])
o5[indices_outliers, ] %>% nrow() # exploratively
o5[indices_outliers, ] %>% filter(L3 <4) %>% nrow()

# essence: PoRT det viol both in uncat & cat data, tho cat is nicer to interpret if precise enough as here
#          kbsd could not pinpoint that underlying the low EDP there is our viol stratum, tho detection dep on choice of EDP threshold (quantile)..


## KBSD cat ----
table(data1_cat$L3)
data1_cat <- data1_cat %>% 
  mutate(L3 = case_when(L3 == "(-Inf,2]" ~ 1, L3 == "(2,3]" ~ 2, L3 == "(3,4]" ~3, 
                        L3 == "(4,5]" ~4, L3 =="(5,6]" ~ 5, L3 == "(6,7]" ~6,
                        L3 == "(7,8]"~7, L3 == "(8, Inf]"~8))
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
                            L4 = sd(o5$L4), L5 = sd(o5$L5), L6 = sd(o5$L6),
                            L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9),
                            L10 = sd(o5$L10), A=0.5*sd(o5$A)),
             plot.out = F)
res5_plot <- kbsd(data = o5,
                  int_data_list = list(o5_1, o5_2),
                  disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3=sd(o5$L3), 
                                 L4 = sd(o5$L4), L5 = sd(o5$L5), L6 = sd(o5$L6),
                                 L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9),
                                 L10 = sd(o5$L10), A=0.5*sd(o5$A)))
res5_plot
# check what strata have low support for IV=1
shift1 <- res5[res5$shift == 1,]
outliers1 <- shift1$diagnostic < quantile(shift1$diagnostic, probs = 0.25)  # create indices for the "outliers"
l_values1 <- data1_cat[outliers1, "L3"]   # original cat L3 values
mfv(l_values1)  # most with few support for IV=1 are from "(7,8]" -> makes no sense bc viol for <6!

# check where few support for IV=2 (A=0)
shift2 <- res5[res5$shift == 2,]
outliers2 <- shift2$diagnostic < quantile(shift2$diagnostic, probs = 0.25)
l_values2 <- data1_cat[outliers2, "L3"]
mfv(l_values2)  # not expected, not planned

# identify obs with few support
indices_outliers <- unique(res5[res5$shift == 1 & 
                                  res5$diagnostic < quantile(res5[res5$shift ==1, "diagnostic"], 0.05), "observation"])
o5[indices_outliers, ] %>% nrow() # exploratively
o5[indices_outliers, ] %>% filter(L3 <4) %>% nrow()



# 10 Confounders With Middle-Gap ----

set.seed(29092025)
L3_1 <- rnorm(500, 3, 1)
L3_2 <- rnorm(500, 7, 1)
L3_1_2 <- c(L3_1, L3_2)
plot(density(L3_1_2))
sem1 <- DAG.empty() +
  node("L1", distr = "rbern", prob = 0.1) +
  node("L2", distr = "rbern", prob = 0.2) +
  node("L3", distr = "rconst", const = L3_1_2) +  # L3 follows a bimodal distribution with
  #        gap in the middle (i.e. values around 5 are rare & will receive A=1 v rarely),
  #        where vals around 5 have v low prob & vals to left & right have higher prob
  node("L4", distr = "rbern", prob = 0.1) +
  node("L5", distr = "rbern", prob = 0.2) +
  node("L6", distr = "rbern", prob = 0.1) +
  node("L7", distr = "rbern", prob = 0.2) +
  node("L8", distr = "rbern", prob = 0.1) +
  node("L9", distr = "rbern", prob = 0.2) +
  node("L10", distr = "rbern", prob = 0.1) +
  node("A", distr = "rbern", prob = plogis(0.2*L1 + 0.3*L2 + 0.1*L4 + 0.3*L5 - 2*L3*(L3 > 4 & L3 < 6)))
# treatment for L3~5 is highly unlikely
# all conf L_i=1 is v unlikely -> treatment not too likely, too but not extreme either
dag1 <- set.DAG(sem1)
plotDAG(dag1)
data1 <- sim(dag1, n = 1000, rndseed = 30092025)

data1 %>% filter(L3> 4 &L3 < 6 & A==1) %>% nrow()/
  data1 %>% filter(L3> 4 &L3 < 6) %>% nrow()  # viol P(A=1)~0 for g>=1, any b, any a bc sample prop = 16%, maybe esp if all other L_i=0!

# table to check for all combos of binary conf if there are any with extreme P(A): fun defined in setup.R
binary_vars <- paste0("L", c(1,2,4:10))
tab <- make_strata_table(data1, A = "A", binary_vars = binary_vars)
tab %>% filter((proba_exp <= 0.1 | proba_exp >= 0.9) & sample_prop >= 0.01) %>% print(n=64)
# all viol with binary conf only with v small a values -> will only be detected for a=0.01 & prob not v relevant pos viol

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
        port(A = "A", cov.quanti = c("L3"),
             cov.quali = c("L1", "L2", "L4", "L5", "L6", "L7", "L8", "L9", "L10"),
             data = data1, alpha = a, beta = b, gamma = g)
    }
  }
}
lst5
#sink("port_10_bimodal_uncat.txt")
# gamma = 1-5: always found for all a & b


## PoRT: continuous var L3 categorised ----
data1_cat <- data1
source("data/port_utils.R")
data1_cat$L3 <- cut(data1_cat$L3, breaks = c(-Inf, 2, 3, 4, 5, 6, 7, 8, Inf))
lst5_cat <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst5_cat[[paste0("gamma=", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A", cov.quanti = NULL,
             cov.quali = c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10"),
             data = data1_cat, alpha = a, beta = b, gamma = g)
    }
  }
}
lst5_cat
#sink("port_10_bimodal_cat.txt")
# gamma = 1-5: also always found for all a & b, as stratum [4,5], [5,6]


## KBSD ----
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
res5_plot  # fewer support for IV=1 & NB: lower EDP range with 0-100 bc more dims <=> more diff for obs to be close
table(data1$A)  # which alr indicated here by fewer obs in A=1


# acc to viol, few support for IV=1 (A=1) if would estimate Y|A=1 further, bc P(A=1)~0 for subgroup L3=(5,6]:
shift1 <- res5[res5$shift == 1,]
outliers1 <- shift1$diagnostic < quantile(shift1$diagnostic, probs = .25)  # create indices for the "outliers"
for (i in names(o5)[-11]) {
  l_values1 <- data1[outliers1, i]
  diag_values1 <- shift1[outliers1,]
  plot(l_values1, diag_values1$diagnostic, xlab = paste0("Values of Confounder ", i), ylab = "EDP")
}
# for quantile = 0.05: visible that overall, L3 has fewer obs around 5: shows that if we intervened all on A=1,
# those with L3 around 5 have fewer support which reflects that there just don't exist many obs with L3~5 & A=1 (even more visible for quantile = 0.25)
# no clear trend for any other confounder -> also didn't expect any viol for P(A=1)

# no viol for IV=2 (A=0) expected
shift2 <- res5[res5$shift == 2,]
outliers2 <- shift2$diagnostic < quantile(shift2$diagnostic, probs = 0.05)
l_values2 <- data1[outliers2, "L3"]  # original L3 values
diag_values2 <- shift1[outliers2,] # diag values
plot(l_values2, diag_values2$diagnostic)  # no low EDP for L3~5 at all, bc all with L3~5 got A=0 so there is lots of support for IV=2 (A=0)

# the more in center, the higher the support, except between 4-6 rarer
# also, overall v low EDP due to high dim, so try alternative EDP formulas below

# identify obs with few support
indices_outliers <- unique(res5[res5$shift == 1 & 
                                  res5$diagnostic < quantile(res5[res5$shift ==1, "diagnostic"], 0.05), "observation"])
o5[indices_outliers, ] %>% nrow() # exploratively
o5[indices_outliers, ] %>% filter(L3 > 3 & L3 < 6) %>% nrow()


## KBSD cat ----
table(data1_cat$L3)
data1_cat <- data1_cat %>% 
  mutate(L3 = case_when(L3 == "(-Inf,2]" ~ 1, L3=="(2,3]" ~2, 
                        L3 == "(3,4]" ~3, L3 =="(4,5]" ~ 4, L3 == "(5,6]" ~5,
                        L3 == "(6,7]" ~6, L3 == "(7,8]"~7, L3 == "(8, Inf]"~8))
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
                            L4 = sd(o5$L4), L5 = sd(o5$L5), L6 = sd(o5$L6),
                            L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9),
                            L10 = sd(o5$L10), A=0.5*sd(o5$A)),
             plot.out = F)
res5_plot <- kbsd(data = o5,
                  int_data_list = list(o5_1, o5_2),
                  disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3=sd(o5$L3), 
                                 L4 = sd(o5$L4), L5 = sd(o5$L5), L6 = sd(o5$L6),
                                 L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9),
                                 L10 = sd(o5$L10), A=0.5*sd(o5$A)))
res5_plot
# check what strata have low support for IV=1
shift1 <- res5[res5$shift == 1,]
outliers1 <- shift1$diagnostic < quantile(shift1$diagnostic, probs = 0.25)  # create indices for the "outliers"
l_values1 <- data1_cat[outliers1, "L3"]   # original cat L3 values
mfv(l_values1)  # most with few support for IV=1 are from "(6,7]" -> weird bc viol for 4<L3<6!

# check where few support for IV=2 (A=0)
shift2 <- res5[res5$shift == 2,]
outliers2 <- shift2$diagnostic < quantile(shift2$diagnostic, probs = 0.25)
l_values2 <- data1_cat[outliers2, "L3"]
mfv(l_values2)  # (2,3] not expected, not planned as critical obs for A=0


# identify obs with few support
indices_outliers <- unique(res5[res5$shift == 1 & 
                                  res5$diagnostic < quantile(res5[res5$shift ==1, "diagnostic"], 0.05), "observation"])
o5[indices_outliers, ] %>% nrow() # exploratively
o5[indices_outliers, ] %>% filter(L3 > 3 & L3 < 6) %>% nrow()

indices_outliers0 <- unique(res5[res5$shift == 2 & 
                                   res5$diagnostic < quantile(res5[res5$shift ==2, "diagnostic"], 0.05), "observation"])
o5[indices_outliers0, ]




# 10 Confounders Correlated ----

setwd("C:/Users/victo/OneDrive/Desktop/Uni/Statistik/BA_thesis/positivity-violations-diag")
set.seed(03102025)
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
  node("L10", distr = "rbern", prob = 0.5) +
  node("A", distr = "rbern", prob = plogis(3*L6*L7*L8)) 
# if one of these 3 conf =0, then P(A)=0.5; if all =1, then P(A)=0.95
DAG <- DAG 
DAG <- set.DAG(DAG)
dat <- sim(DAG, n = 1000, rndseed = 23092025)
table(dat$A)  # balanced
hist(dat$L2)

cor(dat) > 0.3  # rn no real correlation between L_i -> induce artificially:
dat$L2 <- 2*dat$L1 + rnorm(1000,0,1)
cor(dat$L2, dat$L1)  # cor of 0.88
hist(dat$L2)
dat$L10 <- dat$L9
cor(dat$L10, dat$L9)  # cor of 1   - but corr does not involve vars of viol.. stupid?

# table to check for all combos of binary conf if there are any with extreme P(A)
binary_vars <- paste0("L", 6:10)
tab <- make_strata_table(dat, A = "A", binary_vars = binary_vars)  # function defined in setup.R
tab %>% filter((proba_exp <=0.1 | proba_exp >= 0.9) & sample_prop >= 0.01)  # defo L6=1 & L7=1 & L8=1 as viol
dat %>% filter(L6==1 & L7==1 & L8==1 & A==1) %>% nrow()/
  dat %>% filter(L6==1 & L7==1 & L8==1) %>% nrow()  # P(A=1|L)~1 <=> P(A=0|L)~0 & to be det for b=0.1, any a as sample prop = 13%


## PoRT: continuous vars uncategorised ----
source("data/port_utils.R")
lst5 <- list()
a_values <- c(0.01, 0.025, 0.05, 0.1)
b_values <- c(0.01, 5/(sqrt(nrow(dat))*log(nrow(dat))), 0.05, 0.1)
g_values <- 1:5
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst5[[paste0("gamma=", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A", cov.quanti = c("L1", "L2", "L3", "L4", "L5"),
             cov.quali = c("L6", "L7", "L8", "L9", "L10"), data = dat, alpha = a, beta = b, gamma = g)
    }
  }
}
lst5
#sink("port_10_corr_uncat.txt")
# g = 1, 2: not yet poss to det
# g=3-5: only det from a=0.01 & b=0.1 (but should for any b!) -> suspicion from 10 Cf 
#  normal does not hold anymomre .> just random that PoRT misses for certain values then?



## PoRT: continuous vars categorised ----
source("data/port_utils.R")
dat_cat <- dat
for (i in names(dat_cat)[-1]) {
  if (typeof(dat_cat[[i]]) == "double") {
    dat_cat[[i]] <- cut(dat_cat[[i]], breaks = c(-Inf, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, Inf))
  }
}
lst5_cat <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst5_cat[[paste0("gamma=", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A", cov.quanti = NULL,
             cov.quali = c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10"),
             data = dat_cat, alpha = a, beta = b, gamma = g)
    }
  }
}
lst5_cat
# sink("port_10_corr_cat.txt")
# lst5_cat
# sink()
# g=1,2: not poss yet to det
# g = 3-5: never det for a = 0.01, only det for a = 0.025 & b=0.01/gruber, then again from a = 0.05 & any b

# essence: PoRT det safely for larger a & b but weird that not alr earlier!


## KBSD ----
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
table(o5$A)  # a few less obs for A=0 could've indicated that there'll be less support for IV=2 (A=0)

# what strata are those with low support for treated (IV = 1)
subset_5_1 <- res5[res5$diagnostic < quantile(res5[res5$shift == 1, "diagnostic"], 0.25) & res5$shift == 1,]
table(o5[subset_5_1$observation, c("L6", "L7", "L8")]) # no real pattern -> also did not expect viol here
table(o5[subset_5_1$observation, c("L9", "L10")])
# many from L9=0 & L10=0 or 1&1 that have low support in IV=1

# what strata are those with low support for treated (IV = 0): should be L8=1 & L6=1 & L7=1
subset_5_2 <- res5[res5$diagnostic < quantile(res5[res5$shift == 2, "diagnostic"], 0.25) & res5$shift == 2,]
table(o5[subset_5_2$observation, ][, c("L6", "L7", "L8")]) # highest count among IV=0 is L8=1 & L6=1 & L7=1 -> as expected

# identify obs with few support
indices_outliers <- unique(res5[res5$shift == 1 & 
                                  res5$diagnostic < quantile(res5[res5$shift ==1, "diagnostic"], 0.05), "observation"])
o5[indices_outliers, ] # exploratively

indices_outliers0 <- unique(res5[res5$shift == 2 & 
                                   res5$diagnostic < quantile(res5[res5$shift ==2, "diagnostic"], 0.05), "observation"])
o5[indices_outliers0, ] %>% nrow()
o5[indices_outliers0, ] %>% filter(L6==1 & L7==1 & L8==1) %>% nrow()


# essence: both diags det viol, although PoRT fails for smaller a, b values


## KBSD cat ----
table(dat_cat$L1) # give numerical repr to L1-L5 to compute kbsd values
dat_cat <- dat_cat %>%
  mutate(across(all_of(c("L1", "L2", "L3", "L4", "L5")), 
                ~ case_when(
                  . == "(-Inf,-2]" ~ 1, . == "(-2,-1]"    ~ 2, . == "(-1,0]"    ~ 3,
                  . == "(0,1]"     ~ 4, . == "(1,2]"     ~ 5, . == "(2,3]"     ~ 6,
                  . == "(3,4]"     ~ 7, . == "(4,5]"     ~ 8, . == "(5,6]"     ~ 9,
                  . == "(6,7]"     ~ 10, . == "(7,8]"    ~ 11, . == "(8, Inf]"  ~ 12)))
source("kbsd.R")
o5 <- dat_cat[-1]
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
res5_plot  # overall few EDP as high-dim adjustment set
table(o5$A)

# what strata are those with low support for treated (IV = 1)
subset_5_1 <- res5[res5$diagnostic < median(res5[res5$shift == 1, "diagnostic"]) & res5$shift == 1,]
table(o5[subset_5_1$observation, c("L6", "L7", "L8")])
table(o5[subset_5_1$observation, c("L9", "L10")]) # no viol expected; most from L9=0 & L10=1 tho that have low support in IV=1

# what strata are those with low support for treated (IV = 0): should be L8=1 & L6=1 & L7=1
subset_5_2 <- res5[res5$diagnostic < median(res5[res5$shift == 2, "diagnostic"]) & res5$shift == 2,]
table(o5[subset_5_2$observation, ][, c("L6", "L7", "L8")]) # highest count among IV=0 is L8=1 & L6=1 & L7=1 -> as expected


# identify obs with few support
indices_outliers <- unique(res5[res5$shift == 1 & 
                                  res5$diagnostic < quantile(res5[res5$shift ==1, "diagnostic"], 0.05), "observation"])
o5[indices_outliers, ] # exploratively

indices_outliers0 <- unique(res5[res5$shift == 2 & 
                                   res5$diagnostic < quantile(res5[res5$shift ==2, "diagnostic"], 0.05), "observation"])
o5[indices_outliers0, ] %>% nrow()
o5[indices_outliers0, ] %>% filter(L6==1 & L7==1 & L8==1) %>% nrow()


# essence: kbsd also identified the viol with cat data
