library(simcausal)
library(tidyverse)

# 4) 10 Confounders ----

set.seed(22092025)
#  do not use for loop anymore: returned weird obs -> rnorm obs were all similar 
# for (i in 1:5) {
#   DAG <- DAG + node(paste0("L", i), distr = "rnorm", mean = i, sd = 1)
# }
# for (i in 6:10) {
#   DAG <- DAG + node(paste0("L", i), distr = "rbern", prob = 0.5)
# }

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
  node("L10", distr = "rbern", prob = 0.5) +
  node("A", distr = "rbern", prob = plogis(2*L6*L7*L8))  # A only depends on these conf;
# should help characterising problematic strata (nonlinear function of bernoulli L)
DAG <- DAG 
DAG <- set.DAG(DAG)
dat <- sim(DAG, n = 1000)
table(dat$A)  # good that not too imbalanced

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
# sample size of L8=1 & L6=1 & L7=1: 0.13, proba_exp = 123/136 = 0.904 i.e. to detect with any alpha & beta = 0.1


## PoRT: continuous vars uncategorised ----
source("data/port_utils.R")
lst5 <- list()
a_values <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.1)
b_values <- c(0.001, 0.01, 5/(sqrt(nrow(dat))*log(nrow(dat))), 0.05, 0.1)
g_values <- 1:10
for (a in a_values) {
  for (b in b_values) {
    for (g in g_values) {
      lst5[[paste0("gamma" = g, "\n", "alpha = ", a, ", beta = ", b)]] <-
        port(A = "A", cov.quanti = c("L1", "L2", "L3", "L4", "L5"),
             cov.quali = c("L6", "L7", "L8", "L9", "L10"), data = dat, alpha = a, beta = b, gamma = g)
    }
  }
}
lst5
sink(file = "port_10_uncat.txt")
print(lst5)
sink()
# gamma = 1,2: predefined critical stratum not yet poss to cover bc intersection of 3
#  -> so other viol among cont conf, but only for alpha=0.01/0.02
#     (where small groups poss -> reason: greedy categorisation esp for small alphas)
#  -> the smaller you alow the subgroups to be, the extremer the pos viol can be defined WITH CONT CONF
# gamma = 3-10: from now on strata by 3 conf, so should have L6=1 & L7=1 & L8=1
#  -> included! BUT never for a=0.02/0.03/0.04 (& beta=0.1)!, only for a=0.01/0.05/0.1


port("A", cov.quali = c("L6", "L7", "L8", "L9", "L10"),
     cov.quanti = c("L3", "L2"), data = dat, alpha = 0.01, beta = 0.1, gamma = 4)
# important finding: order of specifying CONTINUOUS covars as argument matters -> returns diff subgroups!  ~~~~~~~~~~ seems to work after all ??~~~~~~~~~~~
# e.g. gamma = 3:
#   diff subgroups for c("L1", "L2", "L3", "L4", "L5"), c("L3", "L1", "L2", "L4", "L5"), c("L2", "L3", "L1", "L4", "L5")
#   also for L2, L3 (L8=1 & L6=1 & L7=1 is undetected) and L3, L2 (detected)
# e.g. gamma = 4:
#   undetected if order in continuous vars is L2, L1 -> detected if order is L1, L2
#   same for L2, L3 (undetected) and L3, L2 (detected)
# how is this possible?

## PoRT: continuous vars categorised ----
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
                                                                   data = dat_cat, alpha = a, beta = b, gamma = 10) 
  }
}
lst5_cat
# gamma = 1: no critical subgroup
# gamma = 2: one critical subgroup for alpha = 0.01, but v small
# gamma = 3-8: now detected for all alpha, beta = 0.1 as wanted! shows importance of categorising cont conf!! also computationally faster
# gamma = 9, 10: for all alpha & beta = 0.1, except for alpha = 0.01 & beta = 0.1
# problem of changing order in cov.quanti irrelevant here bc now all vars in cov.quali due to categorisation



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

# what strata are those with low support for treated (IV = 1)
subset_5_1 <- res5[res5$diagnostic < median(res5[res5$shift == 1, "diagnostic"]) & res5$shift == 1,]
o5[subset_5_1$observation, ]
# did not built in viol here on purpose so would have to search in detail which underlying covar values

# what strata are those with low support for treated (IV = 0): should be L8=1 & L6=1 & L7=1
subset_5_2 <- res5[res5$diagnostic < median(res5[res5$shift == 2, "diagnostic"]) & res5$shift == 2,]
o5[subset_5_2$observation, ]
table(o5[subset_5_2$observation, ][, c("L6", "L7", "L8")]) # highest count among untreated is L8=1 & L6=1 & L7=1 -> as expected


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


### alternative metrics for disthalf vec ----
# prob similar results as for other metrics/plots
mad(o5$A)  # mad unsuitable
IQR(o5$A)  # IQR ok here, treatment distr not as imbalanced as in 2 Conf setting
summary(o5$A)


### formula for EDP for high-dim covar set ----

# type = "minval" instead of default type = "Rfast"
res5 <- kbsd(data = o5, int_data_list = list(o5_1, o5_2), type = "minval",
             disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3 = sd(o5$L3), L4 = sd(o5$L4), L5 = sd(o5$L5),
                            L6 = sd(o5$L6), L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9), L10 = sd(o5$L10),
                            A=0.5*sd(o5$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
# type = "harmonicmean" instead of default type = "Rfast"
res5 <- kbsd(data = o5, int_data_list = list(o5_1, o5_2), type = "harmonicmean",
             disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3 = sd(o5$L3), L4 = sd(o5$L4), L5 = sd(o5$L5),
                            L6 = sd(o5$L6), L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9), L10 = sd(o5$L10),
                            A=0.5*sd(o5$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)



