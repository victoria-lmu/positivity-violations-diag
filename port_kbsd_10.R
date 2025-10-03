source("setup.R")

# 4) 10 Confounders ----

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
# if one of them =0, then P(A)=0.5; if all =1, then P(A)=0.88
DAG <- DAG 
DAG <- set.DAG(DAG)
dat <- sim(DAG, n = 1000)
table(dat$A)  # balanced

# table to check for all combos of binary conf if there are any with extreme P(A)
make_strata_table <- function(dat, A = "A", binary_vars){
  
  # create all combos of binary_vars of length >= 1
  subsets <- unlist(lapply(1:length(binary_vars), function(x) combn(binary_vars, x, simplify = FALSE)),
                    recursive = FALSE)
  
  # summary stats for each combo
  results <- map(subsets, function(vars) {
    dat %>%
      group_by(across(all_of(vars))) %>%
      summarise(n          = n(),
                n_treat    = sum(.data[[A]] == 1),
                n_control  = sum(.data[[A]] == 0),
                proba_exp  = mean(.data[[A]] == 1),
                sample_prop = n()/nrow(dat),
                .groups = "drop") %>%
      mutate(vars = paste(vars, collapse = ","))
  })
  
  # bind all results together
  bind_rows(results) %>%
    arrange(vars, proba_exp)
}
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
#sink("port_10_uncat.txt")
# gamma = 1,2: irrelevant for defined critical stratum bc not yet poss to cover as intersection of 3
#  -> so other viol among cont conf, but only for a=0.01 (g=1) and a=0.01/0.02 (g=2) -> v small, prob meaningless viol
#  -> the smaller you allow the subgroups to be, the extremer the pos viol can be defined WITH CONT CONF
# gamma = 3-5: from now on strata by 3 conf, so should have L6=1 & L7=1 & L8=1
# -> included for a = 0.01/0.01/0.1 & b=0.1, not det for a = 0.025 & b=0.1
# -> included for a = 0.05/0.1 & b=0.05, not det for a=0.01/0.025 & b=0.05


## PoRT: continuous vars categorised ----
# check if problem of specifying order of cov.quanti persists with categorisation
dat_cat <- dat
for (i in names(dat_cat)[-1]) {
  if (typeof(dat_cat[[i]]) == "double") {
    dat_cat[[i]] <- cut(dat_cat[[i]], breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8))
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
#sink("port_10_cat.txt")
# gamma = 1: no critical subgroup
# gamma = 2: v small subgroups for a=0.01 only
# g= 3-5: det for a = 0.025/0.05/0.1 & b=0.1 & undet for a = 0.01, b = 0.1  -> weird that not det for a = 0.01
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

# bottomline: PoRT directly returns strata (can check if sensible in context -> don't need to know viol in advance),
#             for kbsd you should know where to look for viol, else explorative by tracing back the strata (good for overview tho)


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
# did not work with minval????????????????????????????????????????


# type = "harmonicmean" instead of default type = "Rfast"
res5_hm_plot <- kbsd(data = o5, int_data_list = list(o5_1, o5_2), type = "harmonicmean",
                     disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3 = sd(o5$L3), L4 = sd(o5$L4), L5 = sd(o5$L5),
                     L6 = sd(o5$L6), L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9), L10 = sd(o5$L10),
                     A=0.5*sd(o5$A)),  # use 1 SD for L_i, 0.5 SD for A
                     plot.out = T)
ggsave("kbsd_10_hm.png", width = 6, height = 3)
res5_hm <- kbsd(data = o5, int_data_list = list(o5_1, o5_2), type = "harmonicmean",
                     disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3 = sd(o5$L3), L4 = sd(o5$L4), L5 = sd(o5$L5),
                                    L6 = sd(o5$L6), L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9), L10 = sd(o5$L10),
                                    A=0.5*sd(o5$A)),  # use 1 SD for L_i, 0.5 SD for A
                     plot.out = F)
sink("kbsd_10_hm.txt")
print(res5_hm)
sink()

