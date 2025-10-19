source("setup.R")

# 20 Confounders ----
set.seed(23092025)
DAG20 <- DAG.empty() +
  node("L1", distr = "rnorm", mean = 1, sd = 1) +
  node("L2", distr = "rnorm", mean = 2, sd = 1) +
  node("L3", distr = "rnorm", mean = 3, sd = 1) +
  node("L4", distr = "rnorm", mean = 4, sd = 1) +
  node("L5", distr = "rnorm", mean = 5, sd = 1) +
  node("L6", distr = "rnorm", mean = 6, sd = 1) +
  node("L7", distr = "rnorm", mean = 7, sd = 1) +
  node("L8", distr = "rnorm", mean = 8, sd = 1) +
  node("L9", distr = "rnorm", mean = 9, sd = 1) +
  node("L10", distr = "rnorm", mean = 10, sd = 1) +
  node("L11", distr = "rbern", prob = 0.5) +
  node("L12", distr = "rbern", prob = 0.5) +
  node("L13", distr = "rbern", prob = 0.5) +
  node("L14", distr = "rbern", prob = 0.5) +
  node("L15", distr = "rbern", prob = 0.5) +
  node("L16", distr = "rbern", prob = 0.5) +
  node("L17", distr = "rbern", prob = 0.5) +
  node("L18", distr = "rbern", prob = 0.5) +
  node("L19", distr = "rbern", prob = 0.5) +
  node("L20", distr = "rbern", prob = 0.5) +
  node("A", distr = "rbern", prob = plogis(L1 - L2 + L5*L19*L20))
# if one of L19 and L20 = 0, then P(A)~0.26, if both =1 then P(A)~1
DAG20 <- set.DAG(DAG20)
dat20 <- sim(DAG20, rndseed = 12082025, n = 1000)[-1]
table(dat20$A)  # good that not too imbalanced
dat20 %>% filter(L19==1 & L20 == 1 & A==1) %>% nrow()/
  dat20 %>% filter(L19==1 & L20 == 1) %>% nrow()  # P(A=1)~1 -> to be det for g>=2, b=0.1, any a (sample prop= 22.7%)

# table to check for all combos of binary conf if there are any with extreme P(A): fun defined in setup.R
binary_vars <- paste0("L", c(1:10))
tab <- make_strata_table(dat20, A = "A", binary_vars = binary_vars)
tab %>% filter((proba_exp <= 0.1 | proba_exp >= 0.9) & sample_prop >= 0.01) %>% print(n=64)
# no viol among binary vars -> all too small sample prop


## PoRT: continuous vars uncategorised ----
source('data/port_utils.R')
gruber <- 5/(sqrt(nrow(dat2))*log(nrow(dat2)))
port("A", cov.quanti = c("L1","L2","L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10"),
     cov.quali = c("L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19", "L20"),
     data = dat20, alpha = 0.01, beta = 0.1, gamma = 3)
# no longer work with a=0.02/0.03/0.04!
# but g=3: undet for a=0.01/0.025/0.05, det for a= 0.1
# (g= 2,3: undet for all b = 0.1 & a=0.01/0.02/gruber/0.05 (too small alpha lets algo focus on smaller strata?), det: a=0.1, 0.04, 0.03)

a_values <- c(0.01, 0.025, 0.05, 0.1)
b_values <- c(0.01, 5/(sqrt(nrow(dat20))*log(nrow(dat20))), 0.05, 0.1)
g_values <- 1:20
lst20 <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst20[[paste0("gamma = ", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A",cov.quanti = c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10"),
             cov.quali = c("L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19", "L20"),
             data = dat20, alpha = a, beta = b, gamma = g)
    }
  }
}
lst20
#sink("port_20_uncat.txt")
# g = 1: not poss yet bc intersection of 2
# g = 2-5: only det for a=0.1 & b=0.1, for no smaller a & b=0.1 even though should!


## PoRT: continuous vars categorised ----

# straightforward (preciser) categorisation
dat20_cat <- dat20
for (i in names(dat20)[1:10]) {
  dat20_cat[[i]] <- cut(dat20[[i]], breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14))
}
dat20_cat
g_values <- 1:5
lst20_cat <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst20_cat[[paste0("gamma = ", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A", cov.quanti = NULL, cov.quali = c("L1", "L2", "L3", "L4", "L5",
                                                       "L6", "L7", "L8", "L9", "L10",
                                                       "L11", "L12", "L13", "L14", "L15",
                                                       "L16", "L17", "L18", "L19", "L20"),
             data = dat20_cat, alpha = a, beta = b, gamma = g)
    }
  }
}
lst20_cat
#sink("output_port/port_20_cat.txt")
# g = 1: only other viol with cont vars & only with v small a=0.01/0.025
# g = 2-5: det from a>=0.025 & b=0.1! alr shows pos eff of cat CONT conf on detection of viol involving CAT conf


## kbsd ----
source("kbsd.R")
o6 <- dat20
o6_1 <- o6
o6_1$A <- 1
o6_2 <- o6
o6_2$A <- 0
res6 <- kbsd(data = o6,
             int_data_list = list(o6_1, o6_2),
             disthalf_vec=c(L1=sd(o6$L1), L2 = sd(o6$L2), L3 = sd(o6$L3), L4 = sd(o6$L4), L5 = sd(o6$L5),
                            L6 = sd(o6$L6), L7 = sd(o6$L7), L8 = sd(o6$L8), L9 = sd(o6$L9), L10 = sd(o6$L10),
                            L11=sd(o6$L11), L12 = sd(o6$L12), L13 = sd(o6$L13), L14 = sd(o6$L14), L15 = sd(o6$L15),
                            L16 = sd(o6$L16), L17 = sd(o6$L17), L18 = sd(o6$L18), L19 = sd(o6$L19), L20 = sd(o6$L20),
                            A=0.5*sd(o6$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
res6_plot <- kbsd(data = o6,
             int_data_list = list(o6_1, o6_2),
             disthalf_vec=c(L1=sd(o6$L1), L2 = sd(o6$L2), L3 = sd(o6$L3), L4 = sd(o6$L4), L5 = sd(o6$L5),
                            L6 = sd(o6$L6), L7 = sd(o6$L7), L8 = sd(o6$L8), L9 = sd(o6$L9), L10 = sd(o6$L10),
                            L11=sd(o6$L11), L12 = sd(o6$L12), L13 = sd(o6$L13), L14 = sd(o6$L14), L15 = sd(o6$L15),
                            L16 = sd(o6$L16), L17 = sd(o6$L17), L18 = sd(o6$L18), L19 = sd(o6$L19), L20 = sd(o6$L20),
                            A=0.5*sd(o6$A)))
res6_plot
# overall v low EDP which makes sense as all obs are v spread out over 20 dims
# -> many defined strata so that within each prob just one treatment level -> no support for the other
#    was warned for in kbsd paper
# -> but interesting that among IV=1 (bringing all to A=1), particularly low support (many with EDP = 0)

# what strata are those with low support for treated (IV = 1)
subset_6_1 <- res6[res6$diagnostic <= median(res6[res6$shift == 1, "diagnostic"]) & res6$shift == 1,]  # many obs with v low EDP
table(o6[subset_6_1$observation, c("L19", "L20")])
# those with few support in IV=1 are L19=0 & L20=0

# what strata are those with low support for treated (IV = 0): should be L8=1 & L6=1 & L7=1
subset_6_2 <- res6[res6$diagnostic < median(res6[res6$shift == 2, "diagnostic"]) & res6$shift == 2,]
table(o6[subset_6_2$observation, c("L19", "L20")]) # highest count among IV=0 is L19=1 & L20=1 -> viol found

# essence: imp insight that cat of cont vars impacts detection of viol in CATEGORICAL vars!
#          -> here, for broader cat higher chance of det with all alpha > 0.01, without cat only for alpha <= 0.03
#          kbsd also flagged those from L19=1 & L20=1 as with few support for A=0 so overall det by both diags


### alternative EDP calc for high-dim covar set ----

# type = "minval" instead of default type = "Rfast"
# requires specification of EDP minvalue accepted for each covar (more imp <=> lower/no minval, less imp <=> 0.5)
res6 <- kbsd(data = o6, int_data_list = list(o6_1, o6_2), type = "minval",
             disthalf_vec=c(L1=sd(o6$L1), L2 = sd(o6$L2), L3 = sd(o6$L3), L4 = sd(o6$L4), L5 = sd(o6$L5),
                            L6 = sd(o6$L6), L7 = sd(o6$L7), L8 = sd(o6$L8), L9 = sd(o6$L9), L10 = sd(o6$L10),
                            L11=sd(o6$L11), L12 = sd(o6$L12), L13 = sd(o6$L13), L14 = sd(o6$L14), L15 = sd(o6$L15),
                            L16 = sd(o6$L16), L17 = sd(o6$L17), L18 = sd(o6$L18), L19 = sd(o6$L19), L20 = sd(o6$L20),
                            A=0.5*sd(o6$A)),  # use 1 SD for L_i, 0.5 SD for A
             minval_vec = , plot.out = T)
# type = "harmonicmean" instead of default type = "Rfast"
res6 <- kbsd(data = o6, int_data_list = list(o6_1, o6_2), type = "harmonicmean",
             disthalf_vec=c(L1=sd(o6$L1), L2 = sd(o6$L2), L3 = sd(o6$L3), L4 = sd(o6$L4), L5 = sd(o6$L5),
                            L6 = sd(o6$L6), L7 = sd(o6$L7), L8 = sd(o6$L8), L9 = sd(o6$L9), L10 = sd(o6$L10),
                            L11=sd(o6$L11), L12 = sd(o6$L12), L13 = sd(o6$L13), L14 = sd(o6$L14), L15 = sd(o6$L15),
                            L16 = sd(o6$L16), L17 = sd(o6$L17), L18 = sd(o6$L18), L19 = sd(o6$L19), L20 = sd(o6$L20),
                            A=0.5*sd(o6$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = T)


## KBSD cat ----
table(dat20_cat$L1)
table(dat20_cat$L2)  # all have same categorisation -> give numerical repr to compute kbsd values
dat20_cat <- dat20_cat %>%
  mutate(across(all_of(c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10")), 
                ~ case_when(
                  . == "(-3,-2]" ~ 1, . == "(-2,-1]"    ~ 2, . == "(-1,0]"    ~ 3,
                  . == "(0,1]"     ~ 4, . == "(1,2]"     ~ 5, . == "(2,3]"     ~ 6,
                  . == "(3,4]"     ~ 7, . == "(4,5]"     ~ 8, . == "(5,6]"     ~ 9,
                  . == "(6,7]"     ~ 10, . == "(7,8]"    ~ 11, . == "(8,9]"  ~ 12,
                  . == "(9,10]"     ~ 13, . == "(10,11]"    ~ 13, . == "(11,12]"  ~ 15,
                  . == "(12,13]"    ~ 16, . == "(13,14]"  ~ 17)))
source("kbsd.R")
o5 <- dat20_cat
o5_1 <- o5
o5_1$A <- 1
o5_2 <- o5
o5_2$A <- 0
res5 <- kbsd(data = o5,
             int_data_list = list(o5_1, o5_2),
             disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3 = sd(o5$L3), L4 = sd(o5$L4), L5 = sd(o5$L5),
                            L6 = sd(o5$L6), L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9), L10 = sd(o5$L10),
                            L11=sd(o5$L11), L12 = sd(o5$L12), L13 = sd(o5$L13), L14 = sd(o5$L14), L15 = sd(o5$L15),
                            L16 = sd(o5$L16), L17 = sd(o5$L17), L18 = sd(o5$L18), L19 = sd(o5$L19), L20 = sd(o5$L20),
                            A=0.5*sd(o5$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
res5_plot <- kbsd(data = o5,
                  int_data_list = list(o5_1, o5_2),
                  disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3 = sd(o5$L3), L4 = sd(o5$L4), L5 = sd(o5$L5),
                                 L6 = sd(o5$L6), L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9), L10 = sd(o5$L10),
                                 L11=sd(o5$L11), L12 = sd(o5$L12), L13 = sd(o5$L13), L14 = sd(o5$L14), L15 = sd(o5$L15),
                                 L16 = sd(o5$L16), L17 = sd(o5$L17), L18 = sd(o5$L18), L19 = sd(o5$L19), L20 = sd(o5$L20),
                                 A=0.5*sd(o5$A)))
res5_plot  # overall v few EDP as high-dim adjustment set
table(o5$A)  # a few less obs for A=1 matches that IV=1 has less support

# what strata are those with low support for treated (IV = 1)
subset_5_1 <- res5[res5$diagnostic < median(res5[res5$shift == 1, "diagnostic"]) & res5$shift == 1,]
table(o5[subset_5_1$observation, c("L19", "L20")]) # no viol expected; most from L9=0 & L10=0 tho that have low support in IV=1

# what strata are those with low support for treated (IV = 0): should be L19=1 & L20=1
subset_5_2 <- res5[res5$diagnostic < median(res5[res5$shift == 2, "diagnostic"]) & res5$shift == 2,]
table(o5[subset_5_2$observation, ][, c("L19", "L20")]) # highest count among IV=0 is L19=1 & L20=1 -> as expected

# essence: kbsd det viol also for cat data



# 20 Confounders Binary ----

set.seed(23092025)
DAG20 <- DAG.empty() +
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
  node("L11", distr = "rbern", prob = 0.5) +
  node("L12", distr = "rbern", prob = 0.5) +
  node("L13", distr = "rbern", prob = 0.5) +
  node("L14", distr = "rbern", prob = 0.5) +
  node("L15", distr = "rbern", prob = 0.5) +
  node("L16", distr = "rbern", prob = 0.5) +
  node("L17", distr = "rbern", prob = 0.5) +
  node("L18", distr = "rbern", prob = 0.5) +
  node("L19", distr = "rbern", prob = 0.5) +
  node("L20", distr = "rbern", prob = 0.5) +
  node("A", distr = "rbern", prob = plogis(L1 - L2 - 3*L19*L20))
# if one of L19 and L20 = 0, then P(A)~0.26, if both =1 then P(A)~0
DAG20 <- set.DAG(DAG20)
dat20 <- sim(DAG20, rndseed = 12082025, n = 1000)[-1]
table(dat20$A)  # a bit imbalanced
dat20 %>% filter(L19==1 & L20 == 1 & A==1) %>% nrow()/
  dat20 %>% filter(L19==1 & L20 == 1) %>% nrow()  # P(A=1)~0 -> to be det for g>=2, b=0.1, any a (sample prop= 24%)

# table to check for all the binary conf if there are any with extreme P(A): fun defined in setup.R
binary_vars <- paste0("L", c(1:20))
tab <- make_strata_table(dat20, A = "A", binary_vars = binary_vars)
tab %>% filter((proba_exp <= 0.1 | proba_exp >= 0.9) & sample_prop >= 0.01) %>% print(n=64)
# no viol among binary vars -> all too small sample prop


## PoRT ----
source('data/port_utils.R')
a_values <- c(0.01, 0.025, 0.05, 0.1)
b_values <- c(0.01, 5/(sqrt(nrow(dat20))*log(nrow(dat20))), 0.05, 0.1)
g_values <- 1:5
lst20 <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst20[[paste0("gamma = ", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A",cov.quanti = NULL,
             cov.quali = c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10",
                           "L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19", "L20"),
             data = dat20, alpha = a, beta = b, gamma = g)
    }
  }
}
lst20
#sink("port_20_binary.txt")
# g = 1: viol not poss yet bc intersection of 2
# g = 2-5: always det for any a & b=0.1, other large viol are subviol of L19=1 & L20=1 or smaller strata


## kbsd ----
source("kbsd.R")
o6 <- dat20
o6_1 <- o6
o6_1$A <- 1
o6_2 <- o6
o6_2$A <- 0
res6 <- kbsd(data = o6,
             int_data_list = list(o6_1, o6_2),
             disthalf_vec=c(L1=sd(o6$L1), L2 = sd(o6$L2), L3 = sd(o6$L3), L4 = sd(o6$L4), L5 = sd(o6$L5),
                            L6 = sd(o6$L6), L7 = sd(o6$L7), L8 = sd(o6$L8), L9 = sd(o6$L9), L10 = sd(o6$L10),
                            L11=sd(o6$L11), L12 = sd(o6$L12), L13 = sd(o6$L13), L14 = sd(o6$L14), L15 = sd(o6$L15),
                            L16 = sd(o6$L16), L17 = sd(o6$L17), L18 = sd(o6$L18), L19 = sd(o6$L19), L20 = sd(o6$L20),
                            A=0.5*sd(o6$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
res6_plot <- kbsd(data = o6,
                  int_data_list = list(o6_1, o6_2),
                  disthalf_vec=c(L1=sd(o6$L1), L2 = sd(o6$L2), L3 = sd(o6$L3), L4 = sd(o6$L4), L5 = sd(o6$L5),
                                 L6 = sd(o6$L6), L7 = sd(o6$L7), L8 = sd(o6$L8), L9 = sd(o6$L9), L10 = sd(o6$L10),
                                 L11=sd(o6$L11), L12 = sd(o6$L12), L13 = sd(o6$L13), L14 = sd(o6$L14), L15 = sd(o6$L15),
                                 L16 = sd(o6$L16), L17 = sd(o6$L17), L18 = sd(o6$L18), L19 = sd(o6$L19), L20 = sd(o6$L20),
                                 A=0.5*sd(o6$A)))

res6_plot
# overall v low EDP which makes sense as all obs are v spread out over 20 dims
# lower support for IV=1 (A=1) than for IV=2 (A=0)

# what strata are those with low support for treated (IV = 1): should find L19=1 & L20=1
subset_6_1 <- res6[res6$diagnostic <= quantile(res6[res6$shift == 1, "diagnostic"], 0.25) & res6$shift == 1,]  # many obs with v low EDP
table(o6[subset_6_1$observation, c("L19", "L20")])
# indeed most from those with few support in IV=1 are L19=1 & L20=1 -> viol found

# what strata are those with low support for treated (IV = 0):
subset_6_2 <- res6[res6$diagnostic < quantile(res6[res6$shift == 2, "diagnostic"], .25) & res6$shift == 2,]
o6[subset_6_2$observation, ]
table(o6[subset_6_2$observation, c("L19", "L20")]) # highest count among IV=0 is L19=1 & L20=0


# essence: both PoRT & kbsd det viol in case where all conf are binary




# 20 Confounders With Left-Side Gap ----

set.seed(23092025)
L5_1 <- rnorm(100, 3, 1)
L5_2 <- rnorm(900, 6, 1)
L5_1_2 <- c(L5_1, L5_2)
plot(density(L5_1_2))
DAG20 <- DAG.empty() +
  node("L1", distr = "rnorm", mean = 1, sd = 1) +
  node("L2", distr = "rnorm", mean = 2, sd = 1) +
  node("L3", distr = "rnorm", mean = 3, sd = 1) +
  node("L4", distr = "rnorm", mean = 4, sd = 1) +
  node("L5", distr = "rconst", const = L5_1_2) +
  node("L6", distr = "rnorm", mean = 6, sd = 1) +
  node("L7", distr = "rnorm", mean = 7, sd = 1) +
  node("L8", distr = "rnorm", mean = 8, sd = 1) +
  node("L9", distr = "rnorm", mean = 9, sd = 1) +
  node("L10", distr = "rnorm", mean = 10, sd = 1) +
  node("L11", distr = "rbern", prob = 0.5) +
  node("L12", distr = "rbern", prob = 0.5) +
  node("L13", distr = "rbern", prob = 0.5) +
  node("L14", distr = "rbern", prob = 0.5) +
  node("L15", distr = "rbern", prob = 0.5) +
  node("L16", distr = "rbern", prob = 0.5) +
  node("L17", distr = "rbern", prob = 0.5) +
  node("L18", distr = "rbern", prob = 0.5) +
  node("L19", distr = "rbern", prob = 0.5) +
  node("L20", distr = "rbern", prob = 0.5) +
  node("A", distr = "rbern", prob = plogis(L1 - L20 - 2*L5*(L5 < 4)))  # now there should be viol if L5 is v small
# if one of L19 and L20 = 0, then P(A)~0.26, if both =1 then P(A)~0
DAG20 <- set.DAG(DAG20)
dat20 <- sim(DAG20, rndseed = 12082025, n = 1000)[-1]
table(dat20$A)  # balanced
dat20 %>% filter(L5 <4 & A==1) %>% nrow()/
  dat20 %>% filter(L5 <4) %>% nrow()  # P(A=1)~0 -> to be det for g>=1, b=0.1, a <= 0.05 (sample prop= 9.9%)

# table to check for all combos of binary conf if there are any with extreme P(A): fun defined in setup.R
binary_vars <- paste0("L", c(6:10))
tab <- make_strata_table(dat20, A = "A", binary_vars = binary_vars)
tab %>% filter((proba_exp <= 0.1 | proba_exp >= 0.9) & sample_prop >= 0.01) %>% print(n=64)
# no viol among binary vars -> all too small sample prop


## PoRT: continuous vars uncategorised ----
source('data/port_utils.R')
a_values <- c(0.01, 0.025, 0.05, 0.1)
b_values <- c(0.01, 5/(sqrt(nrow(dat20))*log(nrow(dat20))), 0.05, 0.1)
g_values <- 1:5
lst20 <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst20[[paste0("gamma = ", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A",cov.quanti = c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10"),
             cov.quali = c("L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19", "L20"),
             data = dat20, alpha = a, beta = b, gamma = g)
    }
  }
}
lst20
#sink("output_port/port_20_leftgap_uncat.txt")
# g = 1-5: always det when should (a<=0.05 & b=0.1), for smaller b only as intersection with other cont vars
# g = 5, a = 0.1, b = 0.1: largest viol subgroup with 12.7% sample prop which was not expected


## PoRT: continuous vars categorised ----

# straightforward (preciser) categorisation
dat20_cat <- dat20
for (i in names(dat20)[1:10]) {
  dat20_cat[[i]] <- cut(dat20[[i]], breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14))
}
g_values <- 1:5
lst20_cat <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst20_cat[[paste0("gamma = ", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A", cov.quanti = NULL, cov.quali = c("L1", "L2", "L3", "L4", "L5",
                                                       "L6", "L7", "L8", "L9", "L10",
                                                       "L11", "L12", "L13", "L14", "L15",
                                                       "L16", "L17", "L18", "L19", "L20"),
             data = dat20_cat, alpha = a, beta = b, gamma = g)
    }
  }
}
lst20_cat
#sink("output_port/port_20_leftgap_cat.txt")
# g = 1-5: only ever as "L5=(1,2],(2,3],(3,4],(9,10]" which enables det for b>=0.05, but not precise enough


## kbsd ----
source("kbsd.R")
set.seed(23092025)
o5 <- dat20
o5_1 <- o5
o5_1$A <- 1
o5_2 <- o5
o5_2$A <- 0
res5 <- kbsd(data = o5,
             int_data_list = list(o5_1, o5_2),
             disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3 = sd(o5$L3), L4 = sd(o5$L4), L5 = sd(o5$L5),
                             L6 = sd(o5$L6), L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9), L10 = sd(o5$L10),
                             L11=sd(o5$L11), L12 = sd(o5$L12), L13 = sd(o5$L13), L14 = sd(o5$L14), L15 = sd(o5$L15),
                             L16 = sd(o5$L16), L17 = sd(o5$L17), L18 = sd(o5$L18), L19 = sd(o5$L19), L20 = sd(o5$L20),
                             A=0.5*sd(o5$A)), plot.out = F)
res5_plot <- kbsd(data = o5,
                  int_data_list = list(o5_1, o5_2),
                  disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3 = sd(o5$L3), L4 = sd(o5$L4), L5 = sd(o5$L5),
                                  L6 = sd(o5$L6), L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9), L10 = sd(o5$L10),
                                  L11=sd(o5$L11), L12 = sd(o5$L12), L13 = sd(o5$L13), L14 = sd(o5$L14), L15 = sd(o5$L15),
                                  L16 = sd(o5$L16), L17 = sd(o5$L17), L18 = sd(o5$L18), L19 = sd(o5$L19), L20 = sd(o5$L20),
                                  A=0.5*sd(o5$A)))  # use 1 SD for L_i, 0.5 SD for A
res5_plot  # IV=1 has slightly more support, but EDP range in [0,1] diff to interpret due to so many dims
table(dat20$A)

# what strata are those with low support for treated (IV = 1)
shift1 <- res5[res5$shift == 1,]
outliers1 <- shift1$diagnostic < quantile(shift1$diagnostic, probs = .25)  # create indices for the "outliers"
for (i in names(o5)[-21]) {
  l_values1 <- dat20[outliers1, i]
  diag_values1 <- shift1[outliers1,]
  plot(l_values1, diag_values1$diagnostic, xlab = paste0("Values of Confounder ", i), ylab = "EDP")
}  
# no clear pattern except for vars L1, L5, L7 -> det that lowest EDP are for L5<4 as expected! 


# what strata are those with low support for untreated (IV = 0)
shift2 <- res5[res5$shift == 2,]
outliers2 <- shift2$diagnostic < quantile(shift2$diagnostic, probs = .25)
for (i in names(o5)[-21]) {
  l_values2 <- dat20[outliers2, i]
  diag_values2 <- shift2[outliers2,]
  plot(l_values2, diag_values2$diagnostic, xlab = paste0("Values of Confounder ", i), ylab = "EDP")
}
# also low, but more imp only FEW obs with few EDP for L5<4 bc few obs of L5<4 in this subset of having low support (most with L5<4 have high EDP values bc a lot with L5<4 & A=0)

# essence: after cat, PoRT aggregates viol subgroup with other cat, i.e. not precise anymore


## KBSD cat ----
table(dat20_cat$L1)
table(dat20_cat$L9)  # all have same categorisation -> give numerical repr to compute kbsd values
dat20_cat <- dat20_cat %>%
  mutate(across(all_of(c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10")), 
                ~ case_when(
                  . == "(-3,-2]" ~ 1, . == "(-2,-1]"    ~ 2, . == "(-1,0]"    ~ 3,
                  . == "(0,1]"     ~ 4, . == "(1,2]"     ~ 5, . == "(2,3]"     ~ 6,
                  . == "(3,4]"     ~ 7, . == "(4,5]"     ~ 8, . == "(5,6]"     ~ 9,
                  . == "(6,7]"     ~ 10, . == "(7,8]"    ~ 11, . == "(8,9]"  ~ 12,
                  . == "(9,10]"     ~ 13, . == "(10,11]"    ~ 13, . == "(11,12]"  ~ 15,
                  . == "(12,13]"    ~ 16, . == "(13,14]"  ~ 17)))
source("kbsd.R")
o5 <- dat20_cat
o5_1 <- o5
o5_1$A <- 1
o5_2 <- o5
o5_2$A <- 0
res5 <- kbsd(data = o5,
             int_data_list = list(o5_1, o5_2),
             disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3 = sd(o5$L3), L4 = sd(o5$L4), L5 = sd(o5$L5),
                            L6 = sd(o5$L6), L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9), L10 = sd(o5$L10),
                            L11=sd(o5$L11), L12 = sd(o5$L12), L13 = sd(o5$L13), L14 = sd(o5$L14), L15 = sd(o5$L15),
                            L16 = sd(o5$L16), L17 = sd(o5$L17), L18 = sd(o5$L18), L19 = sd(o5$L19), L20 = sd(o5$L20),
                            A=0.5*sd(o5$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
res5_plot <- kbsd(data = o5,
                  int_data_list = list(o5_1, o5_2),
                  disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3 = sd(o5$L3), L4 = sd(o5$L4), L5 = sd(o5$L5),
                                 L6 = sd(o5$L6), L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9), L10 = sd(o5$L10),
                                 L11=sd(o5$L11), L12 = sd(o5$L12), L13 = sd(o5$L13), L14 = sd(o5$L14), L15 = sd(o5$L15),
                                 L16 = sd(o5$L16), L17 = sd(o5$L17), L18 = sd(o5$L18), L19 = sd(o5$L19), L20 = sd(o5$L20),
                                 A=0.5*sd(o5$A)))
res5_plot  # IV=1 with a bit more support
table(o5$A)  # also more obs for that

# what strata are those with low support for treated (IV = 1): should be L3<4
subset_5_1 <- res5[res5$diagnostic < median(res5[res5$shift == 1, "diagnostic"]) & res5$shift == 1,]
table(o5[subset_5_1$observation, c("L3")]) #  as expected; most from L3=[2,4] that have low support in IV=1

# what strata are those with low support for treated (IV = 0)
subset_5_2 <- res5[res5$diagnostic < median(res5[res5$shift == 2, "diagnostic"]) & res5$shift == 2,]
table(o5[subset_5_2$observation, ][, c("L3")]) # no viol expected

# essence: kbsd also det the critical stratum L3<4 after cat



# 20 Confounders With Middle Gap ----
set.seed(23092025)
L5_1 <- rnorm(500, 3, 1)
L5_2 <- rnorm(500, 7, 1)
L5_1_2 <- c(L5_1, L5_2)
plot(density(L5_1_2))
DAG20 <- DAG.empty() +
  node("L1", distr = "rnorm", mean = 1, sd = 1) +
  node("L2", distr = "rnorm", mean = 2, sd = 1) +
  node("L3", distr = "rnorm", mean = 3, sd = 1) +
  node("L4", distr = "rnorm", mean = 4, sd = 1) +
  node("L5", distr = "rconst", const = L5_1_2) +
  node("L6", distr = "rnorm", mean = 6, sd = 1) +
  node("L7", distr = "rnorm", mean = 7, sd = 1) +
  node("L8", distr = "rnorm", mean = 8, sd = 1) +
  node("L9", distr = "rnorm", mean = 9, sd = 1) +
  node("L10", distr = "rnorm", mean = 10, sd = 1) +
  node("L11", distr = "rbern", prob = 0.5) +
  node("L12", distr = "rbern", prob = 0.5) +
  node("L13", distr = "rbern", prob = 0.5) +
  node("L14", distr = "rbern", prob = 0.5) +
  node("L15", distr = "rbern", prob = 0.5) +
  node("L16", distr = "rbern", prob = 0.5) +
  node("L17", distr = "rbern", prob = 0.5) +
  node("L18", distr = "rbern", prob = 0.5) +
  node("L19", distr = "rbern", prob = 0.5) +
  node("L20", distr = "rbern", prob = 0.5) +
  node("A", distr = "rbern", prob = plogis(L1 - L20 - 2*L5*(L5 > 4 & L5 < 6)))  # now there should be viol if L5 is v small
DAG20 <- set.DAG(DAG20)
dat20 <- sim(DAG20, rndseed = 12082025, n = 1000)[-1]
table(dat20$A)  # balanced

dat20 %>% filter(L5 > 4 & L5 < 6 & A==1) %>% nrow()/
  dat20 %>% filter(L5 > 4 & L5 < 6) %>% nrow()  # P(A=1)~0 -> to be det for g>=1, any b & any a (sample prop= 17%)

# table to check for all combos of binary conf if there are any with extreme P(A): fun defined in setup.R
binary_vars <- paste0("L", c(6:10))
tab <- make_strata_table(dat20, A = "A", binary_vars = binary_vars)
tab %>% filter((proba_exp <= 0.1 | proba_exp >= 0.9) & sample_prop >= 0.01) %>% print(n=64)
# no viol among binary vars -> all too small sample prop


## PoRT: continuous vars uncategorised ----
source('data/port_utils.R')
a_values <- c(0.01, 0.025, 0.05, 0.1)
b_values <- c(0.01, 5/(sqrt(nrow(dat20))*log(nrow(dat20))), 0.05, 0.1)
g_values <- 1:5
lst20 <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst20[[paste0("gamma = ", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A",cov.quanti = c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10"),
             cov.quali = c("L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19", "L20"),
             data = dat20, alpha = a, beta = b, gamma = g)
    }
  }
}
lst20
#sink("output_port/port_20_bimodal_uncat.txt")
# g = 1-5: viol det for all a & b, most of other viol only involve cont conf bc
#          due to greedy cat of cont conf (but also usually small sample size)


## PoRT: continuous vars categorised ----

# straightforward (preciser) categorisation
dat20_cat <- dat20
for (i in names(dat20)[1:10]) {
  dat20_cat[[i]] <- cut(dat20[[i]], breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14))
}
g_values <- 1:5
lst20_cat <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst20_cat[[paste0("gamma = ", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A", cov.quanti = NULL, cov.quali = c("L1", "L2", "L3", "L4", "L5",
                                                       "L6", "L7", "L8", "L9", "L10",
                                                       "L11", "L12", "L13", "L14", "L15",
                                                       "L16", "L17", "L18", "L19", "L20"),
             data = dat20_cat, alpha = a, beta = b, gamma = g)
    }
  }
}
lst20_cat
#sink("output_port/port_20_bimodal_cat.txt")
# g = 1-5: also always det but only as "L5=(-1,0],(4,5],(5,6]" -> tried to make stratum as large as poss?


## kbsd ----
source("kbsd.R")
set.seed(23092025)
o5 <- dat20
o5_1 <- o5
o5_1$A <- 1
o5_2 <- o5
o5_2$A <- 0
res5 <- kbsd(data = o5,
             int_data_list = list(o5_1, o5_2),
             disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3 = sd(o5$L3), L4 = sd(o5$L4), L5 = sd(o5$L5),
                            L6 = sd(o5$L6), L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9), L10 = sd(o5$L10),
                            L11=sd(o5$L11), L12 = sd(o5$L12), L13 = sd(o5$L13), L14 = sd(o5$L14), L15 = sd(o5$L15),
                            L16 = sd(o5$L16), L17 = sd(o5$L17), L18 = sd(o5$L18), L19 = sd(o5$L19), L20 = sd(o5$L20),
                            A=0.5*sd(o5$A)), plot.out = F)
res5_plot <- kbsd(data = o5,
                  int_data_list = list(o5_1, o5_2),
                  disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3 = sd(o5$L3), L4 = sd(o5$L4), L5 = sd(o5$L5),
                                 L6 = sd(o5$L6), L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9), L10 = sd(o5$L10),
                                 L11=sd(o5$L11), L12 = sd(o5$L12), L13 = sd(o5$L13), L14 = sd(o5$L14), L15 = sd(o5$L15),
                                 L16 = sd(o5$L16), L17 = sd(o5$L17), L18 = sd(o5$L18), L19 = sd(o5$L19), L20 = sd(o5$L20),
                                 A=0.5*sd(o5$A)))  # use 1 SD for L_i, 0.5 SD for A
res5_plot  # IV=2 (A=0) has slightly more support, but again EDP range in [0,1] diff to interpret due to so many dims
table(dat20$A)  # also more obs for A=0

# what strata are those with low support for treated (IV = 1)
shift1 <- res5[res5$shift == 1,]
outliers1 <- shift1$diagnostic < quantile(shift1$diagnostic, probs = .25)  # create indices for the "outliers"
for (i in names(o5)[-21]) {
  l_values1 <- dat20[outliers1, i]
  diag_values1 <- shift1[outliers1,]
  plot(l_values1, diag_values1$diagnostic, xlab = paste0("Values of Confounder ", i), ylab = "EDP")
}  
# no clear pattern for L5 in [4,6] being critical here


# what strata are those with low support for untreated (IV = 0)
shift2 <- res5[res5$shift == 2,]
outliers2 <- shift2$diagnostic < quantile(shift2$diagnostic, probs = .25)
for (i in names(o5)[-21]) {
  l_values2 <- dat20[outliers2, i]
  diag_values2 <- shift2[outliers2,]
  plot(l_values2, diag_values2$diagnostic, xlab = paste0("Values of Confounder ", i), ylab = "EDP")
}
# no obs with few EDP for L5<4 bc all obs of L5<4 have high EDP values as all have L5<4 & A=0

# essence: for cat case, PoRT is not precise anymore (agg viol stratum with other category, 
# maybe to make stratum as large as poss?), also kbsd did not detect that few support in L5=[4,6]


## KBSD cat ----
table(dat20_cat$L1)
table(dat20_cat$L2)  # all have same categorisation -> give numerical repr to compute kbsd values
dat20_cat <- dat20_cat %>%
  mutate(across(all_of(c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10")), 
                ~ case_when(
                  . == "(-3,-2]" ~ 1, . == "(-2,-1]"    ~ 2, . == "(-1,0]"    ~ 3,
                  . == "(0,1]"     ~ 4, . == "(1,2]"     ~ 5, . == "(2,3]"     ~ 6,
                  . == "(3,4]"     ~ 7, . == "(4,5]"     ~ 8, . == "(5,6]"     ~ 9,
                  . == "(6,7]"     ~ 10, . == "(7,8]"    ~ 11, . == "(8,9]"  ~ 12,
                  . == "(9,10]"     ~ 13, . == "(10,11]"    ~ 13, . == "(11,12]"  ~ 15,
                  . == "(12,13]"    ~ 16, . == "(13,14]"  ~ 17)))
summary(dat20_cat)
source("kbsd.R")
o5 <- dat20_cat
o5_1 <- o5
o5_1$A <- 1
o5_2 <- o5
o5_2$A <- 0
res5 <- kbsd(data = o5,
             int_data_list = list(o5_1, o5_2),
             disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3 = sd(o5$L3), L4 = sd(o5$L4), L5 = sd(o5$L5),
                            L6 = sd(o5$L6), L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9), L10 = sd(o5$L10),
                            L11=sd(o5$L11), L12 = sd(o5$L12), L13 = sd(o5$L13), L14 = sd(o5$L14), L15 = sd(o5$L15),
                            L16 = sd(o5$L16), L17 = sd(o5$L17), L18 = sd(o5$L18), L19 = sd(o5$L19), L20 = sd(o5$L20),
                            A=0.5*sd(o5$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
res5_plot <- kbsd(data = o5,
                  int_data_list = list(o5_1, o5_2),
                  disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3 = sd(o5$L3), L4 = sd(o5$L4), L5 = sd(o5$L5),
                                 L6 = sd(o5$L6), L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9), L10 = sd(o5$L10),
                                 L11=sd(o5$L11), L12 = sd(o5$L12), L13 = sd(o5$L13), L14 = sd(o5$L14), L15 = sd(o5$L15),
                                 L16 = sd(o5$L16), L17 = sd(o5$L17), L18 = sd(o5$L18), L19 = sd(o5$L19), L20 = sd(o5$L20),
                                 A=0.5*sd(o5$A)))
res5_plot  # overall v few EDP as high-dim adjustment set
table(o5$A)  # a few less obs for A=1 matches that IV=1 has less support

# what strata are those with low support for treated (IV = 1)
subset_5_1 <- res5[res5$diagnostic < median(res5[res5$shift == 1, "diagnostic"]) & res5$shift == 1,]
table(o5[subset_5_1$observation, c("L3")]) # most from L3=[3,4] that have low support in IV=1 -> should be [4,6] tho

# what strata are those with low support for treated (IV = 0): should be L19=1 & L20=1
subset_5_2 <- res5[res5$diagnostic < median(res5[res5$shift == 2, "diagnostic"]) & res5$shift == 2,]
table(o5[subset_5_2$observation, ][, c("L3")]) # no viol expected

# essence: after cat, kbsd did not find stratum with viol for A=1




# 20 Confounders Correlated ----

set.seed(23092025)
DAG20 <- DAG.empty() +
  node("L1", distr = "rnorm", mean = 1, sd = 1) +
  node("L2", distr = "rnorm", mean = 2, sd = 1) +
  node("L3", distr = "rnorm", mean = 3, sd = 1) +
  node("L4", distr = "rnorm", mean = 4, sd = 1) +
  node("L5", distr = "rnorm", mean = 5, sd = 1) +
  node("L6", distr = "rnorm", mean = 6, sd = 1) +
  node("L7", distr = "rnorm", mean = 7, sd = 1) +
  node("L8", distr = "rnorm", mean = 8, sd = 1) +
  node("L9", distr = "rnorm", mean = 9, sd = 1) +
  node("L10", distr = "rnorm", mean = 10, sd = 1) +
  node("L11", distr = "rbern", prob = 0.5) +
  node("L12", distr = "rbern", prob = 0.5) +
  node("L13", distr = "rbern", prob = 0.5) +
  node("L14", distr = "rbern", prob = 0.5) +
  node("L15", distr = "rbern", prob = 0.5) +
  node("L16", distr = "rbern", prob = 0.5) +
  node("L17", distr = "rbern", prob = 0.5) +
  node("L18", distr = "rbern", prob = 0.5) +
  node("L19", distr = "rbern", prob = 0.5) +
  node("L20", distr = "rbern", prob = 0.5) +
  node("A", distr = "rbern", prob = plogis(L1 - L2 + L5*L19*L20))
# if one of L19 and L20 = 0, then P(A)~0.26, if both =1 then P(A)~1, i.e. viol for P(A=0|L)~0
DAG20 <- set.DAG(DAG20)
dat20 <- sim(DAG20, rndseed = 12082025, n = 1000)[-1]
table(dat20$A)  # balanced

# induce correlations
sum(cor(dat20) > 0.3)
dat20$L2 <- 2*dat20$L1 + rnorm(1000,0,1)
cor(dat20$L2, dat20$L1)  # cor of 0.9
dat20$L18 <- dat20$L17
cor(dat20$L18, dat20$L17)  # cor of 1

dat20 %>% filter(L19==1 & L20 == 1 & A==1) %>% nrow()/
  dat20 %>% filter(L19==1 & L20 == 1) %>% nrow()  # P(A=1)~1 -> to be det for g>=2, b=0.1, any a (sample prop= 22.7%)

# table to check for all combos of binary conf if there are any with extreme P(A): fun defined in setup.R
binary_vars <- paste0("L", c(10:20))
tab <- make_strata_table(dat20, A = "A", binary_vars = binary_vars)
tab %>% filter((proba_exp <= 0.1 | proba_exp >= 0.9) & sample_prop >= 0.01) %>% arrange(desc(sample_prop))
# almost all involve L19 & L20 as wanted, only 9 that don't but with sample prop =1.1%/1.2% only


## PoRT: continuous vars uncategorised ----
source('data/port_utils.R')
gruber <- 5/(sqrt(nrow(dat2))*log(nrow(dat2)))
port("A", cov.quanti = c("L1","L2","L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10"),
     cov.quali = c("L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19", "L20"),
     data = dat20, alpha = 0.05, beta = 0.1, gamma = 3)
# no longer work with a=0.02/0.03/0.04!
# but g=3: undet for a=0.01/0.025/0.05, det for a= 0.1
# (g= 2,3: undet for all b = 0.1 & a=0.01/0.02/gruber/0.05 (too small alpha lets algo focus on smaller strata?), det: a=0.1, 0.04, 0.03)

a_values <- c(0.01, 0.025, 0.05, 0.1)
b_values <- c(0.01, 5/(sqrt(nrow(dat20))*log(nrow(dat20))), 0.05, 0.1)
g_values <- 1:5
lst20 <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst20[[paste0("gamma = ", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A",cov.quanti = c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10"),
             cov.quali = c("L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19", "L20"),
             data = dat20, alpha = a, beta = b, gamma = g)
    }
  }
}
lst20
#sink("port_20_corr_uncat.txt")
# g = 1: only viol with cont conf (greedy cat), viol not yet poss bc intersection of 2
# g = 2-5: only det from a = 0.1 (& b=0.1) although should be for any a with sample prop = 22.7%
#          i.e. suspected focus on viol of smaller strata for smaller a? 



## PoRT: continuous vars categorised ----

dat20_cat <- dat20
for (i in names(dat20)[1:10]) {
  dat20_cat[[i]] <- cut(dat20[[i]], breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14))
}
dat20_cat
g_values <- 1:5
lst20_cat <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst20_cat[[paste0("gamma = ", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A", cov.quanti = NULL, cov.quali = c("L1", "L2", "L3", "L4", "L5",
                                                       "L6", "L7", "L8", "L9", "L10",
                                                       "L11", "L12", "L13", "L14", "L15",
                                                       "L16", "L17", "L18", "L19", "L20"),
             data = dat20_cat, alpha = a, beta = b, gamma = g)
    }
  }
}
lst20_cat
#sink("output_port/port_20_corr_cat.txt")
# g = 1: our viol not poss to find yet
# g=2: now alr det from a = 0.05 (& b=0.1), so still not for all a as should but better than uncategorised


## kbsd ----
source("kbsd.R")
o6 <- dat20
o6_1 <- o6
o6_1$A <- 1
o6_2 <- o6
o6_2$A <- 0
res6 <- kbsd(data = o6,
             int_data_list = list(o6_1, o6_2),
             disthalf_vec=c(L1=sd(o6$L1), L2 = sd(o6$L2), L3 = sd(o6$L3), L4 = sd(o6$L4), L5 = sd(o6$L5),
                            L6 = sd(o6$L6), L7 = sd(o6$L7), L8 = sd(o6$L8), L9 = sd(o6$L9), L10 = sd(o6$L10),
                            L11=sd(o6$L11), L12 = sd(o6$L12), L13 = sd(o6$L13), L14 = sd(o6$L14), L15 = sd(o6$L15),
                            L16 = sd(o6$L16), L17 = sd(o6$L17), L18 = sd(o6$L18), L19 = sd(o6$L19), L20 = sd(o6$L20),
                            A=0.5*sd(o6$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
res6_plot <- kbsd(data = o6,
             int_data_list = list(o6_1, o6_2),
             disthalf_vec=c(L1=sd(o6$L1), L2 = sd(o6$L2), L3 = sd(o6$L3), L4 = sd(o6$L4), L5 = sd(o6$L5),
                            L6 = sd(o6$L6), L7 = sd(o6$L7), L8 = sd(o6$L8), L9 = sd(o6$L9), L10 = sd(o6$L10),
                            L11=sd(o6$L11), L12 = sd(o6$L12), L13 = sd(o6$L13), L14 = sd(o6$L14), L15 = sd(o6$L15),
                            L16 = sd(o6$L16), L17 = sd(o6$L17), L18 = sd(o6$L18), L19 = sd(o6$L19), L20 = sd(o6$L20),
                            A=0.5*sd(o6$A)))
res6_plot  # slightly less support for A=1
table(dat20$A)

# what strata are those with low support for treated (IV = 1)
subset_6_1 <- res6[res6$diagnostic <= quantile(res6[res6$shift == 1, "diagnostic"], 0.25) & res6$shift == 1,]  # many obs with v low EDP
table(o6[subset_6_1$observation, c("L19", "L20")])
# those with few support in IV=1 are L19=0 & L20=0

# what strata are those with low support for treated (IV = 0): should be L19=1 & L20=1
subset_6_2 <- res6[res6$diagnostic < quantile(res6[res6$shift == 2, "diagnostic"], 0.25) & res6$shift == 2,]
table(o6[subset_6_2$observation, c("L19", "L20")])
# highest count among IV=0 is L19=1 & L20=1 -> viol found

# essence: PoRT & kbsd found viol, but PoRT did not always when should (cat made it slightly better, but still undet for smaller a, however smaller a also not that imp maybe)


## KBSD cat ----
table(dat20_cat$L1)
table(dat20_cat$L2)  # all have same categorisation -> give numerical repr to compute kbsd values
dat20_cat <- dat20_cat %>%
  mutate(across(all_of(c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10")), 
                ~ case_when(
                  . == "(-3,-2]" ~ 1, . == "(-2,-1]"    ~ 2, . == "(-1,0]"    ~ 3,
                  . == "(0,1]"     ~ 4, . == "(1,2]"     ~ 5, . == "(2,3]"     ~ 6,
                  . == "(3,4]"     ~ 7, . == "(4,5]"     ~ 8, . == "(5,6]"     ~ 9,
                  . == "(6,7]"     ~ 10, . == "(7,8]"    ~ 11, . == "(8,9]"  ~ 12,
                  . == "(9,10]"     ~ 13, . == "(10,11]"    ~ 13, . == "(11,12]"  ~ 15,
                  . == "(12,13]"    ~ 16, . == "(13,14]"  ~ 17)))
source("kbsd.R")
o5 <- dat20_cat
o5_1 <- o5
o5_1$A <- 1
o5_2 <- o5
o5_2$A <- 0
res5 <- kbsd(data = o5,
             int_data_list = list(o5_1, o5_2),
             disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3 = sd(o5$L3), L4 = sd(o5$L4), L5 = sd(o5$L5),
                            L6 = sd(o5$L6), L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9), L10 = sd(o5$L10),
                            L11=sd(o5$L11), L12 = sd(o5$L12), L13 = sd(o5$L13), L14 = sd(o5$L14), L15 = sd(o5$L15),
                            L16 = sd(o5$L16), L17 = sd(o5$L17), L18 = sd(o5$L18), L19 = sd(o5$L19), L20 = sd(o5$L20),
                            A=0.5*sd(o5$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
res5_plot <- kbsd(data = o5,
                  int_data_list = list(o5_1, o5_2),
                  disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3 = sd(o5$L3), L4 = sd(o5$L4), L5 = sd(o5$L5),
                                 L6 = sd(o5$L6), L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9), L10 = sd(o5$L10),
                                 L11=sd(o5$L11), L12 = sd(o5$L12), L13 = sd(o5$L13), L14 = sd(o5$L14), L15 = sd(o5$L15),
                                 L16 = sd(o5$L16), L17 = sd(o5$L17), L18 = sd(o5$L18), L19 = sd(o5$L19), L20 = sd(o5$L20),
                                 A=0.5*sd(o5$A)))
res5_plot  # overall v few EDP as high-dim adjustment set
table(o5$A)  # a few less obs for A=1 matches that IV=1 has less support

# what strata are those with low support for treated (IV = 1)
subset_5_1 <- res5[res5$diagnostic < median(res5[res5$shift == 1, "diagnostic"]) & res5$shift == 1,]
table(o5[subset_5_1$observation, c("L19", "L20")]) # no viol expected

# what strata are those with low support for treated (IV = 0): should be L19=1 & L20=1, bc for these values will always get A=1, rarely A=0
subset_5_2 <- res5[res5$diagnostic < median(res5[res5$shift == 2, "diagnostic"]) & res5$shift == 2,]
table(o5[subset_5_2$observation, ][, c("L19", "L20")]) # highest count among IV=0 is L19=1 & L20=1 -> as expected

# essence: kbsd det viol also after cat
