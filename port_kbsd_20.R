source("setup.R")

# 20 Confounders ----
set.seed(23092025)
DAG1 <- DAG.empty() +
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
DAG1 <- set.DAG(DAG1)
dat1 <- sim(DAG1, rndseed = 12082025, n = 1000)[-1]
table(dat1$A)  # good that not too imbalanced
dat1 %>% filter(L19==1 & L20 == 1 & A==1) %>% nrow()/
  dat1 %>% filter(L19==1 & L20 == 1) %>% nrow()  # P(A=1)~1 -> to be det for g>=2, b=0.1, any a (sample prop= 22.7%)

# table to check for all combos of binary conf if there are any with extreme P(A): fun defined in setup.R
binary_vars <- paste0("L", c(1:10))
tab <- make_strata_table(dat1, A = "A", binary_vars = binary_vars)
tab %>% filter((proba_exp <= 0.1 | proba_exp >= 0.9) & sample_prop >= 0.01) %>% print(n=64)
# no viol among binary vars -> all too small sample prop


## PoRT: continuous vars uncategorised ----
source('data/port_utils.R')
gruber <- 5/(sqrt(nrow(dat2))*log(nrow(dat2)))
port("A", cov.quanti = c("L1","L2","L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10"),
     cov.quali = c("L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19", "L20"),
     data = dat1, alpha = 0.01, beta = 0.1, gamma = 3)
# no longer work with a=0.02/0.03/0.04!
# but g=3: undet for a=0.01/0.025/0.05, det for a= 0.1
# (g= 2,3: undet for all b = 0.1 & a=0.01/0.02/gruber/0.05 (too small alpha lets algo focus on smaller strata?), det: a=0.1, 0.04, 0.03)

a_values <- c(0.01, 0.025, 0.05, 0.1)
b_values <- c(0.01, 5/(sqrt(nrow(dat1))*log(nrow(dat1))), 0.05, 0.1)
g_values <- 1:20
lst20 <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst20[[paste0("gamma = ", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A",cov.quanti = c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10"),
             cov.quali = c("L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19", "L20"),
             data = dat1, alpha = a, beta = b, gamma = g)
    }
  }
}
lst20
#sink("port_20_uncat.txt")
# g = 1: not poss yet bc intersection of 2
# g = 2-5: only det for a=0.1 & b=0.1, for no smaller a & b=0.1 even though should!


## PoRT: continuous vars categorised ----

# straightforward (preciser) categorisation
dat1_cat <- dat1
for (i in names(dat1)[1:10]) {
  dat1_cat[[i]] <- cut(dat1[[i]], breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14))
}
dat1_cat
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
             data = dat1_cat, alpha = a, beta = b, gamma = g)
    }
  }
}
lst20_cat
#sink("output_port/port_20_cat.txt")
# g = 1: only other viol with cont vars & only with v small a=0.01/0.025
# g = 2-5: det from a>=0.025 & b=0.1! alr shows pos eff of cat CONT conf on detection of viol involving CAT conf


## kbsd ----
source("kbsd.R")
o1 <- dat1
o1_1 <- o1
o1_1$A <- 1
o1_2 <- o1
o1_2$A <- 0
res1 <- kbsd(data = o1,
             int_data_list = list(o1_1, o1_2),
             disthalf_vec=c(L1=sd(o1$L1), L2 = sd(o1$L2), L3 = sd(o1$L3), L4 = sd(o1$L4), L5 = sd(o1$L5),
                            L6 = sd(o1$L6), L7 = sd(o1$L7), L8 = sd(o1$L8), L9 = sd(o1$L9), L10 = sd(o1$L10),
                            L11=sd(o1$L11), L12 = sd(o1$L12), L13 = sd(o1$L13), L14 = sd(o1$L14), L15 = sd(o1$L15),
                            L16 = sd(o1$L16), L17 = sd(o1$L17), L18 = sd(o1$L18), L19 = sd(o1$L19), L20 = sd(o1$L20),
                            A=0.5*sd(o1$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
res1_plot <- kbsd(data = o1,
             int_data_list = list(o1_1, o1_2),
             disthalf_vec=c(L1=sd(o1$L1), L2 = sd(o1$L2), L3 = sd(o1$L3), L4 = sd(o1$L4), L5 = sd(o1$L5),
                            L6 = sd(o1$L6), L7 = sd(o1$L7), L8 = sd(o1$L8), L9 = sd(o1$L9), L10 = sd(o1$L10),
                            L11=sd(o1$L11), L12 = sd(o1$L12), L13 = sd(o1$L13), L14 = sd(o1$L14), L15 = sd(o1$L15),
                            L16 = sd(o1$L16), L17 = sd(o1$L17), L18 = sd(o1$L18), L19 = sd(o1$L19), L20 = sd(o1$L20),
                            A=0.5*sd(o1$A)))
res1_plot

# overall v low EDP which makes sense as all obs are v spread out over 20 dims
# -> many defined strata so that within each prob just one treatment level -> no support for the other
#    was warned for in kbsd paper
# -> but interesting that among IV=1 (bringing all to A=1), particularly low support (many with EDP = 0)

# what strata are those with low support for treated (IV = 1)
subset_1_1 <- res1[res1$diagnostic <= median(res1[res1$shift == 1, "diagnostic"]) & res1$shift == 1,]  # many obs with v low EDP
table(o1[subset_1_1$observation, c("L19", "L20")])
# those with few support in IV=1 are L19=0 & L20=0

# what strata are those with low support for treated (IV = 0): should be L8=1 & L6=1 & L7=1
subset_1_2 <- res1[res1$diagnostic < median(res1[res1$shift == 2, "diagnostic"]) & res1$shift == 2,]
table(o1[subset_1_2$observation, c("L19", "L20")]) # highest count among IV=0 is L19=1 & L20=1 -> viol found

# essence: imp insight that cat of cont vars impacts detection of viol in CATEGORICAL vars!
#          -> here, for broader cat higher chance of det with all alpha > 0.01, without cat only for alpha <= 0.03
#          kbsd also flagged those from L19=1 & L20=1 as with few support for A=0 so overall det by both diags


# check if viol stratum found:
outliers_ind0 <- unique(res1[res1$shift == 2 & res1$diagnostic < quantile(res1[res1$shift == 2, "diagnostic"], 0.05), "observation"])
o1[outliers_ind0,] %>% nrow()
o1[outliers_ind0,] %>% filter(L19==1 & L20==1) %>% nrow()  # majority of outliers are from critical stratum


### alternative EDP calc for high-dim covar set ----

# type = "minval" instead of default type = "Rfast"
# requires specification of EDP minvalue accepted for each covar (more imp <=> lower/no minval, less imp <=> 0.5)
res1 <- kbsd(data = o1, int_data_list = list(o1_1, o1_2), type = "minval",
             disthalf_vec=c(L1=sd(o1$L1), L2 = sd(o1$L2), L3 = sd(o1$L3), L4 = sd(o1$L4), L5 = sd(o1$L5),
                            L6 = sd(o1$L6), L7 = sd(o1$L7), L8 = sd(o1$L8), L9 = sd(o1$L9), L10 = sd(o1$L10),
                            L11=sd(o1$L11), L12 = sd(o1$L12), L13 = sd(o1$L13), L14 = sd(o1$L14), L15 = sd(o1$L15),
                            L16 = sd(o1$L16), L17 = sd(o1$L17), L18 = sd(o1$L18), L19 = sd(o1$L19), L20 = sd(o1$L20),
                            A=0.5*sd(o1$A)),  # use 1 SD for L_i, 0.5 SD for A
             minval_vec = , plot.out = T)

# type = "harmonicmean" instead of default type = "Rfast"
res1_plot <- kbsd(data = o1, int_data_list = list(o1_1, o1_2), type = "harmonicmean",
             disthalf_vec=c(L1=sd(o1$L1), L2 = sd(o1$L2), L3 = sd(o1$L3), L4 = sd(o1$L4), L5 = sd(o1$L5),
                            L6 = sd(o1$L6), L7 = sd(o1$L7), L8 = sd(o1$L8), L9 = sd(o1$L9), L10 = sd(o1$L10),
                            L11=sd(o1$L11), L12 = sd(o1$L12), L13 = sd(o1$L13), L14 = sd(o1$L14), L15 = sd(o1$L15),
                            L16 = sd(o1$L16), L17 = sd(o1$L17), L18 = sd(o1$L18), L19 = sd(o1$L19), L20 = sd(o1$L20),
                            A=0.5*sd(o1$A)))
res1 <- kbsd(data = o1, int_data_list = list(o1_1, o1_2), type = "harmonicmean",
             disthalf_vec=c(L1=sd(o1$L1), L2 = sd(o1$L2), L3 = sd(o1$L3), L4 = sd(o1$L4), L5 = sd(o1$L5),
                            L6 = sd(o1$L6), L7 = sd(o1$L7), L8 = sd(o1$L8), L9 = sd(o1$L9), L10 = sd(o1$L10),
                            L11=sd(o1$L11), L12 = sd(o1$L12), L13 = sd(o1$L13), L14 = sd(o1$L14), L15 = sd(o1$L15),
                            L16 = sd(o1$L16), L17 = sd(o1$L17), L18 = sd(o1$L18), L19 = sd(o1$L19), L20 = sd(o1$L20),
                            A=0.5*sd(o1$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
#saveRDS(res1, "hm_20_uncat_uncorr.RDS")
# checking if viol stratum L19=1 & L20=1 included in the outliers
outliers_ind0 <- unique(res1[res1$shift == 2 & res1$diagnostic < quantile(res1[res1$shift == 2, "diagnostic"], 0.05), "observation"])
o1[outliers_ind0,] %>% nrow()
o1[outliers_ind0,] %>% filter(L19==1 & L20==1) %>% nrow()  # majority of outliers are from critical stratum




## KBSD cat ----
table(dat1_cat$L1)
table(dat1_cat$L2)  # all have same categorisation -> give numerical repr to compute kbsd values
dat1_cat <- dat1_cat %>%
  mutate(across(all_of(c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10")), 
                ~ case_when(
                  . == "(-3,-2]" ~ 1, . == "(-2,-1]"    ~ 2, . == "(-1,0]"    ~ 3,
                  . == "(0,1]"     ~ 4, . == "(1,2]"     ~ 5, . == "(2,3]"     ~ 6,
                  . == "(3,4]"     ~ 7, . == "(4,5]"     ~ 8, . == "(5,6]"     ~ 9,
                  . == "(6,7]"     ~ 10, . == "(7,8]"    ~ 11, . == "(8,9]"  ~ 12,
                  . == "(9,10]"     ~ 13, . == "(10,11]"    ~ 13, . == "(11,12]"  ~ 15,
                  . == "(12,13]"    ~ 16, . == "(13,14]"  ~ 17)))
source("kbsd.R")
o1 <- dat1_cat
o1_1 <- o1
o1_1$A <- 1
o1_2 <- o1
o1_2$A <- 0
res1 <- kbsd(data = o1,
             int_data_list = list(o1_1, o1_2),
             disthalf_vec=c(L1=sd(o1$L1), L2 = sd(o1$L2), L3 = sd(o1$L3), L4 = sd(o1$L4), L5 = sd(o1$L5),
                            L6 = sd(o1$L6), L7 = sd(o1$L7), L8 = sd(o1$L8), L9 = sd(o1$L9), L10 = sd(o1$L10),
                            L11=sd(o1$L11), L12 = sd(o1$L12), L13 = sd(o1$L13), L14 = sd(o1$L14), L15 = sd(o1$L15),
                            L16 = sd(o1$L16), L17 = sd(o1$L17), L18 = sd(o1$L18), L19 = sd(o1$L19), L20 = sd(o1$L20),
                            A=0.5*sd(o1$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
res1_plot <- kbsd(data = o1,
                  int_data_list = list(o1_1, o1_2),
                  disthalf_vec=c(L1=sd(o1$L1), L2 = sd(o1$L2), L3 = sd(o1$L3), L4 = sd(o1$L4), L5 = sd(o1$L5),
                                 L6 = sd(o1$L6), L7 = sd(o1$L7), L8 = sd(o1$L8), L9 = sd(o1$L9), L10 = sd(o1$L10),
                                 L11=sd(o1$L11), L12 = sd(o1$L12), L13 = sd(o1$L13), L14 = sd(o1$L14), L15 = sd(o1$L15),
                                 L16 = sd(o1$L16), L17 = sd(o1$L17), L18 = sd(o1$L18), L19 = sd(o1$L19), L20 = sd(o1$L20),
                                 A=0.5*sd(o1$A)))
res1_plot  # overall v few EDP as high-dim adjustment set
table(o1$A)  # a few less obs for A=1 matches that IV=1 has less support

# check if viol stratum found:
outliers_ind0 <- unique(res1[res1$shift == 2 & res1$diagnostic < quantile(res1[res1$shift == 2, "diagnostic"], 0.05), "observation"])
o1[outliers_ind0,] %>% nrow()
o1[outliers_ind0,] %>% filter(L19==1 & L20==1) %>% nrow()  # still majority of outliers are from critical stratum

# essence: kbsd det viol also for cat data

# hm
res1 <- kbsd(data = o1, int_data_list = list(o1_1, o1_2), type = "harmonicmean",
             disthalf_vec=c(L1=sd(o1$L1), L2 = sd(o1$L2), L3 = sd(o1$L3), L4 = sd(o1$L4), L5 = sd(o1$L5),
                            L6 = sd(o1$L6), L7 = sd(o1$L7), L8 = sd(o1$L8), L9 = sd(o1$L9), L10 = sd(o1$L10),
                            L11=sd(o1$L11), L12 = sd(o1$L12), L13 = sd(o1$L13), L14 = sd(o1$L14), L15 = sd(o1$L15),
                            L16 = sd(o1$L16), L17 = sd(o1$L17), L18 = sd(o1$L18), L19 = sd(o1$L19), L20 = sd(o1$L20),
                            A=0.5*sd(o1$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
#saveRDS(res1, "hm_20_cat_uncorr.RDS")
# checking if viol stratum L19=1 & L20=1 included in the outliers
outliers_ind0 <- unique(res1[res1$shift == 2 & res1$diagnostic < quantile(res1[res1$shift == 2, "diagnostic"], 0.05), "observation"])
o1[outliers_ind0,] %>% nrow()
o1[outliers_ind0,] %>% filter(L19==1 & L20==1) %>% nrow()  # majority of outliers are from critical stratum




# 20 Confounders Correlated ----

set.seed(23092025)
DAG2 <- DAG.empty() +
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
DAG2 <- set.DAG(DAG2)
dat2 <- sim(DAG2, rndseed = 12082025, n = 1000)[-1]
table(dat2$A)  # balanced

# induce correlations
sum(cor(dat2) > 0.3)
dat2$L2 <- 2*dat2$L1 + rnorm(1000,0,1)
cor(dat2$L2, dat2$L1)  # cor of 0.9
dat2$L18 <- dat2$L17
cor(dat2$L18, dat2$L17)  # cor of 1

dat2 %>% filter(L19==1 & L20 == 1 & A==1) %>% nrow()/
  dat2 %>% filter(L19==1 & L20 == 1) %>% nrow()  # P(A=1)~1 -> to be det for g>=2, b=0.1, any a (sample prop= 22.7%)

# table to check for all combos of binary conf if there are any with extreme P(A): fun defined in setup.R
binary_vars <- paste0("L", c(10:20))
tab <- make_strata_table(dat2, A = "A", binary_vars = binary_vars)
tab %>% filter((proba_exp <= 0.1 | proba_exp >= 0.9) & sample_prop >= 0.01) %>% arrange(desc(sample_prop))
# almost all involve L19 & L20 as wanted, only 9 that don't but with sample prop =1.1%/1.2% only


## PoRT: continuous vars uncategorised ----
source('data/port_utils.R')
gruber <- 5/(sqrt(nrow(dat2))*log(nrow(dat2)))
port("A", cov.quanti = c("L1","L2","L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10"),
     cov.quali = c("L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19", "L20"),
     data = dat2, alpha = 0.05, beta = 0.1, gamma = 3)
# no longer work with a=0.02/0.03/0.04!
# but g=3: undet for a=0.01/0.025/0.05, det for a= 0.1
# (g= 2,3: undet for all b = 0.1 & a=0.01/0.02/gruber/0.05 (too small alpha lets algo focus on smaller strata?), det: a=0.1, 0.04, 0.03)

a_values <- c(0.01, 0.025, 0.05, 0.1)
b_values <- c(0.01, 5/(sqrt(nrow(dat2))*log(nrow(dat2))), 0.05, 0.1)
g_values <- 1:5
lst2 <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst2[[paste0("gamma = ", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A",cov.quanti = c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10"),
             cov.quali = c("L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19", "L20"),
             data = dat2, alpha = a, beta = b, gamma = g)
    }
  }
}
lst2
#sink("port_20_corr_uncat.txt")
# g = 1: only viol with cont conf (greedy cat), viol not yet poss bc intersection of 2
# g = 2-5: only det from a = 0.1 (& b=0.1) although should be for any a with sample prop = 22.7%
#          i.e. suspected focus on viol of smaller strata for smaller a? 



## PoRT: continuous vars categorised ----

dat2_cat <- dat2
for (i in names(dat2)[1:10]) {
  dat2_cat[[i]] <- cut(dat2[[i]], breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14))
}
dat2_cat
g_values <- 1:5
lst2_cat <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst2_cat[[paste0("gamma = ", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A", cov.quanti = NULL, cov.quali = c("L1", "L2", "L3", "L4", "L5",
                                                       "L6", "L7", "L8", "L9", "L10",
                                                       "L11", "L12", "L13", "L14", "L15",
                                                       "L16", "L17", "L18", "L19", "L20"),
             data = dat2_cat, alpha = a, beta = b, gamma = g)
    }
  }
}
lst2_cat
#sink("output_port/port_20_corr_cat.txt")
# g = 1: our viol not poss to find yet
# g=2: now alr det from a = 0.05 (& b=0.1), so still not for all a as should but better than uncategorised


## kbsd ----
source("kbsd.R")
o2 <- dat2
o2_1 <- o2
o2_1$A <- 1
o2_2 <- o2
o2_2$A <- 0
res2 <- kbsd(data = o2,
             int_data_list = list(o2_1, o2_2),
             disthalf_vec=c(L1=sd(o2$L1), L2 = sd(o2$L2), L3 = sd(o2$L3), L4 = sd(o2$L4), L5 = sd(o2$L5),
                            L6 = sd(o2$L6), L7 = sd(o2$L7), L8 = sd(o2$L8), L9 = sd(o2$L9), L10 = sd(o2$L10),
                            L11=sd(o2$L11), L12 = sd(o2$L12), L13 = sd(o2$L13), L14 = sd(o2$L14), L15 = sd(o2$L15),
                            L16 = sd(o2$L16), L17 = sd(o2$L17), L18 = sd(o2$L18), L19 = sd(o2$L19), L20 = sd(o2$L20),
                            A=0.5*sd(o2$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
res2_plot <- kbsd(data = o2,
                  int_data_list = list(o2_1, o2_2),
                  disthalf_vec=c(L1=sd(o2$L1), L2 = sd(o2$L2), L3 = sd(o2$L3), L4 = sd(o2$L4), L5 = sd(o2$L5),
                                 L6 = sd(o2$L6), L7 = sd(o2$L7), L8 = sd(o2$L8), L9 = sd(o2$L9), L10 = sd(o2$L10),
                                 L11=sd(o2$L11), L12 = sd(o2$L12), L13 = sd(o2$L13), L14 = sd(o2$L14), L15 = sd(o2$L15),
                                 L16 = sd(o2$L16), L17 = sd(o2$L17), L18 = sd(o2$L18), L19 = sd(o2$L19), L20 = sd(o2$L20),
                                 A=0.5*sd(o2$A)))
res2_plot  # slightly less support for A=1
table(dat2$A)

# what strata are those with low support for treated (IV = 1)
subset_2_1 <- res2[res2$diagnostic <= quantile(res2[res2$shift == 1, "diagnostic"], 0.25) & res2$shift == 1,]  # many obs with v low EDP
table(o2[subset_2_1$observation, c("L19", "L20")])
# those with few support in IV=1 are L19=0 & L20=0

# what strata are those with low support for treated (IV = 0): should be L19=1 & L20=1
subset_2_2 <- res2[res2$diagnostic < quantile(res2[res2$shift == 2, "diagnostic"], 0.25) & res2$shift == 2,]
table(o2[subset_2_2$observation, c("L19", "L20")])
# highest count among IV=0 is L19=1 & L20=1 -> viol found

# check if viol stratum found:
outliers_ind0 <- unique(res2[res2$shift == 2 & res2$diagnostic < quantile(res2[res2$shift == 2, "diagnostic"], 0.05), "observation"])
o2[outliers_ind0,] %>% nrow()
o2[outliers_ind0,] %>% filter(L19==1 & L20==1) %>% nrow()  # as for uncorr: majority of outliers are from critical stratum

# essence: PoRT & kbsd found viol, but PoRT did not always when should (cat made it slightly better, but still undet for smaller a, however smaller a also not that imp maybe)

# hm 
res2 <- kbsd(data = o2, int_data_list = list(o2_1, o2_2), type = "harmonicmean",
             disthalf_vec=c(L1=sd(o2$L1), L2 = sd(o2$L2), L3 = sd(o2$L3), L4 = sd(o2$L4), L5 = sd(o2$L5),
                            L6 = sd(o2$L6), L7 = sd(o2$L7), L8 = sd(o2$L8), L9 = sd(o2$L9), L10 = sd(o2$L10),
                            L11=sd(o2$L11), L12 = sd(o2$L12), L13 = sd(o2$L13), L14 = sd(o2$L14), L15 = sd(o2$L15),
                            L16 = sd(o2$L16), L17 = sd(o2$L17), L18 = sd(o2$L18), L19 = sd(o2$L19), L20 = sd(o2$L20),
                            A=0.5*sd(o2$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
#saveRDS(res2, "hm_20_uncat_corr.RDS")
# checking if viol stratum L19=1 & L20=1 included in the outliers
outliers_ind0 <- unique(res2[res2$shift == 2 & res2$diagnostic < quantile(res2[res2$shift == 2, "diagnostic"], 0.05), "observation"])
o2[outliers_ind0,] %>% nrow()
o2[outliers_ind0,] %>% filter(L19==1 & L20==1) %>% nrow()  # majority of outliers are from critical stratum




## KBSD cat ----
table(dat2_cat$L1)
table(dat2_cat$L2)  # all have same categorisation -> give numerical repr to compute kbsd values
dat2_cat <- dat2_cat %>%
  mutate(across(all_of(c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10")), 
                ~ case_when(
                  . == "(-3,-2]" ~ 1, . == "(-2,-1]"    ~ 2, . == "(-1,0]"    ~ 3,
                  . == "(0,1]"     ~ 4, . == "(1,2]"     ~ 5, . == "(2,3]"     ~ 6,
                  . == "(3,4]"     ~ 7, . == "(4,5]"     ~ 8, . == "(5,6]"     ~ 9,
                  . == "(6,7]"     ~ 10, . == "(7,8]"    ~ 11, . == "(8,9]"  ~ 12,
                  . == "(9,10]"     ~ 13, . == "(10,11]"    ~ 13, . == "(11,12]"  ~ 15,
                  . == "(12,13]"    ~ 16, . == "(13,14]"  ~ 17)))
source("kbsd.R")
o2 <- dat2_cat
o2_1 <- o2
o2_1$A <- 1
o2_2 <- o2
o2_2$A <- 0
res5 <- kbsd(data = o2,
             int_data_list = list(o2_1, o2_2),
             disthalf_vec=c(L1=sd(o2$L1), L2 = sd(o2$L2), L3 = sd(o2$L3), L4 = sd(o2$L4), L5 = sd(o2$L5),
                            L6 = sd(o2$L6), L7 = sd(o2$L7), L8 = sd(o2$L8), L9 = sd(o2$L9), L10 = sd(o2$L10),
                            L11=sd(o2$L11), L12 = sd(o2$L12), L13 = sd(o2$L13), L14 = sd(o2$L14), L15 = sd(o2$L15),
                            L16 = sd(o2$L16), L17 = sd(o2$L17), L18 = sd(o2$L18), L19 = sd(o2$L19), L20 = sd(o2$L20),
                            A=0.5*sd(o2$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
res5_plot <- kbsd(data = o2,
                  int_data_list = list(o2_1, o2_2),
                  disthalf_vec=c(L1=sd(o2$L1), L2 = sd(o2$L2), L3 = sd(o2$L3), L4 = sd(o2$L4), L5 = sd(o2$L5),
                                 L6 = sd(o2$L6), L7 = sd(o2$L7), L8 = sd(o2$L8), L9 = sd(o2$L9), L10 = sd(o2$L10),
                                 L11=sd(o2$L11), L12 = sd(o2$L12), L13 = sd(o2$L13), L14 = sd(o2$L14), L15 = sd(o2$L15),
                                 L16 = sd(o2$L16), L17 = sd(o2$L17), L18 = sd(o2$L18), L19 = sd(o2$L19), L20 = sd(o2$L20),
                                 A=0.5*sd(o2$A)))
res5_plot  # overall v few EDP as high-dim adjustment set
table(o2$A)  # a few less obs for A=1 matches that IV=1 has less support

# what strata are those with low support for treated (IV = 1)
subset_2_1 <- res5[res5$diagnostic < median(res5[res5$shift == 1, "diagnostic"]) & res5$shift == 1,]
table(o2[subset_2_1$observation, c("L19", "L20")]) # no viol expected

# what strata are those with low support for treated (IV = 0): should be L19=1 & L20=1, bc for these values will always get A=1, rarely A=0
subset_2_2 <- res5[res5$diagnostic < median(res5[res5$shift == 2, "diagnostic"]) & res5$shift == 2,]
table(o2[subset_2_2$observation, ][, c("L19", "L20")]) # highest count among IV=0 is L19=1 & L20=1 -> as expected

# check if viol stratum found:
outliers_ind0 <- unique(res2[res2$shift == 2 & res2$diagnostic < quantile(res2[res2$shift == 2, "diagnostic"], 0.05), "observation"])
o2[outliers_ind0,] %>% nrow()
o2[outliers_ind0,] %>% filter(L19==1 & L20==1) %>% nrow()  # same as for uncorr, same as for uncat: majority of outliers are from critical stratum

# essence: kbsd det viol also after cat

# hm
res2 <- kbsd(data = o2, int_data_list = list(o2_1, o2_2), type = "harmonicmean",
             disthalf_vec=c(L1=sd(o2$L1), L2 = sd(o2$L2), L3 = sd(o2$L3), L4 = sd(o2$L4), L5 = sd(o2$L5),
                            L6 = sd(o2$L6), L7 = sd(o2$L7), L8 = sd(o2$L8), L9 = sd(o2$L9), L10 = sd(o2$L10),
                            L11=sd(o2$L11), L12 = sd(o2$L12), L13 = sd(o2$L13), L14 = sd(o2$L14), L15 = sd(o2$L15),
                            L16 = sd(o2$L16), L17 = sd(o2$L17), L18 = sd(o2$L18), L19 = sd(o2$L19), L20 = sd(o2$L20),
                            A=0.5*sd(o2$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
#saveRDS(res2, "hm_20_cat_corr.RDS")
# checking if viol stratum L19=1 & L20=1 included in the outliers
outliers_ind0 <- unique(res2[res2$shift == 2 & res2$diagnostic < quantile(res2[res2$shift == 2, "diagnostic"], 0.05), "observation"])
o2[outliers_ind0,] %>% nrow()
o2[outliers_ind0,] %>% filter(L19==1 & L20==1) %>% nrow()  # majority of outliers are from critical stratum




# 20 Confounders With Middle Gap ----
set.seed(23092025)
L5_1 <- rnorm(500, 3, 1)
L5_2 <- rnorm(500, 7, 1)
L5_1_2 <- c(L5_1, L5_2)
plot(density(L5_1_2))
DAG3 <- DAG.empty() +
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
DAG3 <- set.DAG(DAG3)
dat3 <- sim(DAG3, rndseed = 12082025, n = 1000)[-1]
table(dat3$A)  # balanced

dat3 %>% filter(L5 >= 4 & L5 <= 6 & A==1) %>% nrow()/
  dat3 %>% filter(L5 >= 4 & L5 <= 6) %>% nrow()  # P(A=1)~0 -> to be det for g>=1, any b & any a (sample prop= 17%)

# table to check for all combos of binary conf if there are any with extreme P(A): fun defined in setup.R
binary_vars <- paste0("L", c(6:10))
tab <- make_strata_table(dat3, A = "A", binary_vars = binary_vars)
tab %>% filter((proba_exp <= 0.1 | proba_exp >= 0.9) & sample_prop >= 0.01) %>% print(n=64)
# no viol among binary vars -> all too small sample prop


## PoRT: continuous vars uncategorised ----
source('data/port_utils.R')
a_values <- c(0.01, 0.025, 0.05, 0.1)
b_values <- c(0.01, 5/(sqrt(nrow(dat3))*log(nrow(dat3))), 0.05, 0.1)
g_values <- 1:5
lst3 <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst3[[paste0("gamma = ", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A",cov.quanti = c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10"),
             cov.quali = c("L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19", "L20"),
             data = dat3, alpha = a, beta = b, gamma = g)
    }
  }
}
lst3
#sink("output_port/port_20_bimodal_uncat.txt")
# g = 1-5: viol det for all a & b, most of other viol only involve cont conf bc
#          due to greedy cat of cont conf (but also usually small sample size)


## PoRT: continuous vars categorised ----

# straightforward (preciser) categorisation
dat3_cat <- dat3
for (i in names(dat3)[1:10]) {
  dat3_cat[[i]] <- cut(dat3[[i]], breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14))
}
g_values <- 1:5
lst3_cat <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst3_cat[[paste0("gamma = ", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A", cov.quanti = NULL, cov.quali = c("L1", "L2", "L3", "L4", "L5",
                                                       "L6", "L7", "L8", "L9", "L10",
                                                       "L11", "L12", "L13", "L14", "L15",
                                                       "L16", "L17", "L18", "L19", "L20"),
             data = dat3_cat, alpha = a, beta = b, gamma = g)
    }
  }
}
lst3_cat
#sink("output_port/port_20_bimodal_cat.txt")
# g = 1-5: also always det but only as "L5=(-1,0],(4,5],(5,6]" -> tried to make stratum as large as poss?


## kbsd ----
source("kbsd.R")
set.seed(23092025)
o3 <- dat3
o3_1 <- o3
o3_1$A <- 1
o3_2 <- o3
o3_2$A <- 0
res3 <- kbsd(data = o3,
             int_data_list = list(o3_1, o3_2),
             disthalf_vec=c(L1=sd(o3$L1), L2 = sd(o3$L2), L3 = sd(o3$L3), L4 = sd(o3$L4), L5 = sd(o3$L5),
                            L6 = sd(o3$L6), L7 = sd(o3$L7), L8 = sd(o3$L8), L9 = sd(o3$L9), L10 = sd(o3$L10),
                            L11=sd(o3$L11), L12 = sd(o3$L12), L13 = sd(o3$L13), L14 = sd(o3$L14), L15 = sd(o3$L15),
                            L16 = sd(o3$L16), L17 = sd(o3$L17), L18 = sd(o3$L18), L19 = sd(o3$L19), L20 = sd(o3$L20),
                            A=0.5*sd(o3$A)), plot.out = F)
res3_plot <- kbsd(data = o3,
                  int_data_list = list(o3_1, o3_2),
                  disthalf_vec=c(L1=sd(o3$L1), L2 = sd(o3$L2), L3 = sd(o3$L3), L4 = sd(o3$L4), L5 = sd(o3$L5),
                                 L6 = sd(o3$L6), L7 = sd(o3$L7), L8 = sd(o3$L8), L9 = sd(o3$L9), L10 = sd(o3$L10),
                                 L11=sd(o3$L11), L12 = sd(o3$L12), L13 = sd(o3$L13), L14 = sd(o3$L14), L15 = sd(o3$L15),
                                 L16 = sd(o3$L16), L17 = sd(o3$L17), L18 = sd(o3$L18), L19 = sd(o3$L19), L20 = sd(o3$L20),
                                 A=0.5*sd(o3$A)))  # use 1 SD for L_i, 0.5 SD for A
res3_plot  # IV=2 (A=0) has slightly more support, but again EDP range in [0,1] diff to interpret due to so many dims
table(dat3$A)  # also more obs for A=0

# what strata are those with low support for treated (IV = 1)
shift1 <- res3[res3$shift == 1,]
outliers1 <- shift1$diagnostic < quantile(shift1$diagnostic, probs = .25)  # create indices for the "outliers"
for (i in names(o3)[-21]) {
  l_values1 <- dat3[outliers1, i]
  diag_values1 <- shift1[outliers1,]
  plot(l_values1, diag_values1$diagnostic, xlab = paste0("Values of Confounder ", i), ylab = "EDP")
}  
# no clear pattern for L5 in [4,6] being critical here


# what strata are those with low support for untreated (IV = 0)
shift2 <- res3[res3$shift == 2,]
outliers2 <- shift2$diagnostic < quantile(shift2$diagnostic, probs = .25)
for (i in names(o3)[-21]) {
  l_values2 <- dat3[outliers2, i]
  diag_values2 <- shift2[outliers2,]
  plot(l_values2, diag_values2$diagnostic, xlab = paste0("Values of Confounder ", i), ylab = "EDP")
}
# no obs with few EDP for L5<4 bc all obs of L5<4 have high EDP values as all have L5<4 & A=0

# check if viol stratum found:
outliers_ind1 <- unique(res3[res3$shift == 1 & res3$diagnostic < quantile(res3[res3$shift == 1, "diagnostic"], 0.05), "observation"])
o3[outliers_ind1,] %>% nrow()
o3[outliers_ind1,] %>% filter(L5 > 4 & L5 <6) %>% nrow()  # as for uncorr: majority of outliers are from critical stratum

# essence: for cat case, PoRT is not precise anymore (agg viol stratum with other category, 
# maybe to make stratum as large as poss?), also kbsd did not detect that few support in L5=[4,6]


# hm 
res3 <- kbsd(data = o3, int_data_list = list(o3_1, o3_2), type = "harmonicmean",
             disthalf_vec=c(L1=sd(o3$L1), L2 = sd(o3$L2), L3 = sd(o3$L3), L4 = sd(o3$L4), L5 = sd(o3$L5),
                            L6 = sd(o3$L6), L7 = sd(o3$L7), L8 = sd(o3$L8), L9 = sd(o3$L9), L10 = sd(o3$L10),
                            L11=sd(o3$L11), L12 = sd(o3$L12), L13 = sd(o3$L13), L14 = sd(o3$L14), L15 = sd(o3$L15),
                            L16 = sd(o3$L16), L17 = sd(o3$L17), L18 = sd(o3$L18), L19 = sd(o3$L19), L20 = sd(o3$L20),
                            A=0.5*sd(o3$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = T)
#ggsave("bp20_hm_uncat_int.png", width=2.5, height=6)
#saveRDS(res3, "hm_20_uncat_int.RDS")
# checking if viol stratum included in the outliers
outliers_ind0 <- unique(res3[res3$shift == 1 & res3$diagnostic < quantile(res3[res3$shift == 1, "diagnostic"], 0.05), "observation"])
o3[outliers_ind0,] %>% nrow()
o3[outliers_ind0,] %>% filter(L5 > 4 & L5 < 6) %>% nrow()  # majority of outliers are from critical stratum



## KBSD cat ----
table(dat3_cat$L1)
table(dat3_cat$L2)  # all have same categorisation -> give numerical repr to compute kbsd values
dat3_cat <- dat3_cat %>%
  mutate(across(all_of(c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10")), 
                ~ case_when(
                  . == "(-3,-2]" ~ 1, . == "(-2,-1]"    ~ 2, . == "(-1,0]"    ~ 3,
                  . == "(0,1]"     ~ 4, . == "(1,2]"     ~ 5, . == "(2,3]"     ~ 6,
                  . == "(3,4]"     ~ 7, . == "(4,5]"     ~ 8, . == "(5,6]"     ~ 9,
                  . == "(6,7]"     ~ 10, . == "(7,8]"    ~ 11, . == "(8,9]"  ~ 12,
                  . == "(9,10]"     ~ 13, . == "(10,11]"    ~ 13, . == "(11,12]"  ~ 15,
                  . == "(12,13]"    ~ 16, . == "(13,14]"  ~ 17)))
summary(dat3_cat)
source("kbsd.R")
o3 <- dat3_cat
o3_1 <- o3
o3_1$A <- 1
o3_2 <- o3
o3_2$A <- 0
res3 <- kbsd(data = o3,
             int_data_list = list(o3_1, o3_2),
             disthalf_vec=c(L1=sd(o3$L1), L2 = sd(o3$L2), L3 = sd(o3$L3), L4 = sd(o3$L4), L5 = sd(o3$L5),
                            L6 = sd(o3$L6), L7 = sd(o3$L7), L8 = sd(o3$L8), L9 = sd(o3$L9), L10 = sd(o3$L10),
                            L11=sd(o3$L11), L12 = sd(o3$L12), L13 = sd(o3$L13), L14 = sd(o3$L14), L15 = sd(o3$L15),
                            L16 = sd(o3$L16), L17 = sd(o3$L17), L18 = sd(o3$L18), L19 = sd(o3$L19), L20 = sd(o3$L20),
                            A=0.5*sd(o3$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
res3_plot <- kbsd(data = o3,
                  int_data_list = list(o3_1, o3_2),
                  disthalf_vec=c(L1=sd(o3$L1), L2 = sd(o3$L2), L3 = sd(o3$L3), L4 = sd(o3$L4), L5 = sd(o3$L5),
                                 L6 = sd(o3$L6), L7 = sd(o3$L7), L8 = sd(o3$L8), L9 = sd(o3$L9), L10 = sd(o3$L10),
                                 L11=sd(o3$L11), L12 = sd(o3$L12), L13 = sd(o3$L13), L14 = sd(o3$L14), L15 = sd(o3$L15),
                                 L16 = sd(o3$L16), L17 = sd(o3$L17), L18 = sd(o3$L18), L19 = sd(o3$L19), L20 = sd(o3$L20),
                                 A=0.5*sd(o3$A)))
res3_plot  # overall v few EDP as high-dim adjustment set
table(o3$A)  # a few less obs for A=1 matches that IV=1 has less support

# what strata are those with low support for treated (IV = 1)
subset_3_1 <- res3[res3$diagnostic < median(res3[res3$shift == 1, "diagnostic"]) & res3$shift == 1,]
table(o3[subset_3_1$observation, c("L3")]) # most from L3=[3,4] that have low support in IV=1 -> should be [4,6] tho

# what strata are those with low support for treated (IV = 0): should be L19=1 & L20=1
subset_3_2 <- res3[res3$diagnostic < median(res3[res3$shift == 2, "diagnostic"]) & res3$shift == 2,]
table(o3[subset_3_2$observation, ][, c("L3")]) # no viol expected

# check if viol stratum found:
outliers_ind1 <- unique(res3[res3$shift == 1 & res3$diagnostic < quantile(res3[res3$shift == 1, "diagnostic"], 0.05), "observation"])
o3[outliers_ind1,] %>% nrow()
o3[outliers_ind1,] %>% filter(L5 > 4 & L5 <6) %>% nrow()  # as for uncorr: majority of outliers are from critical stratum

# essence: after cat, kbsd did not find stratum with viol for A=1

# hm 
res3 <- kbsd(data = o3, int_data_list = list(o3_1, o3_2), type = "harmonicmean",
             disthalf_vec=c(L1=sd(o3$L1), L2 = sd(o3$L2), L3 = sd(o3$L3), L4 = sd(o3$L4), L5 = sd(o3$L5),
                            L6 = sd(o3$L6), L7 = sd(o3$L7), L8 = sd(o3$L8), L9 = sd(o3$L9), L10 = sd(o3$L10),
                            L11=sd(o3$L11), L12 = sd(o3$L12), L13 = sd(o3$L13), L14 = sd(o3$L14), L15 = sd(o3$L15),
                            L16 = sd(o3$L16), L17 = sd(o3$L17), L18 = sd(o3$L18), L19 = sd(o3$L19), L20 = sd(o3$L20),
                            A=0.5*sd(o3$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
#saveRDS(res3, "hm_20_cat_int.RDS")
# checking if viol stratum L19=1 & L20=1 included in the outliers
outliers_ind0 <- unique(res3[res3$shift == 2 & res3$diagnostic < quantile(res3[res3$shift == 2, "diagnostic"], 0.05), "observation"])
o3[outliers_ind0,] %>% nrow()
o3[outliers_ind0,] %>% filter(L19==1 & L20==1) %>% nrow()  # majority of outliers are from critical stratum



# 20 Confounders With Left-Side Gap ----

set.seed(23092025)
L5_1 <- rnorm(100, 3, 1)
L5_2 <- rnorm(900, 6, 1)
L5_1_2 <- c(L5_1, L5_2)
plot(density(L5_1_2))
DAG4 <- DAG.empty() +
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
DAG4 <- set.DAG(DAG4)
dat4 <- sim(DAG4, rndseed = 12082025, n = 1000)[-1]
table(dat4$A)  # balanced
dat4 %>% filter(L5 <4 & A==1) %>% nrow()/
  dat4 %>% filter(L5 <4) %>% nrow()  # P(A=1)~0 -> to be det for g>=1, b=0.1, a <= 0.05 (sample prop= 9.9%)

# table to check for all combos of binary conf if there are any with extreme P(A): fun defined in setup.R
binary_vars <- paste0("L", c(6:10))
tab <- make_strata_table(dat4, A = "A", binary_vars = binary_vars)
tab %>% filter((proba_exp <= 0.1 | proba_exp >= 0.9) & sample_prop >= 0.01) %>% print(n=64)
# no viol among binary vars -> all too small sample prop


## PoRT: continuous vars uncategorised ----
source('data/port_utils.R')
a_values <- c(0.01, 0.025, 0.05, 0.1)
b_values <- c(0.01, 5/(sqrt(nrow(dat4))*log(nrow(dat4))), 0.05, 0.1)
g_values <- 1:5
lst4 <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst4[[paste0("gamma = ", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A",cov.quanti = c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10"),
             cov.quali = c("L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19", "L20"),
             data = dat4, alpha = a, beta = b, gamma = g)
    }
  }
}
lst4
#sink("output_port/port_20_leftgap_uncat.txt")
# g = 1-5: always det when should (a<=0.05 & b=0.1), for smaller b only as intersection with other cont vars
# g = 5, a = 0.1, b = 0.1: largest viol subgroup with 12.7% sample prop which was not expected


## PoRT: continuous vars categorised ----

# straightforward (preciser) categorisation
dat4_cat <- dat4
for (i in names(dat4)[1:10]) {
  dat4_cat[[i]] <- cut(dat4[[i]], breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14))
}
g_values <- 1:5
lst4_cat <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst4_cat[[paste0("gamma = ", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A", cov.quanti = NULL, cov.quali = c("L1", "L2", "L3", "L4", "L5",
                                                       "L6", "L7", "L8", "L9", "L10",
                                                       "L11", "L12", "L13", "L14", "L15",
                                                       "L16", "L17", "L18", "L19", "L20"),
             data = dat4_cat, alpha = a, beta = b, gamma = g)
    }
  }
}
lst4_cat
#sink("output_port/port_20_leftgap_cat.txt")
# g = 1-5: only ever as "L5=(1,2],(2,3],(3,4],(9,10]" which enables det for b>=0.05, but not precise enough


## kbsd ----
source("kbsd.R")
set.seed(23092025)
o4 <- dat4
o4_1 <- o4
o4_1$A <- 1
o4_2 <- o4
o4_2$A <- 0
res4 <- kbsd(data = o4,
             int_data_list = list(o4_1, o4_2),
             disthalf_vec=c(L1=sd(o4$L1), L2 = sd(o4$L2), L3 = sd(o4$L3), L4 = sd(o4$L4), L5 = sd(o4$L5),
                            L6 = sd(o4$L6), L7 = sd(o4$L7), L8 = sd(o4$L8), L9 = sd(o4$L9), L10 = sd(o4$L10),
                            L11=sd(o4$L11), L12 = sd(o4$L12), L13 = sd(o4$L13), L14 = sd(o4$L14), L15 = sd(o4$L15),
                            L16 = sd(o4$L16), L17 = sd(o4$L17), L18 = sd(o4$L18), L19 = sd(o4$L19), L20 = sd(o4$L20),
                            A=0.5*sd(o4$A)), plot.out = F)
res4_plot <- kbsd(data = o4,
                  int_data_list = list(o4_1, o4_2),
                  disthalf_vec=c(L1=sd(o4$L1), L2 = sd(o4$L2), L3 = sd(o4$L3), L4 = sd(o4$L4), L5 = sd(o4$L5),
                                 L6 = sd(o4$L6), L7 = sd(o4$L7), L8 = sd(o4$L8), L9 = sd(o4$L9), L10 = sd(o4$L10),
                                 L11=sd(o4$L11), L12 = sd(o4$L12), L13 = sd(o4$L13), L14 = sd(o4$L14), L15 = sd(o4$L15),
                                 L16 = sd(o4$L16), L17 = sd(o4$L17), L18 = sd(o4$L18), L19 = sd(o4$L19), L20 = sd(o4$L20),
                                 A=0.5*sd(o4$A)))  # use 1 SD for L_i, 0.5 SD for A
res4_plot  # IV=1 has slightly more support, but EDP range in [0,1] diff to interpret due to so many dims
table(dat4$A)

# what strata are those with low support for treated (IV = 1)
shift1 <- res4[res4$shift == 1,]
outliers1 <- shift1$diagnostic < quantile(shift1$diagnostic, probs = .25)  # create indices for the "outliers"
for (i in names(o4)[-21]) {
  l_values1 <- dat4[outliers1, i]
  diag_values1 <- shift1[outliers1,]
  plot(l_values1, diag_values1$diagnostic, xlab = paste0("Values of Confounder ", i), ylab = "EDP")
}  
# no clear pattern except for vars L1, L5, L7 -> det that lowest EDP are for L5<4 as expected! 

# what strata are those with low support for untreated (IV = 0)
shift2 <- res4[res4$shift == 2,]
outliers2 <- shift2$diagnostic < quantile(shift2$diagnostic, probs = .25)
for (i in names(o4)[-21]) {
  l_values2 <- dat4[outliers2, i]
  diag_values2 <- shift2[outliers2,]
  plot(l_values2, diag_values2$diagnostic, xlab = paste0("Values of Confounder ", i), ylab = "EDP")
}
# also low, but more imp only FEW obs with few EDP for L5<4 bc few obs of L5<4 in this subset of having low support (most with L5<4 have high EDP values bc a lot with L5<4 & A=0)

outliers_ind1 <- unique(res4[res4$shift == 1 & res4$diagnostic < quantile(res4[res4$shift == 1, "diagnostic"], 0.05), "observation"])
o4[outliers_ind1,] %>% nrow()
o4[outliers_ind1,] %>% filter(L5 <=4) %>% nrow()  # as for uncorr: majority of outliers are from critical stratum

# essence: after cat, PoRT aggregates viol subgroup with other cat, i.e. not precise anymore

# hm 
res4 <- kbsd(data = o4, int_data_list = list(o4_1, o4_2), type = "harmonicmean",
             disthalf_vec=c(L1=sd(o4$L1), L2 = sd(o4$L2), L3 = sd(o4$L3), L4 = sd(o4$L4), L5 = sd(o4$L5),
                            L6 = sd(o4$L6), L7 = sd(o4$L7), L8 = sd(o4$L8), L9 = sd(o4$L9), L10 = sd(o4$L10),
                            L11=sd(o4$L11), L12 = sd(o4$L12), L13 = sd(o4$L13), L14 = sd(o4$L14), L15 = sd(o4$L15),
                            L16 = sd(o4$L16), L17 = sd(o4$L17), L18 = sd(o4$L18), L19 = sd(o4$L19), L20 = sd(o4$L20),
                            A=0.5*sd(o4$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = T)
#ggsave("bp20_hm_uncat_ext.png", width=2.5, height=6)
#saveRDS(res4, "hm_20_uncat_ext.RDS")
# checking if viol stratum L19=1 & L20=1 included in the outliers
outliers_ind0 <- unique(res4[res4$shift == 1 & res4$diagnostic < quantile(res4[res4$shift == 1, "diagnostic"], 0.05), "observation"])
o4[outliers_ind0,] %>% nrow()
o4[outliers_ind0,] %>% filter(L5 <4) %>% nrow()  # majority of outliers are from critical stratum



## KBSD cat ----
table(dat4_cat$L1)
table(dat4_cat$L9)  # all have same categorisation -> give numerical repr to compute kbsd values
dat4_cat <- dat4_cat %>%
  mutate(across(all_of(c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10")), 
                ~ case_when(
                  . == "(-3,-2]" ~ 1, . == "(-2,-1]"    ~ 2, . == "(-1,0]"    ~ 3,
                  . == "(0,1]"     ~ 4, . == "(1,2]"     ~ 5, . == "(2,3]"     ~ 6,
                  . == "(3,4]"     ~ 7, . == "(4,5]"     ~ 8, . == "(5,6]"     ~ 9,
                  . == "(6,7]"     ~ 10, . == "(7,8]"    ~ 11, . == "(8,9]"  ~ 12,
                  . == "(9,10]"     ~ 13, . == "(10,11]"    ~ 13, . == "(11,12]"  ~ 15,
                  . == "(12,13]"    ~ 16, . == "(13,14]"  ~ 17)))
source("kbsd.R")
o4 <- dat4_cat
o4_1 <- o4
o4_1$A <- 1
o4_2 <- o4
o4_2$A <- 0
res4 <- kbsd(data = o4,
             int_data_list = list(o4_1, o4_2),
             disthalf_vec=c(L1=sd(o4$L1), L2 = sd(o4$L2), L3 = sd(o4$L3), L4 = sd(o4$L4), L5 = sd(o4$L5),
                            L6 = sd(o4$L6), L7 = sd(o4$L7), L8 = sd(o4$L8), L9 = sd(o4$L9), L10 = sd(o4$L10),
                            L11=sd(o4$L11), L12 = sd(o4$L12), L13 = sd(o4$L13), L14 = sd(o4$L14), L15 = sd(o4$L15),
                            L16 = sd(o4$L16), L17 = sd(o4$L17), L18 = sd(o4$L18), L19 = sd(o4$L19), L20 = sd(o4$L20),
                            A=0.5*sd(o4$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
res4_plot <- kbsd(data = o4,
                  int_data_list = list(o4_1, o4_2),
                  disthalf_vec=c(L1=sd(o4$L1), L2 = sd(o4$L2), L3 = sd(o4$L3), L4 = sd(o4$L4), L5 = sd(o4$L5),
                                 L6 = sd(o4$L6), L7 = sd(o4$L7), L8 = sd(o4$L8), L9 = sd(o4$L9), L10 = sd(o4$L10),
                                 L11=sd(o4$L11), L12 = sd(o4$L12), L13 = sd(o4$L13), L14 = sd(o4$L14), L15 = sd(o4$L15),
                                 L16 = sd(o4$L16), L17 = sd(o4$L17), L18 = sd(o4$L18), L19 = sd(o4$L19), L20 = sd(o4$L20),
                                 A=0.5*sd(o4$A)))
res4_plot  # IV=1 with a bit more support
table(o4$A)  # also more obs for that

# what strata are those with low support for treated (IV = 1): should be L3<4
subset_4_1 <- res4[res4$diagnostic < median(res4[res4$shift == 1, "diagnostic"]) & res4$shift == 1,]
table(o4[subset_4_1$observation, c("L3")]) #  as expected; most from L3=[2,4] that have low support in IV=1

# what strata are those with low support for treated (IV = 0)
subset_4_2 <- res4[res4$diagnostic < median(res4[res4$shift == 2, "diagnostic"]) & res4$shift == 2,]
table(o4[subset_4_2$observation, ][, c("L3")]) # no viol expected

outliers_ind1 <- unique(res4[res4$shift == 1 & res4$diagnostic < quantile(res4[res4$shift == 1, "diagnostic"], 0.05), "observation"])
o4[outliers_ind1,] %>% nrow()
o4[outliers_ind1,] %>% filter(L5 <8) %>% nrow()  # as for uncorr: majority of outliers are from critical stratum

# essence: kbsd also det the critical stratum L3<4 after cat

# hm 
res4 <- kbsd(data = o4, int_data_list = list(o4_1, o4_2), type = "harmonicmean",
             disthalf_vec=c(L1=sd(o4$L1), L2 = sd(o4$L2), L3 = sd(o4$L3), L4 = sd(o4$L4), L5 = sd(o4$L5),
                            L6 = sd(o4$L6), L7 = sd(o4$L7), L8 = sd(o4$L8), L9 = sd(o4$L9), L10 = sd(o4$L10),
                            L11=sd(o4$L11), L12 = sd(o4$L12), L13 = sd(o4$L13), L14 = sd(o4$L14), L15 = sd(o4$L15),
                            L16 = sd(o4$L16), L17 = sd(o4$L17), L18 = sd(o4$L18), L19 = sd(o4$L19), L20 = sd(o4$L20),
                            A=0.5*sd(o4$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = F)
#saveRDS(res4, "hm_20_cat_ext.RDS")
# checking if viol stratum L19=1 & L20=1 included in the outliers
outliers_ind0 <- unique(res4[res4$shift == 2 & res4$diagnostic < quantile(res4[res4$shift == 2, "diagnostic"], 0.05), "observation"])
o4[outliers_ind0,] %>% nrow()
o4[outliers_ind0,] %>% filter(L19==1 & L20==1) %>% nrow()  # majority of outliers are from critical stratum





# 20 Confounders Binary ----

set.seed(23092025)
DAG5 <- DAG.empty() +
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
DAG5 <- set.DAG(DAG5)
dat5 <- sim(DAG5, rndseed = 12082025, n = 1000)[-1]
table(dat5$A)  # a bit imbalanced
dat5 %>% filter(L19==1 & L20 == 1 & A==1) %>% nrow()/
  dat5 %>% filter(L19==1 & L20 == 1) %>% nrow()  # P(A=1)~0 -> to be det for g>=2, b=0.1, any a (sample prop= 24%)

# table to check for all the binary conf if there are any with extreme P(A): fun defined in setup.R
binary_vars <- paste0("L", c(1:20))
tab <- make_strata_table(dat5, A = "A", binary_vars = binary_vars)
tab %>% filter((proba_exp <= 0.1 | proba_exp >= 0.9) & sample_prop >= 0.01) %>% print(n=64)
# no viol among binary vars -> all too small sample prop


## PoRT ----
source('data/port_utils.R')
a_values <- c(0.01, 0.025, 0.05, 0.1)
b_values <- c(0.01, 5/(sqrt(nrow(dat5))*log(nrow(dat5))), 0.05, 0.1)
g_values <- 1:5
lst5 <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      lst5[[paste0("gamma = ", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A = "A",cov.quanti = NULL,
             cov.quali = c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10",
                           "L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19", "L20"),
             data = dat5, alpha = a, beta = b, gamma = g)
    }
  }
}
lst5
#sink("port_20_binary.txt")
# g = 1: viol not poss yet bc intersection of 2
# g = 2-5: always det for any a & b=0.1, other large viol are subviol of L19=1 & L20=1 or smaller strata


## kbsd ----
source("kbsd.R")
o5 <- dat5
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

res5_plot
# overall v low EDP which makes sense as all obs are v spread out over 20 dims
# lower support for IV=1 (A=1) than for IV=2 (A=0)

# what strata are those with low support for treated (IV = 1): should find L19=1 & L20=1
subset_5_1 <- res5[res5$diagnostic <= quantile(res5[res5$shift == 1, "diagnostic"], 0.25) & res5$shift == 1,]  # many obs with v low EDP
table(o5[subset_5_1$observation, c("L19", "L20")])
# indeed most from those with few support in IV=1 are L19=1 & L20=1 -> viol found

# what strata are those with low support for treated (IV = 0):
subset_5_2 <- res5[res5$diagnostic < quantile(res5[res5$shift == 2, "diagnostic"], .25) & res5$shift == 2,]
o5[subset_5_2$observation, ]
table(o5[subset_5_2$observation, c("L19", "L20")]) # highest count among IV=0 is L19=1 & L20=0

outliers_ind1 <- unique(res5[res5$shift == 1 & res5$diagnostic < quantile(res5[res5$shift == 1, "diagnostic"], 0.05), "observation"])
o5[outliers_ind1,] %>% nrow()
o5[outliers_ind1,] %>% filter(L19==1 & L20==1) %>% nrow()  # as for uncorr: majority of outliers are from critical stratum

# essence: both PoRT & kbsd det viol in case where all conf are binary

# hm 
res5 <- kbsd(data = o5, int_data_list = list(o5_1, o5_2), type = "harmonicmean",
             disthalf_vec=c(L1=sd(o5$L1), L2 = sd(o5$L2), L3 = sd(o5$L3), L4 = sd(o5$L4), L5 = sd(o5$L5),
                            L6 = sd(o5$L6), L7 = sd(o5$L7), L8 = sd(o5$L8), L9 = sd(o5$L9), L10 = sd(o5$L10),
                            L11=sd(o5$L11), L12 = sd(o5$L12), L13 = sd(o5$L13), L14 = sd(o5$L14), L15 = sd(o5$L15),
                            L16 = sd(o5$L16), L17 = sd(o5$L17), L18 = sd(o5$L18), L19 = sd(o5$L19), L20 = sd(o5$L20),
                            A=0.5*sd(o5$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = T)
#ggsave("bp20_hm_binary.png", width=2.5, height=6)
#saveRDS(res5, "hm_20_uncat_binary.RDS")
# checking if viol stratum L19=1 & L20=1 included in the outliers
outliers_ind0 <- unique(res5[res5$shift == 1 & res5$diagnostic < quantile(res5[res5$shift == 1, "diagnostic"], 0.05), "observation"])
o5[outliers_ind0,] %>% nrow()
o5[outliers_ind0,] %>% filter(L19==1 & L20==1) %>% nrow()  # majority of outliers are from critical stratum


