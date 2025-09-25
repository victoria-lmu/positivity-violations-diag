source("setup.R")

# 5) 20 Confounders ----
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
# if one of L19 and L20 = 0, then P(A)~0.26, if both =1 then P(A)=211/227=0.93 (pos viol for b=0.1)
DAG2 <- set.DAG(DAG2)
dat2 <- sim(DAG2, rndseed = 12082025, n = 1000)[-1]
table(dat2$A)  # good that not too imbalanced
dat2 %>% filter(L19==1 & L20 == 1 & A==1) %>% nrow()/
  dat2 %>% filter(L19==1 & L20 == 1) %>% nrow()  # expect to be detected for b=0.1 from g=2 for any a (sample prop= 22.7%)



source('data/port_utils.R')

## PoRT: continuous vars uncategorised ----
gruber <- 5/(sqrt(nrow(dat2))*log(nrow(dat2)))
port("A", cov.quanti = c("L1","L2","L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10"),
     cov.quali = c("L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19", "L20"),
     data = dat2, alpha = 0.05, beta = 0.1, gamma = 3)
# gamma = 2: always b = 0.1, undetected: a=0.01/0.02/gruber/0.05 (too small alpha lets algo focus on smaller strata?), det: a=0.1, 0.04, 0.03
# gamma = 3: undet: a=0.01, 0.02, 0.05, gruber, det: a=0.03, 0.04, 0.1

## PoRT: continuous vars categorised ----

# straightforward (preciser) categorisation
dat2_cat <- dat2
for (i in names(dat2)[1:10]) {
  dat2_cat[[i]] <- cut(dat2[[i]], breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14))
}
dat2_cat
port("A", cov.quanti = NULL,
     cov.quali = c("L1","L2","L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10",
                   "L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19", "L20"),
     data = dat2_cat, alpha = 0.1, beta = 0.1, gamma = 3)
# gamma = 2: detected for all with beta = 0.1 except with a = 0.01 (too small -> as above, but alr better with categorisation)
# gamma = 3: same as gamma=2 

a_values <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.1)
b_values <- c(0.01, 5/(sqrt(nrow(dat2))*log(nrow(dat2))), 0.05, 0.1)
lst6_cat <- list()
for (a in a_values) {
  for (b in b_values) {
    lst6_cat[[paste0("alpha = ", a, ", beta = ", b)]] <- port(A = "A", cov.quanti = NULL,
                                                              cov.quali = c("L1", "L2", "L3", "L4", "L5",
                                                                            "L6", "L7", "L8", "L9", "L10",
                                                                            "L11", "L12", "L13", "L14", "L15",
                                                                            "L16", "L17", "L18", "L19", "L20"),
                                                              data = dat2_cat, alpha = a, beta = b, gamma = 20)
  }
}
lst5_cat

# broader categorisation
dat2_cat <- dat2
for (i in names(dat2)[1:10]) {
  dat2_cat[[i]] <- cut(dat2[[i]], breaks = c(-3, 0, 3, 6, 9, 12, 15))
}
dat2_cat
port("A", cov.quanti = NULL,
     cov.quali = c("L1","L2","L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10",
                   "L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19", "L20"),
     data = dat2_cat, alpha = 0.01, beta = 0.1, gamma = 2)
# gamma = 2: again for all except a=0.01, with beta = 0.1

# essence: imp insight that cat of cont vars impacts detection of viol in CATEGORICAL vars!
#          -> here, for broader cat higher chance of detection with all alpha > 0.01, without cat only for alpha <= 0.03


## kbsd ----
source("kbsd.R")
o6 <- dat2
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
# overall v low EDP which makes sense as all obs are v spread out over 20 dims
# -> many defined strata so that within each prob just one treatment level -> no support for the other
#    was warned for in kbsd paper
# -> but interesting that among IV=1 (bringing all to A=1), particularly low support (many with EDP = 0)


### alternative EDP calc for high-dim covar set ----

# type = "minval" instead of default type = "Rfast"
# requires specification of EDP minvalue accepted for each covar (more imp <=> lower/no minval, less imp <=> 0.5)
res6 <- kbsd(data = o6, int_data_list = list(o6_1, o6_2), type = "minval",
             disthalf_vec=c(L1=sd(o6$L1), L2 = sd(o6$L2), L3 = sd(o6$L3), L4 = sd(o6$L4), L5 = sd(o6$L5),
                            L6 = sd(o6$L6), L7 = sd(o6$L7), L8 = sd(o6$L8), L9 = sd(o6$L9), L10 = sd(o6$L10),
                            L11=sd(o6$L11), L12 = sd(o6$L12), L13 = sd(o6$L13), L14 = sd(o6$L14), L15 = sd(o6$L15),
                            L16 = sd(o6$L16), L17 = sd(o6$L17), L18 = sd(o6$L18), L19 = sd(o6$L19), L20 = sd(o6$L20),
                            A=0.5*sd(o6$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = T)
# type = "harmonicmean" instead of default type = "Rfast"
res6 <- kbsd(data = o6, int_data_list = list(o6_1, o6_2), type = "harmonicmean",
             disthalf_vec=c(L1=sd(o6$L1), L2 = sd(o6$L2), L3 = sd(o6$L3), L4 = sd(o6$L4), L5 = sd(o6$L5),
                            L6 = sd(o6$L6), L7 = sd(o6$L7), L8 = sd(o6$L8), L9 = sd(o6$L9), L10 = sd(o6$L10),
                            L11=sd(o6$L11), L12 = sd(o6$L12), L13 = sd(o6$L13), L14 = sd(o6$L14), L15 = sd(o6$L15),
                            L16 = sd(o6$L16), L17 = sd(o6$L17), L18 = sd(o6$L18), L19 = sd(o6$L19), L20 = sd(o6$L20),
                            A=0.5*sd(o6$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = T)






# 6) 50 Confounders ----
set.seed(23092025)
DAG3 <- DAG.empty() +
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
  node("L11", distr = "rnorm", mean = 11, sd = 1) +
  node("L12", distr = "rnorm", mean = 12, sd = 1) +
  node("L13", distr = "rnorm", mean = 13, sd = 1) +
  node("L14", distr = "rnorm", mean = 14, sd = 1) +
  node("L15", distr = "rnorm", mean = 15, sd = 1) +
  node("L16", distr = "rnorm", mean = 16, sd = 1) +
  node("L17", distr = "rnorm", mean = 17, sd = 1) +
  node("L18", distr = "rnorm", mean = 18, sd = 1) +
  node("L19", distr = "rnorm", mean = 19, sd = 1) +
  node("L20", distr = "rnorm", mean = 20, sd = 1) +
  node("L21", distr = "rnorm", mean = 21, sd = 1) +
  node("L22", distr = "rnorm", mean = 22, sd = 1) +
  node("L23", distr = "rnorm", mean = 23, sd = 1) +
  node("L24", distr = "rnorm", mean = 24, sd = 1) +
  node("L25", distr = "rnorm", mean = 25, sd = 1) +
  node("L26", distr = "rbern", prob = 0.5) +
  node("L27", distr = "rbern", prob = 0.5) +
  node("L28", distr = "rbern", prob = 0.5) +
  node("L29", distr = "rbern", prob = 0.5) +
  node("L30", distr = "rbern", prob = 0.5) +
  node("L31", distr = "rbern", prob = 0.5) +
  node("L32", distr = "rbern", prob = 0.5) +
  node("L33", distr = "rbern", prob = 0.5) +
  node("L34", distr = "rbern", prob = 0.5) +
  node("L35", distr = "rbern", prob = 0.5) +
  node("L36", distr = "rbern", prob = 0.5) +
  node("L37", distr = "rbern", prob = 0.5) +
  node("L38", distr = "rbern", prob = 0.5) +
  node("L39", distr = "rbern", prob = 0.5) +
  node("L40", distr = "rbern", prob = 0.5) +
  node("L41", distr = "rbern", prob = 0.5) +
  node("L42", distr = "rbern", prob = 0.5) +
  node("L43", distr = "rbern", prob = 0.5) +
  node("L44", distr = "rbern", prob = 0.5) +
  node("L45", distr = "rbern", prob = 0.5) +
  node("L46", distr = "rbern", prob = 0.5) +
  node("L47", distr = "rbern", prob = 0.5) +
  node("L48", distr = "rbern", prob = 0.5) +
  node("L49", distr = "rbern", prob = 0.5) +
  node("L50", distr = "rbern", prob = 0.5) +
  node("A", distr = "rbern", prob = plogis(L1 - L2 + L5*L49*L50))
  # simulate so that if one of L49&L50 =0 -> P(A)=0.27, if both =1 -> P(A)~0.98 (expected viol!)
DAG3 <- set.DAG(DAG3)
dat3 <- sim(DAG3, rndseed = 12082025, n = 1000)[-1]
table(dat3$A)  # balanced

dat3 %>% filter(L49==1 & L50==1 & A==1) %>% nrow()/
  dat3 %>% filter(L49==1 & L50==1) %>% nrow()  # viol to be found with b = 0.05/0.1, any alpha as sample prop = 27%


## PoRT: continuous vars uncategorised ----
port("A", cov.quanti = c("L1","L2","L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10",
                         "L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19", "L20",
                         "L21", "L22", "L23", "L24" ,"L25"),
     cov.quali = c("L26", "L27", "L28", "L29", "L30", "L31", "L32", "L33", "L34", "L35", "L36",
                   "L37", "L38", "L39", "L40", "L41", "L42", "L43", "L44", "L45", "L46", "L47",
                   "L48", "L49", "L50"),
     data = dat3, alpha = 0.1, beta = 0.1, gamma = 3)
# gamma = 2: undetected: a=0.02, gruber, 0.03, 0.04, 0.05, det: a = 0.1, 0.01, special: returned "root(3.8%)" for a=0.03
# no pattern -> maybe really due to random fluctuations? only pattern overall is
#    that for constant b, g, and decreasing a (smaller a), more strata bc if continuous can define them v v small
# gamma = 3: undet: 0.02, gruber, 0.03, 0.04, 0.05, det: 0.01, 0.1


## PoRT: continuous vars categorised ----

# straightforward (preciser) categorisation
dat3_cat <- dat3
for (i in names(dat3)[1:25]) {
  dat3_cat[[i]] <- cut(dat3[[i]], breaks = c(-3, 0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30))
}
dat3_cat
port("A", cov.quanti = c("L1","L2","L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10",
                         "L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19", "L20",
                         "L21", "L22", "L23", "L24" ,"L25"),
     cov.quali = c("L26", "L27", "L28", "L29", "L30", "L31", "L32", "L33", "L34", "L35", "L36",
                   "L37", "L38", "L39", "L40", "L41", "L42", "L43", "L44", "L45", "L46", "L47",
                   "L48", "L49", "L50"),
     data = dat3_cat, alpha = 0.1, beta = 0.1, gamma = 3)
# gamma = 2: undet: 0.01, det: 0.02, gruber, 0.03, 0.04, 0.05, 0.1
# gamma = 3: undet: 0.01, det: 0.02, gruber, 0.03, 0.04, 0.05, 0.1

# broader categorisation
dat3_cat <- dat3
for (i in names(dat3)[1:25]) {
  dat3_cat[[i]] <- cut(dat3[[i]], breaks = c(-3, 6, 15, 24, 30))
}
dat3_cat
port("A", cov.quanti = c("L1","L2","L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10",
                         "L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19", "L20",
                         "L21", "L22", "L23", "L24" ,"L25"),
     cov.quali = c("L26", "L27", "L28", "L29", "L30", "L31", "L32", "L33", "L34", "L35", "L36",
                   "L37", "L38", "L39", "L40", "L41", "L42", "L43", "L44", "L45", "L46", "L47",
                   "L48", "L49", "L50"),
     data = dat3_cat, alpha = 0.1, beta = 0.1, gamma = 2)
# gamma = 2: undet: 0.01, det: 0.02-0.1
# gamma = 3: undet: 0.01, det: 0.02-0.1


## kbsd ----
source("kbsd.R")
o6 <- dat3
o6_1 <- o6 %>% mutate(A=1)
o6_2 <- o6 %>% mutate(A=0)
res6 <- kbsd(data = dat3, int_data_list = list(o6_1, o6_2), type = "Rfast",
             disthalf_vec=c(L1=sd(o6$L1), L2 = sd(o6$L2), L3 = sd(o6$L3), L4 = sd(o6$L4), L5 = sd(o6$L5),
                            L6 = sd(o6$L6), L7 = sd(o6$L7), L8 = sd(o6$L8), L9 = sd(o6$L9), L10 = sd(o6$L10),
                            L11=sd(o6$L11), L12 = sd(o6$L12), L13 = sd(o6$L13), L14 = sd(o6$L14), L15 = sd(o6$L15),
                            L16 = sd(o6$L16), L17 = sd(o6$L17), L18 = sd(o6$L18), L19 = sd(o6$L19), L20 = sd(o6$L20),
                            L21=sd(o6$L21), L22 = sd(o6$L22), L23 = sd(o6$L23), L24 = sd(o6$L24), L25 = sd(o6$L25),
                            L26 = sd(o6$L26), L27 = sd(o6$L27), L28 = sd(o6$L28), L29 = sd(o6$L29), L30 = sd(o6$L30),
                            L31=sd(o6$L31), L32 = sd(o6$L32), L33 = sd(o6$L33), L34 = sd(o6$L34), L35 = sd(o6$L35),
                            L36 = sd(o6$L36), L37 = sd(o6$L37), L38 = sd(o6$L38), L39 = sd(o6$L39), L40 = sd(o6$L40),
                            L41=sd(o6$L41), L42 = sd(o6$L42), L43 = sd(o6$L43), L44 = sd(o6$L44), L45 = sd(o6$L45),
                            L46 = sd(o6$L46), L47 = sd(o6$L47), L48 = sd(o6$L48), L49 = sd(o6$L49), L50 = sd(o6$L50),
                            A=0.5*sd(o6$A)),  # use 1 SD for L_i, 0.5 SD for A
             plot.out = T)
# as for setting with 20 confounders, only few EDP always due to high number of dims

### alternative EDP calc for high-dim covar set ----
# type = "minval" instead of default type = "Rfast"
# type = "harmonicmean" instead of default type = "Rfast"
