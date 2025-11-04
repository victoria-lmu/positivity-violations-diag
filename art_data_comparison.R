source("setup.R")

# 1) Data Cleaning ----
setwd("C:/Users/victo/OneDrive/Desktop/Uni/Statistik/BA_thesis/positivity-violations-diag")
sdd <- read.csv("data/sameday_wide_MI.csv") # read in .csv file with imputed data
sameday <- subset(sdd, select=c("sexpregnancy","age_cat2","education_cat1","mstatus_cat1",
                                "year","facility_cat2_rename","HIVdiagEnrol1","preEAAA1",
                                "whostudybase","cd4studybase","tb_atARtinit_adapted",
                                "bmi","hb","alt","creat","phone","sameDayEn14",
                                "outcome12m","outcome18m","outcome24m","X_mi_m"))

# select according to DAG
colnames(sameday)  <- c("sex","age","education","marital","year","facility","TimeHIVToEnrol",  # rename diff for clarity *V*
                        "underTreatAll","who","cd4","tb","bmi","hb","alt","creat","phone",
                        "SAMEDAY","outcome12m","outcome18m","outcome24m","X_mi_m")
sameday <- sameday[sameday$X_mi_m!=0,]  # keep only 10 out of all 11 entries per person *V*
sameday$sex        <- factor(sameday$sex)
sameday$age        <- factor(sameday$age)
sameday$education  <- factor(sameday$education)
sameday$maritial   <- factor(sameday$marital)
sameday$year   <- factor(sameday$year)
sameday$facility   <- factor(sameday$facility)
sameday$TimeHIVToEnrol   <- factor(sameday$TimeHIVToEnrol)
sameday$underTreatAll   <- factor(sameday$underTreatAll)
sameday$who   <- factor(sameday$who)
sameday$tb   <- factor(sameday$tb)
sameday$phone   <- factor(sameday$phone)
sameday$C.12   <- "uncensored"
sameday$C.12[is.na(sameday$outcome12m)] <- "censored"
sameday$C.12 <- as.factor(sameday$C.12)
sameday$C.18   <- "uncensored"
sameday$C.18[is.na(sameday$outcome18m)] <- "censored"
sameday$C.18 <- as.factor(sameday$C.18)
sameday$C.24   <- "uncensored"
sameday$C.24[is.na(sameday$outcome24m)] <- "censored"
sameday$C.24 <- as.factor(sameday$C.24)
sameday <- sameday[,c(colnames(sameday)[1:17],
                      "C.12","outcome12m","C.18","outcome18m","C.24","outcome24m","X_mi_m")]
# order

# generate 10 seperate (imputed) datasets -> using index from last column in sameday
M <- 10
for(m in 1:M){
  impname <- (paste("sameday",m,sep=""))
  assign(impname, sameday[sameday$X_mi_m==m,-24])
}



# for some IDs, var values change over the 10 obs -> bc imputed diff 10x using Multiple Imputation
sdd %>%
  group_by(id) %>%
  summarise(across(
    .cols = everything(), 
    .fns  = n_distinct, 
    .names = "{.col}_nunique")) %>% View()  # but have to subtract 1 from all (bc includes first of 11 entries with NA value!!)
# 12 vars: marital, edu, sex_pregnant, TimeHIVToEnrol, underTreatAll, cd4, who clinical stage, bmi, hb, alt, creat, phone
# (formerly: Marital, education, sex_pregnant, hivdiagenrol, preeeaa1, cd4studybase,
# whostudyase, cd4cat6, whocat3, bmi, bmicat3, hb, hbcat1, alt, altcat1, creat, creatcat2, phone, all VLelevation)

# e.g. obs 28: edu cat, marital, alt, creat change
sdd %>%
  filter(id == 28) %>%
  select(c("sexpregnancy","age_cat2","education_cat1","mstatus_cat1",
           "year","facility_cat2_rename","HIVdiagEnrol1","preEAAA1",
           "whostudybase","cd4studybase","tb_atARtinit_adapted",
           "bmi","hb","alt","creat","phone","sameDayEn14"))

summary(sameday1)
# underTreatAll = Indicator for Timing of HIV diag (0=before treat-all policy, 1=during/under policy)

table(sameday1$SAMEDAY)
# twice as many with same-day AR than early ART

# remove variables created from tMLE (censoring indicator & estimated outcome based on probs of unfavourable outcome)
art <- sameday1 %>% select(!c("C.12", "C.18", "C.24",
                             "outcome12m", "outcome18m", "outcome24m"))
summary(art)



# 2) Search for Positivity Violations ----

# 1. possibly viol for predictors of same-day ART (TimeHIVToEnrol, pregnant, PHC)

art %>% filter(TimeHIVToEnrol==3 & SAMEDAY==1) %>% nrow()/art %>% filter(TimeHIVToEnrol==3) %>% nrow()
art %>% filter(sex=="pregnant" & SAMEDAY==1) %>% nrow()/art %>% filter(sex=="pregnant") %>% nrow()
art %>% filter(facility %in% c("12-Dwaleni", "13-Gege" ,"14-Magubheleni", "15-Mahlandle",
                               "16-Mashobeni", "17-SOS", "18-Tfokotani", "19-Zombodze") & SAMEDAY==1) %>% nrow()/
art %>% filter(facility %in% c("12-Dwaleni", "13-Gege" ,"14-Magubheleni", "15-Mahlandle",
                                 "16-Mashobeni", "17-SOS", "18-Tfokotani", "19-Zombodze")) %>% nrow()
# indeed higher probs of same-day AR for these predictors, but not that extreme that could be considered as pos viol P(A=0|L)


# 2. Web Appendix 1 suggests CD4, unmarried, higher edu, year bc increased experience of medical personnel facilitate same-day ART
#    and treatment readiness & counselling which couldn't be measured tho (=> unmeasured conf)


# 3. additionally create table to check for extreme P(A) among binary conf
#    use customised fun defined in setup.R: 
binary_vars <- c("phone", "tb", "underTreatAll")
make_strata_table(art, A = "SAMEDAY", binary_vars = binary_vars) %>% arrange(proba_exp) %>% print(n = 125)
# among intersections of binary vars only, the group with viol is too small to be pos viol



## PoRT: continuous vars uncategorised ----
source("data/port_utils.R")
res <- list()
a_values <- c(0.01, 0.025, 0.05, 0.1)
gruber <- 5/(sqrt(nrow(art))*log(nrow(art)))
b_values <- c(0.01, gruber, 0.05, 0.1)
g_values <- 1:5

for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      res[[paste0("gamma = ", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A= "SAMEDAY", cov.quanti = c("cd4", "bmi", "hb", "alt", "creat"),
             cov.quali = c("sex", "age", "education", "marital", "year", "facility",
                           "TimeHIVToEnrol", "underTreatAll", "who", "tb", "phone"),  # 16 conf in total
             data = art, alpha = a, beta = b, gamma = g)
    }
  }
}
res
#sink(file = "output_port_uncat.txt")


# g = 1: largest viol subgroup: 16-Mashobeni (5.6%)
# a = 0.01: strata with cd4, hb, alt, bmi, 16-Mashobeni, but v small except 16-Mashobeni which until a <= 0.05, all with high P(A=1)

# g = 2 (most representative for clinical practices acc to Danelian et al):
#       largest sg: sex=pregnant & bmi< 30.2 & bmi>=25.1 (10.1%), 
#                   alt>=10.6 & facility=15-Mahlandle,16-Mashobeni (11%)
#                   creat< 96.2 & creat>=17 & facility=15-Mahlandle,16-Mashobeni (11.3%)
# a = 0.01: only one with P(A=1)=0: facility=11-Nhlangano-HC,19-Zombodze & creat< 64.7 & creat>=61.9 (1.2%)
#           -> shown only until b<=0.05 & not clear why bc with P(A=1)=0 still a viol under b=0.1
#           -> but generally, a=0.01 really small and 46 subgroups that maybe not relevant, esp not that still uncat
# a = 0.025/0.05/0.1 (for beta = 0.1): sex=pregnant & bmi< 30.2 & bmi>=25.1 -> why not alr for a=0.01?
#                                 (matches with sex as pred for A=1 in kersch paper, but always in combo bc pregnant alone doesn't have extreme P(A))
# a=0.01 & b=0.01/gruber/0.05 + a = 0.025, b=0.1: TimeHIVToEnrol=3 & other covar -> matches with TimeHIVToEnrol >= 90d as predictor!
# best choice a=0.05/0.1 & b=0.05/0.1: 6 subgroups in total

# g = 3,4: same largest sg as g=2, same thing with only facility=11-Nhlangano-HC,19-Zombodze & creat< 64.7 & creat>=61.9 (1.2%) having P(A)=0
# best choice a=0.05/0.1 & b=0.05/0.1: 8 subgroups in total

# g = 5: same largest sd as g=2, same thing with only facility=11-Nhlangano-HC,19-Zombodze & creat< 64.7 & creat>=61.9 (1.2%) having P(A)=0
# best choice a=0.05/0.1 & b=0.05/0.1: 12 subgroups in total

# g = 6-16 (no longer in output): same largest sd as g=2, same thing with only facility=11,19 & creat< 64.7 & creat>=61.9 (1.2%) having P(A)=0
# best choice a=0.05/0.1 & b=0.05/0.1: 12 subgroups in total, other subgroups mainly same (not more than 5 conf for def a stratum)


# essence: - a=0.01/0.025 (esp in combo with b=0.01) too detailed (recommend to only employ with knowledge/purpose; if care about small strata)
#          - with g=2, a=0.05, b=gruber/0.05/0.1 as rec in lit:
#            fac=15,16 & 142<= cd4 < 467 ; fac=16, sex=pregnant & 25.1<= bmi< 30.2
#          - with a=0.05/0.1:
#            15-, 16-, 17-, 18-, pregnant, bmi in 25-30, alt >= 10, creat, hb, age = 16-24/25-49, cd4, 2016, who=1, phone=Yes
#          - sometimes even if subgroup prop >a, only report for larger a not for smaller (not over all a consistently -> alr obs in simulations!)
#          - obv, for larger g, the subgroups are def as intersections of more conf
# -> assumption: are predictors for same-day ART often strata with large P(A=1)?
#    yes, but seldomly alone (only 16), else intersected with other covars bc only then P(A=1) is "extreme enough"
#    - for any a, sex=pregnant as viol in combo with other vars with high P(A=1)
#    - only for a < 0.05, TimeHIVToEnrol =3 (diag >=90d) & always in combo with other var
#    - only fac 15-18 reported as with extreme P(A=1), except 12 once for a=0.01,b=0.05
#    -> a=0.05/0.1 & b=0.05/0.1 from literature are sufficient for detecting the suspected viol groups
#       pregnant, PHC (IF PHC = 16-, but then again not all PHC?), but not TimeHIVToEnrol =3
# interesting that all viol are with v high P(A=1), i.e. v high P(same-day), except for the one with P(A=1)=0!

## PoRT: continuous vars categorised (acc to paper, can assume that domain experts decided on sensible thresholds there) ----
art_cat <- art
art_cat$cd4 <- cut(art$cd4, breaks = c(-Inf, 100, 200, 350, 500, Inf))
art_cat$bmi <- cut(art$bmi, breaks = c(0, 18.5, 25, Inf), right = F)
art_cat$hb <- cut(art$hb, breaks = c(0, 9.5, Inf))  # cat in middle of [9,10] to be exhaustive, but in paper not exhaustive (<=9, >=10)
art_cat$alt <- cut(art$alt, breaks = c(0, 42.5, Inf))  # same here
art_cat$creat <- cut(art$creat, breaks = c(0, 120, Inf))

res_cat <- list()
for (g in g_values) {
  for (a in a_values) {
    for (b in b_values) {
      res_cat[[paste0("gamma = ", g, ", alpha = ", a, ", beta = ", b)]] <-
        port(A= "SAMEDAY", cov.quanti = NULL,
             cov.quali = c("cd4", "bmi", "hb", "alt", "creat",  # all qualitative now
                           "sex", "age", "education", "marital", "year", "facility",
                           "TimeHIVToEnrol", "underTreatAll", "who", "tb", "phone"),  # 16 conf in total
             data = art_cat, alpha = a, beta = b, gamma = g)
    }
  } 
}
res_cat
# sink(file = "output_port_cat.txt")

# g = 1: only 16- reported for any a, b = 0.1, largest viol stratum = 16- (5.6%)

# g = 2: first time sex = pregnant included with "sex=pregnant & year=2016 ", largest viol stratum as g=1 (5.6%)
# a = 0.01: many combos involving 16-, year = 2016, pregnant
# a = 0.025: 16- in combo with year=2016, age = 25-49, cd4, edu
# a = 0.05: only facility=16-Mashobeni
# a = 0.1: no viol, i.e. strata with size >= 10% of sample have sufficient obs for both same-day and early ART

# g = 3: largest viol stratum = facility=15-,16-,17-,18- & year=2016 & alt=(0,42.5] (7.2%)

# g = 4-16: largest is creat=(0,120] & facility=15-,17-,18- & cd4=(100,200],(350,500],(500, Inf] & bmi=[25,Inf) (11.4%)

# essence: much fewer subgroups also within a=0.05/0.1 after categorisation (maximum is 12 now, vs 46 above)
#          - with g=2, a=0.05, b=gruber/0.05/0.1: only fac=16
#          - with a=0.05/0.1:
#            creat=(0,120] & facility=15-,17-,18- & cd4 > 100 & bmi=[25,Inf)
#            fac=16-, who=1 & year=2016 & cd4=(200,350],(350,500]
#            alt, year =2016
#          -> assumption: are predictors for same-day ART often strata with large P(A=1)?
#             yes, but seldomly alone (except for 16- which is a PHC), sex=pregnant & HIVTimeToEnrol only in combo with other covar
#             but TimeHIVToEnrol=3 only for a <0.05 and this time also TimeHIVToEnrol=2,3 (i.e. also diag >=1d)
#          - no longer any groups with P(A=1)=0
#          - for a=0.05/0.1 & b=0.05/0.1 as recommended in literature, detect:
#            largest subgroups for g=1/2,3,4
#            facility = 15,16,17,18 & year=2016 & bmi=[0,18.5),[18.5,25)
#            who=1 & year=2016 & cd4=(200,350],(350,500]
#            -> i.e. sufficient to detect assumed viol strata (predictors for same-day ART IF PHC = 16-) except for TimeHIVToEnrol
# overall trend as always: small a <=> smaller subgroups allowed <=> seems like more viol
# but prob not significant for estimation bc due to small size wouldn't have sign impact



## KBSD Rfast ----
source("kbsd.R")

# kbsd cannot deal with factors (only cont + binary, bc binary seen as cont),
# bc SD not computable for factors -> encode to cont ("multiclass") var with numerical repr:
# sex, age, edu, year, facility, timeHIVto Enrol, undertreatall, who, tb, phone

art_num <- art %>%
  mutate(sex = case_when(sex == "male" ~ 0, sex == "female (non-preg)" ~ 1, sex == "pregnant" ~ 2))
table(art$sex, art_num$sex)

art_num <- art_num %>%
  mutate(age = case_when(age == "16-24" ~ 0, age == "25-49" ~ 1, age == "50+" ~ 2))
table(art$age, art_num$age)

art_num <- art_num %>%
  mutate(education = case_when(education == "0" ~ 0, education == "1" ~ 1, education == "2" ~ 2, education == "3" ~3))
table(art$education, art_num$education)

art_num <- art_num %>%
  mutate(year = case_when(year == "2014" ~ 2014, year == "2015" ~ 2015, year == "2016" ~ 2016))
table(art$year, art_num$year)

art_num <- art_num %>%
  mutate(facility = case_when(facility == "11-Nhlangano-HC" ~ 11, facility == "12-Dwaleni" ~ 12,
                              facility == "13-Gege" ~ 13, facility == "14-Magubheleni" ~ 14,
                              facility == "15-Mahlandle" ~ 15, facility == "16-Mashobeni" ~ 16,
                              facility == "17-SOS" ~ 17, facility == "18-Tfokotani" ~ 18,
                              facility == "19-Zombodze" ~19))
table(art$facility, art_num$facility)

art_num <- art_num %>%
  mutate(TimeHIVToEnrol = case_when(TimeHIVToEnrol == "1" ~ 1, TimeHIVToEnrol == "2" ~ 2, TimeHIVToEnrol == "3" ~ 3))
table(art$TimeHIVToEnrol, art_num$TimeHIVToEnrol)

art_num <- art_num %>%
  mutate(underTreatAll = case_when(underTreatAll == "0" ~ 0, underTreatAll == "1" ~ 1))
table(art$underTreatAll, art_num$underTreatAll)

art_num <- art_num %>%
  mutate(who = case_when(who == "1" ~ 1, who == "2" ~ 2, who == "3" ~3, who == "4" ~4))
table(art$who, art_num$who)

art_num <- art_num %>%
  mutate(tb = case_when(tb == "0" ~ 0, tb == "1" ~ 1))
table(art$tb, art_num$tb)

art_num <- art_num %>%
  mutate(phone = case_when(phone == "No" ~ 0, phone == "Yes" ~ 1))
table(art$phone, art_num$phone)

summary(art_num)


## check for correlation ----
abs(cor(art_num)) > 0.3  # who & tb/cd4, sex & age/bmi/hb/creat, underTreatAll & TimeHIVToEnrol
cor(art_num$who, art_num$tb)
cor(art_num$who, art_num$cd4)  # the more advanced the stage, the fewer cd4 cells makes sense

cor(art_num$sex, art_num$age)
cor(art_num$sex, art_num$bmi)  # pregnancy associated with higher BMI makes sense
cor(art_num$sex, art_num$hb)
cor(art_num$sex, art_num$creat)

cor(art_num$underTreatAll, art_num$TimeHIVToEnrol) # strong neg corr: if initiated ART
# under treat all policy, then usually less time between diag&enrol bc policy facilitates faster enrolment & ART start?

cor(art_num$facility, art_num$SAMEDAY)  # also almost cor=0.3



# marital status, education, sex_pregnant, TimeHIVToEnrol, underTreatAll, cd4, who clinical stage, bmi, hb, alt, creat, phone

# all 16 confounders ---
# o1, o2 with intervened-on obs must only have the variables used lateron in disthalf_vec!!!
# imp to select all vars even if actually all cols, bc art_num df has weird indices
o1 <- art_num %>% mutate(SAMEDAY=1) %>%
  select(sex, facility, TimeHIVToEnrol, age, year, cd4, bmi, tb, phone, who, marital, alt, creat, education, underTreatAll, hb, SAMEDAY)
o2 <- art_num %>% mutate(SAMEDAY=0) %>% 
  select(sex, facility, TimeHIVToEnrol, age, year, cd4, bmi, tb, phone, who, marital, alt, creat, education, underTreatAll, hb, SAMEDAY)
res <- kbsd(data = art_num, int_data_list = list(o1, o2), type = "Rfast",
                 disthalf_vec = c(sex =sd(art_num$sex), facility =sd(art_num$facility), TimeHIVToEnrol=sd(art_num$TimeHIVToEnrol),
                                  age =sd(art_num$age), year =sd(art_num$year), cd4 = sd(art_num$cd4), bmi= sd(art_num$bmi), 
                                  tb =sd(art_num$tb), phone=sd(art_num$phone), who = sd(art_num$who),
                                  marital = sd(art_num$marital), alt =sd(art_num$alt), creat = sd(art_num$creat), 
                                  education = sd(art_num$education), underTreatAll = sd(art_num$underTreatAll),
                                  hb=sd(art_num$hb), SAMEDAY =0.5*sd(art_num$SAMEDAY)), plot.out = F)
res_plot <- kbsd(data = art_num, int_data_list = list(o1, o2), type = "Rfast",
     disthalf_vec = c(sex =sd(art_num$sex), facility =sd(art_num$facility), TimeHIVToEnrol=sd(art_num$TimeHIVToEnrol),
                      age =sd(art_num$age), year =sd(art_num$year), cd4 = sd(art_num$cd4), bmi= sd(art_num$bmi), 
                      tb =sd(art_num$tb), phone=sd(art_num$phone), who = sd(art_num$who),
                      marital = sd(art_num$marital), alt =sd(art_num$alt), creat = sd(art_num$creat), 
                      education = sd(art_num$education), underTreatAll = sd(art_num$underTreatAll),
                      hb=sd(art_num$hb), SAMEDAY =0.5*sd(art_num$SAMEDAY)))
#ggsave("art_bp.png", width = 5, height = 6)
# order in disthalf_vec has to be equal to col order in o1/o2?

# overall v few EDP bc 16 confounders <=> 16 dimensions
# fewer support for early ART than for same-day 
table(art$SAMEDAY)  # which makes sense bc almost twice as many same-day as early


# check who/what strata have few support among IV=1 (A=1) -> "problem": can just check univariately ---
# (within each var, we know what cat has few support then, but diff to agg to critical stratum then/to know which intersections have lowest support)
shift1 <- res[res$shift==1,]
outliers_ind1 <- shift1$diagnostic < quantile(shift1$diagnostic, 0.05)#0.05
outliers_val1 <- shift1[outliers_ind1, "diagnostic"]  # diag values of obs with few support
# check overall what the obs are that have few EDP:
for (i in names(art_num)[-17]) {
  print(paste0(i, ": ", mfv(art_num[outliers_ind1, i])))
}
# quantile=0.25: most from sex=1, from age=1, edu=2, marital=1, year=2015, fac=11, TimeHIVToEnrol=2, 
#  underTreatAll=1, who=1, cd4= [38,62], tb=0, bmi=22.7, hb=[11.7 12.8 13.8], alt=24, creat=27, phone=1
# quantile =0.05: only diff in sex=0, some clinical vars

# essence: obs that got sameday ART/A=1 most rarely of all cat (det bc few obs with similar covar vals were found), are
#   pred   non-pregnant women (men for quantile=0.05),
#   pred   from a SHC,
#   pred   with 1-89d between diag & enrolment,
#          25-49 yr olds (few in this medium-age group got same-day ART -> if only look at it univariately; combined with pregnancy/PHC prob more likely)
#          secondary edu, married, from year 2015, under treat all policy, first clinical stage, without tb, with phone
# but must always be careful: having more with low EDP from a certain cat if that cat is dominant anyway 
# is kind of normal? e.g. underTreatAll as almost 10x more with =1, so there'll prob more often be obs with =1 in general!
art_num[outliers_ind1,] %>% nrow()  # q=0.05: 67 obs
art_num[outliers_ind1,] %>% filter(facility %in% 11) %>% nrow()  # q=0.05: 43(majority) from SHC
art_num[outliers_ind1,] %>% select(TimeHIVToEnrol) %>% table()  # q=0.05: most from 1-89d 
art_num[outliers_ind1,] %>% select(tb) %>% table()

art_num[outliers_ind1,] %>% filter(facility == 11 & sex == 0)


# what obs/strata have few support in IV=2 (A=0) ---
shift2 <- res[res$shift==2,]
outliers_ind2 <- shift2$diagnostic < quantile(shift2$diagnostic, 0.05)  # indices of obs with few support
outliers_val2 <- shift2[outliers_ind2, "diagnostic"]  # diag values of obs with few support
# check overall what the obs are that have few EDP:
for (i in names(art_num)[-17]) {
  print(paste0(i, ": ", mfv(art_num[outliers_ind2, i])))
}
# quantile=0.25: most from sex=1, age=1, edu=2, marital=1, year=2015, fac=18, timeHIVToEnrol=1, underTreatAll=1, who=1,
#   cd4 = 74/189, tb=0, bmi=23.4, hb=11.1/12.6, alt=3, creat=27, phone=1
# quantile=0.05: diff for fac=16, timeHIVToEnrol=3, some clinical vars
# sex (pregnancy)
l_sex <- art_num[outliers_ind2, "sex"]
table(l_sex)   # expected higher count for pregnant (2), but makes sense that not few support for A=0 
# PHC
l_phc <- art_num[outliers_ind2, "facility"]
plot(l_phc, outliers_val2)
table(l_phc)  # most from 17, 18 (PHC) -> makes sense bc A=0 (early ART) rarer if in PHC, bc have poss to do sameday ART
# TimeHIVToEnrol
l_TimeHIVToEnrol <- art_num[outliers_ind2, "TimeHIVToEnrol"]
table(l_TimeHIVToEnrol)  # most had HIV diag & enrolment on same day, followed by >=90d in between
# if on same day, v unlikely to get early ART <=> few obs with HIV diag & enrolment on same day in A=0?
# but 2nd largest group had >=90d in between <=> few obs there makes sense bc if diag was longer time ago, more ready for A=1
art_num[outliers_ind2,] %>% nrow()  # q=0.05: 67 obs
art_num[outliers_ind2,] %>% filter(facility %in% 12:19) %>% nrow()  # q=0.05: 57(majority) from SHC
art_num[outliers_ind2,] %>% select(TimeHIVToEnrol) %>% table()  # q=0.05: most from <90d
art_num[outliers_ind2,] %>% select(sex) %>% table()

art_num[outliers_ind2,] %>% filter(facility != 11 & underTreatAll==1)


# essence: obs that got early ART/A=0 most rarely of all cat (det bc few obs with similar covar vals were found), are
#          non-pregnant women,
#          fac = 18/PHC (for quantile=0.05: fac=16/PHC)
#          with diag on same day as enrolment (q=0.05: >=90d in between),
#          25-49 yr olds, secondary edu, married, from year 2015,
#          under treat all policy, first clinical stage, without tb, with phone


## KBSD with categorised cont confounders ----

art_num$cd4 <- cut(art_num$cd4, breaks = c(-Inf, 100, 200, 350, 500, Inf))
art_num$bmi <- cut(art_num$bmi, breaks = c(0, 18.5, 25, Inf), right = F)
art_num$hb <- cut(art_num$hb, breaks = c(0, 9.5, Inf))  # cat in middle of [9,10] to be exhaustive, but in paper not exhaustive (<=9, >=10)
art_num$alt <- cut(art_num$alt, breaks = c(0, 42.5, Inf))  # same here
art_num$creat <- cut(art_num$creat, breaks = c(0, 120, Inf))
art_num

art_num <- art_num %>%
  mutate(cd4 = case_when(cd4 == "(-Inf,100]" ~ 1, cd4 == "(100,200]" ~ 2,
                         cd4 == "(200,350]" ~ 3, cd4 == "(350,500]" ~ 4,
                         cd4 == "(500, Inf]" ~ 5))
art_num <- art_num %>%
  mutate(bmi = case_when(bmi == "[0,18.5)" ~ 1, bmi == "[18.5,25)" ~ 2,
                         bmi == "[25,Inf)" ~ 3))
art_num <- art_num %>%
  mutate(hb = case_when(hb == "(0,9.5]" ~ 0, hb == "(9.5,Inf]" ~ 1))
art_num <- art_num %>%
  mutate(alt = case_when(alt == "(0,42.5]" ~ 0, alt == "(42.5,Inf]" ~ 1))
art_num <- art_num %>%
  mutate(creat = case_when(creat == "(0,120]" ~ 0, creat == "(120,Inf]" ~ 1))

any(is.na(art_num))

o1 <- art_num %>% mutate(SAMEDAY=1) %>%
  select(sex, facility, TimeHIVToEnrol, age, year, cd4, bmi, tb, phone, who, marital, alt, creat, education, underTreatAll, hb, SAMEDAY)
o2 <- art_num %>% mutate(SAMEDAY=0) %>% 
  select(sex, facility, TimeHIVToEnrol, age, year, cd4, bmi, tb, phone, who, marital, alt, creat, education, underTreatAll, hb, SAMEDAY)
res_hm_cat <- kbsd(data = art_num, int_data_list = list(o1, o2),
                   disthalf_vec = c(sex =sd(art_num$sex), facility =sd(art_num$facility), TimeHIVToEnrol=sd(art_num$TimeHIVToEnrol),
                                    age =sd(art_num$age), year =sd(art_num$year), cd4 = sd(art_num$cd4), bmi= sd(art_num$bmi), 
                                    tb =sd(art_num$tb), phone=sd(art_num$phone), who = sd(art_num$who),
                                    marital = sd(art_num$marital), alt =sd(art_num$alt), creat = sd(art_num$creat), 
                                    education = sd(art_num$education), underTreatAll = sd(art_num$underTreatAll),
                                    hb=sd(art_num$hb), SAMEDAY =0.5*sd(art_num$SAMEDAY)), plot.out = F)
#saveRDS(res_hm_cat, "output_kbsd_cat.RDS")
res_hm_plot <- kbsd(data = art_num, int_data_list = list(o1, o2),
                    disthalf_vec = c(sex =sd(art_num$sex), facility =sd(art_num$facility), TimeHIVToEnrol=sd(art_num$TimeHIVToEnrol),
                                     age =sd(art_num$age), year =sd(art_num$year), cd4 = sd(art_num$cd4), bmi= sd(art_num$bmi), 
                                     tb =sd(art_num$tb), phone=sd(art_num$phone), who = sd(art_num$who),
                                     marital = sd(art_num$marital), alt =sd(art_num$alt), creat = sd(art_num$creat), 
                                     education = sd(art_num$education), underTreatAll = sd(art_num$underTreatAll),
                                     hb=sd(art_num$hb), SAMEDAY =0.5*sd(art_num$SAMEDAY)))
#ggsave("art_bp_cat.png", width = 5, height = 6)



## KBSD harmonic mean ----

o1 <- art_num %>% mutate(SAMEDAY=1) %>%
  select(sex, facility, TimeHIVToEnrol, age, year, cd4, bmi, tb, phone, who, marital, alt, creat, education, underTreatAll, hb, SAMEDAY)
o2 <- art_num %>% mutate(SAMEDAY=0) %>% 
  select(sex, facility, TimeHIVToEnrol, age, year, cd4, bmi, tb, phone, who, marital, alt, creat, education, underTreatAll, hb, SAMEDAY)
res_hm <- kbsd(data = art_num, int_data_list = list(o1, o2), type = "harmonicmean",
            disthalf_vec = c(sex =sd(art_num$sex), facility =sd(art_num$facility), TimeHIVToEnrol=sd(art_num$TimeHIVToEnrol),
                             age =sd(art_num$age), year =sd(art_num$year), cd4 = sd(art_num$cd4), bmi= sd(art_num$bmi), 
                             tb =sd(art_num$tb), phone=sd(art_num$phone), who = sd(art_num$who),
                             marital = sd(art_num$marital), alt =sd(art_num$alt), creat = sd(art_num$creat), 
                             education = sd(art_num$education), underTreatAll = sd(art_num$underTreatAll),
                             hb=sd(art_num$hb), SAMEDAY =0.5*sd(art_num$SAMEDAY)), plot.out = F)
#saveRDS(res_hm, "output_kbsd_hm.RDS)
res_hm_plot <- kbsd(data = art_num, int_data_list = list(o1, o2), type = "harmonicmean",
     disthalf_vec = c(sex =sd(art_num$sex), facility =sd(art_num$facility), TimeHIVToEnrol=sd(art_num$TimeHIVToEnrol),
                      age =sd(art_num$age), year =sd(art_num$year), cd4 = sd(art_num$cd4), bmi= sd(art_num$bmi), 
                      tb =sd(art_num$tb), phone=sd(art_num$phone), who = sd(art_num$who),
                      marital = sd(art_num$marital), alt =sd(art_num$alt), creat = sd(art_num$creat), 
                      education = sd(art_num$education), underTreatAll = sd(art_num$underTreatAll),
                      hb=sd(art_num$hb), SAMEDAY =0.5*sd(art_num$SAMEDAY)))
#ggsave("art_bp_hm_no_fac16.png", width = 5, height = 6)

# double-check what strata have few support among IV=1 (A=1) -> "problem": can just check univariately ---
shift1 <- res_hm[res_hm$shift==1,]
outliers_ind1 <- shift1$diagnostic < quantile(shift1$diagnostic, 0.05)  # indices of obs with few support
outliers_val1 <- shift1[shift1$diagnostic < quantile(shift1$diagnostic, 0.05), "diagnostic"]  # diag values of obs with few support
# check overall what the obs are that have few EDP:
for (i in names(art_num)[-17]) {
  print(paste0(i, ": ", mfv(art_num[outliers_ind1, i])))
}
# diff values for some vars now: TimeHIVToEnrol=1, cd4=74, bmi=23 instead of 22.7, hb=10.7, alt=3&25, creat=54
# few support for diag & enrolment on same day makes sense -> would rather not start ART immediately (those with diag longer time ago would rather start fast after enrolment)


# double-check what strata have few support among IV=2 (A=0) ---
shift2 <- res_hm[res_hm$shift==2,]
outliers_ind2 <- shift2$diagnostic < quantile(shift2$diagnostic, 0.05)  # indices of obs with few support
outliers_val2 <- shift2[shift2$diagnostic < quantile(shift2$diagnostic, 0.05), "diagnostic"]  # diag values of obs with few support
# check overall what the obs are that have few EDP:
for (i in names(art_num)[-17]) {
  print(paste0(i, ": ", mfv(art_num[outliers_ind2, i])))
}
# diff vals for: fac=11, cd4=74, tb=0, bmi=22.5&23, hb=10.7, alt=25, creat= 17&54
# few support for fac=11 among A=0 does not make sense.. but prob just bc many from fac=11
# so some or most have lots of support for A=0, but some also few support bc with other covar vals!

r <- (res_plot + ylab("Effective Data Points (EDP)") + ylim(0,10) )| 
  (res_hm_plot + ylab("Effective Data Points (EDP)") + ylim(0,250))
#ggsave("output_kbsd_both.png")



## KBSD hm with categorised continuous confounders (cd4, bmi, hb, alt, creat) ----

o1 <- art_num %>% mutate(SAMEDAY=1) %>%
  select(sex, facility, TimeHIVToEnrol, age, year, cd4, bmi, tb, phone, who, marital, alt, creat, education, underTreatAll, hb, SAMEDAY)
o2 <- art_num %>% mutate(SAMEDAY=0) %>% 
  select(sex, facility, TimeHIVToEnrol, age, year, cd4, bmi, tb, phone, who, marital, alt, creat, education, underTreatAll, hb, SAMEDAY)
res_hm_cat <- kbsd(data = art_num, int_data_list = list(o1, o2), type = "harmonicmean",
               disthalf_vec = c(sex =sd(art_num$sex), facility =sd(art_num$facility), TimeHIVToEnrol=sd(art_num$TimeHIVToEnrol),
                                age =sd(art_num$age), year =sd(art_num$year), cd4 = sd(art_num$cd4), bmi= sd(art_num$bmi), 
                                tb =sd(art_num$tb), phone=sd(art_num$phone), who = sd(art_num$who),
                                marital = sd(art_num$marital), alt =sd(art_num$alt), creat = sd(art_num$creat), 
                                education = sd(art_num$education), underTreatAll = sd(art_num$underTreatAll),
                                hb=sd(art_num$hb), SAMEDAY =0.5*sd(art_num$SAMEDAY)), plot.out = F)
#saveRDS(res_hm_cat, "output_kbsd_hm_cat.RDS")
res_hm_plot <- kbsd(data = art_num, int_data_list = list(o1, o2), type = "harmonicmean",
                    disthalf_vec = c(sex =sd(art_num$sex), facility =sd(art_num$facility), TimeHIVToEnrol=sd(art_num$TimeHIVToEnrol),
                                     age =sd(art_num$age), year =sd(art_num$year), cd4 = sd(art_num$cd4), bmi= sd(art_num$bmi), 
                                     tb =sd(art_num$tb), phone=sd(art_num$phone), who = sd(art_num$who),
                                     marital = sd(art_num$marital), alt =sd(art_num$alt), creat = sd(art_num$creat), 
                                     education = sd(art_num$education), underTreatAll = sd(art_num$underTreatAll),
                                     hb=sd(art_num$hb), SAMEDAY =0.5*sd(art_num$SAMEDAY)))
#ggsave("art_bp_hm_cat.png", width = 5, height = 6)
