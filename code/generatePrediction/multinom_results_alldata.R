rm(list=ls())
nbiom <- 3
df.re <- 3
nb = 2000; nt = 8000
source("code/00_loadlibs.R")
setwd(wd)

## load data 
dt1c <- data
## backward ----
dt.pre.mult <- dt1c

dt.pre.mult.train <- dt1c %>% filter((osler_id %in% unique(dt1c$osler_id))) #length(unique(dt.pre.mult.train$osler_id))
dt.pre.mult.train <- dt.pre.mult.train %>% filter(!is.na(endpoint)) %>% filter(endpoint != 0)
#dt.pre.mult.train <- dt.pre.mult.train %>% filter(interval != 0)


dt.pre.mult.train <- dt.pre.mult.train  %>% mutate(sex = factor(sex, levels = c("Female", "Male")),
                                                   smoking = factor(smoking, levels = c("Never Smoker", "Former Smoker", "Current Smoker")),
                                                   charlson_cat = factor(charlson_cat, levels = c("0", "1-2" , "3-4", ">5")),
                                                   demo5 = factor(demo5, levels =   c("Young_White", "Old_White" , "Young_Black", "Old_Black", "Latinx/Other")),
                                                   bmi_cat = factor(bmi_cat, levels = c("<30", ">=30")))
mult.pre <- multinom(factor(Wt) ~ resp_bl + temp_bl + pulse_bl + sao2fio2_bl +
                       crp_bl + alc_bl + dd_bl + gfr_bl +
                       charlson_cat + demo5 + sex + bmi_cat +
                       smoking + ns(interval, knots = c(2,4, 6 ,8,12,16)),
                     data = dt.pre.mult.train, maxit = 1000)


save(mult.pre, file = "model/alldata_multinom_model.RData")

load("model/alldata_multinom_model.RData")
t(summary(mult.pre)$coefficients)
t(summary(mult.pre)$standard.errors)


## forward------
dt.pre.mult.train <- dt.pre.mult.train %>%
  group_by(osler_id) %>% 
  arrange(interval) %>% 
  mutate(prev1_sao2 = lag(sao2_norm),
         prev1_temp = lag(temp_norm),
         prev1_pulse = lag(pulse_norm))%>% 
  mutate(prev2_sao2 = lag(sao2_norm, 2),
         prev2_temp = lag(temp_norm, 2),
         prev2_pulse = lag(pulse_norm, 2))
dt.pre.mult.train <- dt.pre.mult.train %>% 
  group_by(osler_id) %>% 
  mutate(prev1_interval = lag(interval),
         prev2_interval = lag(interval, 2)) %>% 
  mutate(slope_sao2 = (prev1_sao2-prev2_sao2)/(prev1_interval-prev2_interval),
         slope_temp = (prev1_temp-prev2_temp)/(prev1_interval-prev2_interval),
         slope_pulse = (prev1_pulse-prev2_pulse)/(prev1_interval-prev2_interval))
dt.pre.mult.train <- dt.pre.mult.train  %>% mutate(sex = factor(sex, levels = c("Female", "Male")),
                                                   smoking = factor(smoking, levels = c("Never Smoker", "Former Smoker", "Current Smoker")),
                                                   charlson_cat = factor(charlson_cat, levels = c("0", "1-2" , "3-4", ">5")),
                                                   demo5 = factor(demo5, levels =   c("Young_White", "Old_White" , "Young_Black", "Old_Black", "Latinx/Other")),
                                                   bmi_cat = factor(bmi_cat, levels = c("<30", ">=30")))

mult.pre.forw <- multinom(factor(Wt) ~ resp_bl + temp_bl + pulse_bl + sao2fio2_bl +
                                  crp_bl + alc_bl + dd_bl + gfr_bl +
                                  charlson_cat + demo5 + sex + bmi_cat +
                                  smoking + ns(interval, knots = c(4, 6 ,8,12,16)) +
                                  prev1_sao2 + prev1_temp + prev1_pulse +
                                  slope_sao2 + slope_temp + slope_pulse,
                                na.action = na.omit,
                                data = dt.pre.mult.train, maxit = 1000)
save(mult.pre.forw, file = "model/alldata_multinom_model_forw.RData")

## paste to excel --- ?
load("model/alldata_multinom_model.RData")
t(summary(mult.pre)$coefficients)
t(summary(mult.pre)$standard.errors)

load("model/alldata_multinom_model_forw.RData")
t(summary(mult.pre.forw)$coefficients)
t(summary(mult.pre.forw)$standard.errors)
