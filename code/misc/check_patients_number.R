rm(list=ls())
setwd("S:/COVID_LDA/code/zwang_mgbowring/prediction_model")
source("00_loadlibs.R")
dt_joint <- readRDS("data/data_for_joint.rds")
vitals <- readRDS("S:/COVID_LDA/code/zwang_mgbowring/prediction_model/data/vitals.rds")
# length(unique(vitals$osler_id))
# length(unique(dt_joint$osler_id))
# 
# tmp1 <- dt_joint %>%  filter(complete.cases(sao2_fio2_ratio_m, temp_c_m, pulse_m))
# tmp2 <- vitals %>%  filter(complete.cases(sao2_fio2_ratio_m, temp_c_m, pulse_m))
# length(unique(tmp1$osler_id))
# length(unique(tmp2$osler_id))
# 
# dt.pre <- dt_joint %>% filter(status %in% c(1,2,3)) #length(unique(dt.pre$osler_id)) #1762 patients
# length(unique(dt.pre$osler_id))
# dt.pre$long <- ifelse(dt.pre$rev_day < -20,'long','short')
# id.check.pre <- dt.pre %>% dplyr::select(osler_id, long, status) %>% filter(long == "long") %>% unique()#28 patients
# dt.pre <- dt.pre %>% dplyr::filter(!(osler_id %in% id.check.pre$osler_id)) ## 1734 patients
# length(unique(dt.pre$osler_id))


imp <- readRDS("paperPrep/data/processed_data/baseline_char_imputed.rds")
dtimp <- complete(imp, action=1)

dt_joint <- dt_joint %>%
  dplyr::select(osler_id, interval, rev_day, status, sao2_fio2_ratio_m, pulse_m, temp_c_m, outcome)%>% 
  left_join(dtimp)
## create X categories for models
dt_joint <- dt_joint %>% mutate(charlson_cat = ifelse(charlson == 0, "0",
                                                      ifelse(charlson %in% c(1,2), "1-2", 
                                                             ifelse(charlson %in% c(3,4),"3-4",
                                                                    ifelse(charlson >= 5, ">5", NA)))),
                                race_new = ifelse(race == "Black", race, 
                                                  ifelse(race == "White or Caucasian", "White",
                                                         ifelse(!(race %in% c("White or Caucasian", "Black")), "Latinx/Other", NA))),
                                bmi_cat = ifelse(bmi < 30, "<30", ifelse(bmi >=30, ">=30", NA)), 
                                age_cat = ifelse(age <= 74, '<75',
                                                 ifelse(age > 74, '>=75', NA)),
                                demo5 = ifelse(race_new == "Latinx/Other", "Latinx/Other",
                                               ifelse(race_new == "Black" & age_cat  == "<75", "Young_Black",
                                                      ifelse(race_new == "Black" & age_cat == ">=75", "Old_Black",
                                                             ifelse(race_new == "White" & age_cat == "<75", "Young_White",
                                                                    ifelse(race_new == "White"&age_cat == ">=75", "Old_White", NA))))))
dt.joint <- dt_joint %>% 
  filter(complete.cases(dd_bl, crp_bl, alc_bl, gfr_bl, 
                        resp_bl, sao2fio2_bl, pulse_bl, temp_bl, charlson_cat, demo5, 
                        sex, bmi_cat, smoking, rev_day, status))
dt.joint <- dt.joint %>%  filter(complete.cases(sao2_fio2_ratio_m, temp_c_m, pulse_m)) 
length(unique(dt.joint$osler_id))
## filter for pre-ventilation
dt.pre <- dt.joint %>% filter(status %in% c(1,2,3)) #length(unique(dt.pre$osler_id)) #1762 patients
dt.pre$long <- ifelse(dt.pre$rev_day < -20,'long','short')
id.check.pre <- dt.pre %>% dplyr::select(osler_id, long, status) %>% filter(long == "long") %>% unique()#28 patients
dt.pre <- dt.pre %>% dplyr::filter(!(osler_id %in% id.check.pre$osler_id)) ## 1734 patients

dt.pre <- dt.pre %>% 
  dplyr::rename(sao2_fio2_ratio_o = sao2_fio2_ratio_m,
                temp_c_o = temp_c_m,
                pulse_o = pulse_m) %>% 
  mutate(sao2_fio2_ratio_m = sao2_fio2_ratio_o/100,
         temp_c_m = temp_c_o/10,
         pulse_m =  pulse_o/10)

library(bestNormalize)
dt.pre <- dt.pre %>% 
  filter(complete.cases(sao2_fio2_ratio_o, temp_c_o, pulse_o)) %>% 
  mutate(sao2_norm = orderNorm(sao2_fio2_ratio_o)$x.t,
         temp_norm = orderNorm(temp_c_o)$x.t,
         pulse_norm = orderNorm(pulse_o)$x.t)
length(unique(dt.pre$osler_id))

tmp <- dt.pre %>% dplyr::select(osler_id, status) %>% unique()
table(tmp$status)

status3 <- dt.pre %>% filter(osler_id %in% tmp$osler_id[tmp$status==3])
tmp3 <- status3 %>% dplyr::select(osler_id, outcome) %>% unique()
table(tmp3$outcome)
