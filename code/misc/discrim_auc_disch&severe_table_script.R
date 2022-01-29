##  AUC for two competing outcomes
## given data on day 0, day 2, 4, and 8 
rm(list=ls())
setwd("S:/COVID_LDA/code/zwang_mgbowring/prediction_model")
source("00_loadlibs.R")
dt_joint <- readRDS("data/data_for_joint.rds")

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

dt.pre <- dt.pre  %>% mutate(sex = factor(sex, levels = c("Female", "Male")),
                             smoking = factor(smoking, levels = c("Never Smoker", "Former Smoker", "Current Smoker")),
                             charlson_cat = factor(charlson_cat, levels = c("0", "1-2" , "3-4", ">5")),
                             demo5 = factor(demo5, levels =   c("Young_White", "Old_White" , "Young_Black", "Old_Black", "Latinx/Other")),
                             bmi_cat = factor(bmi_cat, levels = c("<30", ">=30")))

load("S:/COVID_LDA/code/zwang_mgbowring/data/vitals_data0809.RData")
tmp  <- vitals_merge %>% dplyr::select(osler_id, interval, group2, first_vented,last_interval, outcome2) %>% unique()
dt.pre.mult.tmp <- dt.joint %>% 
  left_join(tmp) %>% 
  mutate(first_vented_interval = ifelse(first_vented == 'Inf', NA, first_vented),
         endpoint = ifelse(status == 3, first_vented_interval,
                           ifelse(status != 3, last_interval, NA))) %>% 
  filter(group2 == 'Not vented' | (!is.na(first_vented_interval)& first_vented_interval!=0 & interval <= first_vented_interval)) %>% 
  mutate(Wt = ifelse(interval < endpoint, 0, status),
         Wt = ifelse(!is.na(first_vented_interval)&first_vented_interval == interval, 3, Wt))
dt.endpoints <- dt.pre.mult.tmp %>% 
  dplyr::group_by(osler_id) %>% 
  dplyr::summarise(endpoint.u = min(endpoint), status.u = min(status))%>% 
  dplyr::select(osler_id, endpoint.u,status.u) %>% 
  unique() %>% 
  filter(endpoint.u != 0)
dt.pre.mult <- dt.pre.mult.tmp %>% dplyr::filter(!(osler_id %in% c(id.check.pre$osler_id))) 
#length(unique(dt.pre.mult$osler_id))
id.c<-unique(dt.pre$osler_id)
id.cm<-unique(dt.pre.mult$osler_id)
length(id.c[!(id.c %in% id.cm)])

## every patient needs to have an event
## fix dt1
dt1<-dt.pre.mult
sum(dt1$Wt != 0)
length(unique(dt1$osler_id[dt1$Wt !=0]))
length(unique(dt1$osler_id))
miso_id <- unique(dt1$osler_id)[!(unique(dt1$osler_id) %in% unique(dt1$osler_id[dt1$Wt !=0]))] #2 patients don't have outcome
dt1c <- dt1
for(i in 1:length(miso_id)){
  tmp0 <- dt1 %>% filter(osler_id == miso_id[i])
  new <- tmp0[1,]
  new$interval <- max(tmp0$interval)+1
  new$rev_day <- 0
  new <- new %>% 
    mutate(Wt = ifelse(is.na(first_vented_interval) & group2=="Not vented" & outcome2 == "Discharged or AC", 1,
                       ifelse(is.na(first_vented_interval) & group2=="Not vented" & outcome2 == "Died", 2,
                              ifelse(!is.na(first_vented_interval), 3, NA))))
  dt1c <- rbind(new, dt1c)
}
sum(dt1c$Wt != 0)
length(unique(dt1c$osler_id[dt1c$Wt !=0]))
length(unique(dt1c$osler_id))

dt.pre.mult <- dt1c



idlist1 <- unique(dt.pre$osler_id[dt.pre$status == 1])
idlist2 <- unique(dt.pre$osler_id[dt.pre$status == 2])
idlist3 <- unique(dt.pre$osler_id[dt.pre$status == 3])

set.seed(2022) 
test_inx1 = createFolds(idlist1, k=5, list = T) 
test_inx2 = createFolds(idlist2, k=5, list = T) 
test_inx3 = createFolds(idlist3, k=5, list = T) 

test_ids1 <-test_ids2 <-test_ids3 <- list()
for(i in 1:length(test_inx1)){
  test.inx <- test_inx1[[i]]
  test_ids1[[i]] = idlist1[test.inx]
}
for(i in 1:length(test_inx2)){
  test.inx <- test_inx2[[i]]
  test_ids2[[i]] = idlist2[test.inx]
}
for(i in 1:length(test_inx1)){
  test.inx <- test_inx3[[i]]
  test_ids3[[i]] = idlist3[test.inx]
}



path<-"paperPrep/finalizedCode/figs/calibration/"
load(paste0(path,"id_remove_2.RData")); id_remove2 <- id.remove
load(paste0(path,"id_remove_4.RData")); id_remove4 <- id.remove
load(paste0(path,"id_remove_8.RData")); id_remove8 <- id.remove
id_remove2  <- as.data.frame(id_remove2);id_remove4  <- as.data.frame(id_remove4);id_remove8  <- as.data.frame(id_remove8)
colnames(id_remove2) <- colnames(id_remove4) <- colnames(id_remove8) <- c("cv", "position", "osler_id")

i <- 2 #given data until day 2
## filter out patients have an event prior to day 4, keep only event occurred after day 4.
dt.endpoints <- dt.pre.mult.tmp %>% 
  dplyr::group_by(osler_id) %>% 
  dplyr::summarise(endpoint.u = min(endpoint), status.u = min(status))%>% 
  dplyr::select(osler_id, endpoint.u, status.u) %>% 
  unique() %>% 
  filter(endpoint.u > i)

#load("paperPrep/finalizedCode/figs/calibration/bkwr_di=2_M.RData")
#load("paperPrep/finalizedCode/figs/calibration/forw_di=2_M.RData")
 load("paperPrep/finalizedCode/figs/calibration/bkwr_di=0.RData")
 load("paperPrep/finalizedCode/figs/calibration/forw_di=0.RData")


f.prob.save <- forw.save2[[5]]
b.prob.save <- bkwr.save2[[5]]
id_remove <- id_remove2

b.data.roc <- f.data.roc <- NULL
for(v in 1:5){
  message(v)
  test_ids <- c(test_ids1[[v]],test_ids2[[v]],test_ids3[[v]])
  test_ids_4d <- test_ids[test_ids %in% dt.endpoints$osler_id]
  ## filter out patients missing all biomarkers at every day 
  testids_nonmis <- dt.pre %>%  ## keep patients with at least one marker oberved at any day
    filter(osler_id %in% test_ids_4d) %>%
    filter(complete.cases(temp_c_m, sao2_fio2_ratio_m, pulse_m)) %>% dplyr::select(osler_id) %>% unique()
  test_ids_4d <- test_ids_4d[test_ids_4d %in% testids_nonmis$osler_id]
  
  b.prob.save.cv <- f.prob.save.cv <- list()
  b.prob.save.cv <- b.prob.save
  f.prob.save.cv <- f.prob.save

  tmp.data.b <- NULL
  for(s in 1:length(test_ids_4d)){
    sid <- test_ids_4d[s]
    if(sid %in% id_remove$osler_id){
      tmp.data.b <- rbind(tmp.data.b, NULL)
    }else{
      data <- as.matrix(b.prob.save.cv[[v]][[s]])
      colnames(data) <- c('E1', 'E2', 'E3')
      if(i == 0){data <- cbind(data, day = 1:20)}else{data <- cbind(data, day = (i+1):20)}
      data[is.na(data)] <- 0
      data <- as.data.frame(data) %>% mutate(cumE1 = cumsum(E1),
                                             cumE2 = cumsum(E2), 
                                             cumE3 = cumsum(E3))
      tmp <- dt1c %>% filter(osler_id == sid)
      endpoint <- min(unique(tmp$endpoint));W <- max(tmp$Wt)
      data <- data %>% mutate(O1 = ifelse(day==endpoint & W==1, 1, 0),
                              O2 = ifelse(day==endpoint & W==2, 1, 0),
                              O3 = ifelse(day==endpoint & W==3, 1, 0))
      data <- data %>% mutate(cumO1 = cumsum(O1),
                              cumO2 = cumsum(O2),
                              cumO3 = cumsum(O3))
      data$osler_id <- sid
      data$cv <- v
      tmp.data.b <- rbind(tmp.data.b, data)
    }
  }
  b.data.roc <- rbind(b.data.roc, tmp.data.b)


  tmp.data.f <- NULL
  for(s in 1:length(test_ids_4d)){
    sid <- test_ids_4d[s]
    if(sid %in% id_remove$osler_id){
      tmp.data.f <- rbind(tmp.data.f, NULL)
    }else{
    data <- as.matrix(f.prob.save.cv[[v]][[s]])
    colnames(data) <- c('E1', 'E2', 'E3')
    data <- cbind(data, day = 0:20)
    data[is.na(data)] <- 0
    data <- as.data.frame(data) %>% mutate(cumE1 = cumsum(E1),
                                           cumE2 = cumsum(E2), 
                                           cumE3 = cumsum(E3))
    tmp <- dt1c %>% filter(osler_id == sid)
    endpoint <- min(unique(tmp$endpoint));W <- max(tmp$Wt)
    data <- data %>% mutate(O1 = ifelse(day==endpoint & W==1, 1, 0),
                            O2 = ifelse(day==endpoint & W==2, 1, 0),
                            O3 = ifelse(day==endpoint & W==3, 1, 0))
    data <- data %>% mutate(cumO1 = cumsum(O1),
                            cumO2 = cumsum(O2),
                            cumO3 = cumsum(O3))
    data$osler_id <- sid
    data$cv <- v
    tmp.data.f <- rbind(tmp.data.f, data)
    }
  }
  f.data.roc <- rbind(f.data.roc, tmp.data.f)
}

f.data.roc <- f.data.roc %>%
  mutate(O1_new = O1, O2_new = O2 + O3,
         E1_new = E1, E2_new = E2 + E3,
         cumO1_new = cumO1, cumO2_new = cumO2 + cumO3,
         cumE1_new = cumE1, cumE2_new = cumE2 + cumE3)
b.data.roc <- b.data.roc %>%
  mutate(O1_new = O1, O2_new = O2 + O3,
         E1_new = E1, E2_new = E2 + E3,
         cumO1_new = cumO1, cumO2_new = cumO2 + cumO3,
         cumE1_new = cumE1, cumE2_new = cumE2 + cumE3)
library(pROC)
b_auc_uet<- f_auc_uet <- NULL
b_auc_uet<- f_auc_uet <- NULL
for(t in (i+1):20){
  # backward
  b.tmp.auc <- NULL
  datasurv <- b.data.roc %>% filter(day == t)
  
  auc.e1 <- auc(datasurv$cumO1_new,datasurv$cumE1_new)
  auc.e2 <- auc(datasurv$cumO2_new,datasurv$cumE2_new)
  
  
  b.tmp.auc <- data.frame(event = c("disch", "severe"), 
                          day = rep(t, 2), 
                          mean = c(auc.e1, auc.e2))

  b_auc_uet <- rbind(b_auc_uet, b.tmp.auc)
  
  # forward
  f.tmp.auc <- NULL
  datasurv <- f.data.roc %>% filter(day == t)
  
  auc.e1 <- auc(datasurv$cumO1_new,datasurv$cumE1_new)
  auc.e2 <- auc(datasurv$cumO2_new,datasurv$cumE2_new)
  
  f.tmp.auc <- data.frame(event = c("disch", "severe"), 
                          day = rep(t, 2), 
                          mean = c(auc.e1, auc.e2))
  f_auc_uet <- rbind(f_auc_uet, f.tmp.auc)
}

pauc<- ggplot()+
  geom_line(data = b_auc_uet,aes(x= day, y = mean, group = event, color = event, linetype = "backward"))+
  #geom_errorbar(data = b_auc_uet %>% filter(day %in% c(4,9,14,19)), 
  # aes(x= day, ymin = lwr, ymax = upr,group = event, color = event , linetype = "backward"), alpha = 0.6, width = 0.4)+
  geom_line(data = f_auc_uet,aes(x= day, y = mean, group = event, color = event, linetype = "forward"))+
  #geom_errorbar(data = f_auc_uet %>% filter(day %in% c(5,10,15,20)), 
  # aes(x= day, ymin = lwr, ymax = upr,group = event,color = event,  linetype = "forward"), alpha = 0.6, width = 0.4)+
  theme_classic()+
  ylim(0.5,1)+
  ylab("AUC")+
  xlab("Days projected")+
  scale_color_manual(breaks = c("disch", "severe"),
                     label = c("Discharge", "Severe"),
                     values = c("#2FB3CA","#F1564F"),name=" ")+
  scale_linetype(name = "")+
  ggtitle(paste0("AUC(u=",i,", e, t)"))+
  xlim(0,20)

f_auc_uet <- f_auc_uet %>% mutate(model = "forw", U = i)
b_auc_uet <- b_auc_uet %>% mutate(model = "bkwr", U = i)

auc_uet_8  =  rbind(f_auc_uet, b_auc_uet)
auc_uet_4  =  rbind(f_auc_uet, b_auc_uet)
auc_uet_2  =  rbind(f_auc_uet, b_auc_uet)
auc_uet_0  =  rbind(f_auc_uet, b_auc_uet)

auc.stats <- rbind(auc_uet_0, auc_uet_2, auc_uet_4, auc_uet_8)
auc.stats.d20 = auc.stats %>% filter(day == 20)
discrim_result <- auc.stats.d20  %>% 
  mutate(event_model = paste0(event, "_", model)) %>% 
  dplyr::select(-event, -model) %>% 
  spread(event_model, mean)
save(discrim_result, file = "paperPrep/finalizedCode/figs/calibration/discrim_result_table.RData")

## paste discrim & calib table together ----
load("paperPrep/finalizedCode/figs/calibration/calib_result_table.RData")
load("paperPrep/finalizedCode/figs/calibration/discrim_result_table.RData")

out <- cbind(calib_result, discrim_result)
out <- out[,-1]

round(out,2)
