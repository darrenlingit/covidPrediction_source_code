rm(list=ls())
source("00_loadlibs.R")
setwd(wd)
dt_joint <- readRDS("data/data_for_joint.rds")
imp <- readRDS("paperPrep/data/processed_data/baseline_char_imputed.rds")
dtimp <- complete(imp, action=1)
dt_joint <- dt_joint %>%
  dplyr::select(osler_id, interval, rev_day, status, sao2_fio2_ratio_m, pulse_m, temp_c_m, outcome)%>% 
  left_join(dtimp)
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

path<-"paperPrep/finalizedCode/figs/calibration/"
load(paste0(path, "forw_di=", 2 ,".RData"))
load(paste0(path, "forw_di=", 4 ,".RData"))
load(paste0(path, "forw_di=", 8 ,".RData"))
load(paste0(path,"id_remove_2.RData")); id_remove2 <- id.remove
load(paste0(path,"id_remove_4.RData")); id_remove4 <- id.remove
load(paste0(path,"id_remove_8.RData")); id_remove8 <- id.remove
load("paperPrep/finalizedCode/figs/calibration/cv_ids.RData")
test_ids1 = cv_ids[[1]]
test_ids2 = cv_ids[[2]]
test_ids3 = cv_ids[[3]]

id_remove2  <- as.data.frame(id_remove2);id_remove4  <- as.data.frame(id_remove4);id_remove8  <- as.data.frame(id_remove8)
colnames(id_remove2) <- colnames(id_remove4) <- colnames(id_remove8) <- c("cv", "position", "osler_id")

## change to correct data
i <- 2
exp.freq.cv <- forw.save2[[3]]
obs.freq.cv <- forw.save2[[4]]
pat_prob_save <- forw.save2[[5]]
id_remove <- id_remove2

## filter out patients have an event prior to day 4, keep only event occurred after day 4.
dt.endpoints <- dt.pre.mult.tmp %>% 
  dplyr::group_by(osler_id) %>% 
  dplyr::summarise(endpoint.u = min(endpoint), status.u = min(status))%>% 
  dplyr::select(osler_id, endpoint.u, status.u) %>% 
  unique() %>% 
  filter(endpoint.u > i)

obs.freq.new <- pprob_d4_new <- list()
pat_prob_save_new <- list()
for(v in 1:5){
  test_ids <- c(test_ids1[[v]],test_ids2[[v]],test_ids3[[v]])
  test_ids_4d <- test_ids[test_ids %in% dt.endpoints$osler_id]
  exp.freq <-  exp.freq.cv[[v]]
  obs.freq <- obs.freq.cv[[v]]
  pprob.save <- pat_prob_save[[v]]
    sNA <- id_remove$position[id_remove$cv == v]
    pprob.save.new <- pat_prob_save[[v]]
    if(length(sNA) !=0){
      for(s in 1:length(sNA)){
        nid <- as.numeric(sNA[s])
        pprob.save.new[[nid]] <- 0
        
        ## remove observed values
        sid <- test_ids_4d[nid]
        sid.end <- dt.endpoints %>% filter(osler_id == sid)
        endpoint.u <- sid.end$endpoint.u
        status.u <- sid.end$status.u
        obs.freq[rownames(obs.freq) == endpoint.u, colnames(obs.freq) == status.u] <- 
          obs.freq[rownames(obs.freq) == endpoint.u, colnames(obs.freq) == status.u] - 1
        obs.freq.new[[v]] <- obs.freq
      }
    }else{
      obs.freq.new[[v]] <- obs.freq.cv[[v]]
    }
    pat_prob_save_new[[v]] <- pprob.save.new
    pprob_d4_new[[v]] <- Reduce("+", pprob.save.new)
}

load("paperPrep/finalizedCode/figs/calibration/bkwr_di=2_M.RData")
load("paperPrep/finalizedCode/figs/calibration/bkwr_di=4_M.RData")
load("paperPrep/finalizedCode/figs/calibration/bkwr_di=8_M.RData")
sum(bkwr.save2[[3]][[1]])+sum(bkwr.save2[[3]][[2]])+
  sum(bkwr.save2[[3]][[3]])+sum(bkwr.save2[[3]][[4]]) + 
  sum(bkwr.save2[[3]][[5]])
sum(bkwr.save4[[3]][[1]])+sum(bkwr.save4[[3]][[2]])+
  sum(bkwr.save4[[3]][[3]])+sum(bkwr.save4[[3]][[4]]) + 
  sum(bkwr.save4[[3]][[5]])
sum(bkwr.save8[[3]][[1]])+sum(bkwr.save8[[3]][[2]])+
  sum(bkwr.save8[[3]][[3]])+sum(bkwr.save8[[3]][[4]]) + 
  sum(bkwr.save8[[3]][[5]])
sum(obs.freq.new[[1]]) + sum(obs.freq.new[[2]]) + 
  sum(obs.freq.new[[3]]) + sum(obs.freq.new[[4]]) + 
  sum(obs.freq.new[[5]])

exp.freq.cv <- pprob_d4_new
obs.freq.cv <- obs.freq.new
pat_prob_save <- pat_prob_save_new

dt.exp <- NULL
dt.pr <- NULL
dt.obs <- NULL
for(k in 1:5){
  #day <- (i+1):20
  day <- 0:20
  cv = rep(k, length(day))
  tmp.exp <- exp.freq.cv[[k]]; tmp.exp <- cbind(tmp.exp,day); tmp.exp <- cbind(tmp.exp,cv);
  tmp.obs <- obs.freq.cv[[k]]; tmp.obs <- cbind(tmp.obs,day); tmp.obs <- cbind(tmp.obs,cv);
  #tmp.pr <- pr.cv[[k]]; tmp.pr <- cbind(tmp.pr, day); tmp.pr <- cbind(tmp.pr,cv);
  dt.exp <- rbind(dt.exp, tmp.exp)
  dt.obs <- rbind(dt.obs, tmp.obs)
  #dt.pr <- rbind(dt.pr, tmp.pr)
}
dt.exp <- as.data.frame(dt.exp) %>% rename(Discharge = prior1, Die = prior2, Vent = prior3)%>% gather(event, exp, 1:3)
#dt.pr <- as.data.frame(dt.pr) %>% rename(Discharge = `1`, Die = `2`, Vent = `3`) %>% gather(event, pr, 1:3)
dt.obs <- as.data.frame(dt.obs) %>% rename(Discharge = `1`, Die = `2`, Vent = `3`) %>% gather(event, obs, 1:3)

dtp <- dt.exp %>% full_join(dt.obs) 

#pdf(paste0("paperPrep/figs/modelchecking/posterior_model_check_d",i,"_withCV.pdf"))
p <- dtp %>% 
  dplyr::group_by(day, event) %>% 
  dplyr::summarize(exp.tot = sum(exp, na.rm = T),
                   obs.tot = sum(obs, na.rm = T)) %>% 
  filter(day > i) %>% 
  ggplot()+
  geom_line(aes(x = day, y = obs.tot, color = event, linetype = 'obs')) + 
  geom_line(aes(x = day, y = exp.tot, color = event, linetype = 'exp')) +
  #facet_grid(cv~.)+
  theme_classic()+
  scale_color_manual(breaks = c("Discharge", "Die", "Vent"),
                     values = c("#2FB3CA","#F1564F","#90C13E"),
                     name=" ")+
  scale_linetype_manual(name = "", 
                        breaks = c("obs", "exp"),
                        values = c(1,2) ,
                        label = c("Observed", "Expected"))+
  ylab("Number of events")+
  xlab("Day of hospitalization")
#dev.off()

exp.freq <- Reduce("+", exp.freq.cv)[(i+2):21,]
obs.freq <- Reduce("+", obs.freq.cv)[(i+2):21,]
tmp <- sum((obs.freq-exp.freq)^2/exp.freq)
calib.dt <- tibble(`Chi-square stat` = round(tmp,2))

po <- p+geom_table_npc(data = calib.dt,
                       label = list(calib.dt),
                       npcx = 20, npcy = 1,
                       #table.theme = ttheme_gtlight(core = list(bg_params = list(fill = c("white")))),
                       table.theme = ttheme_gtminimal(),
                       size =4)

forw.save2 <- list(p, po,exp.freq.cv, obs.freq.cv, pat_prob_save)
save(forw.save2, 
     file = paste0("paperPrep/finalizedCode/figs/calibration/forw_di=", i,"_M.RData"))



