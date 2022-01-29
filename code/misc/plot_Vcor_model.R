## model performance check script 2
## Checking likelihood -  mean and variance structure
rm(list=ls())
nbiom <- 3
df.re <- 3
nb = 2000; nt = 8000
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
                             smoking = factor(smoking, levels = c("Current Smoker", "Former Smoker", "Never Smoker")),
                             charlson_cat = factor(charlson_cat, levels = c("0", "1-2" , "3-4", ">5")),
                             demo5 = factor(demo5, levels =   c("Young_White", "Old_White" , "Young_Black", "Old_Black", "Latinx/Other")),
                             bmi_cat = factor(bmi_cat, levels = c("<30", ">=30")))
# attr(ns(dt.pre$rev_day,4), "knots"); attr(ns(dt.pre$rev_day,4), "Boundary.knots")
# attr(ns(dt.pre$rev_day,2), "knots"); attr(ns(dt.pre$rev_day,2), "Boundary.knots")
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


load("paperPrep/finalizedCode/figs/calibration/cv_ids.RData")
test_ids1 = cv_ids[[1]]
test_ids2 = cv_ids[[2]]
test_ids3 = cv_ids[[3]]

## check mean ---------
## for one pair of testing and training set
# i=2 ## with first 2 days of data observed
maxT = 20
test_ids <- c(test_ids1[[1]],test_ids2[[1]],test_ids3[[1]],
              test_ids1[[2]],test_ids2[[2]],test_ids3[[2]],
              test_ids1[[3]],test_ids2[[3]],test_ids3[[3]],
              test_ids1[[4]],test_ids2[[4]],test_ids3[[4]],
              test_ids1[[5]],test_ids2[[5]],test_ids3[[5]])
load(paste0("paperPrep/finalizedCode/model/alldata_jtgaussian_model.RData"))
load(paste0("paperPrep/finalizedCode/model/alldata_multinom_model.RData"))
dt.endpoints <- dt.pre.mult.tmp %>% 
  dplyr::group_by(osler_id) %>% 
  dplyr::summarise(endpoint.u = min(endpoint), status.u = min(status))%>% 
  dplyr::select(osler_id, endpoint.u, status.u) %>% 
  unique() 
i <- 2
# get observed mean for 15 days-----
obs.mat1 <- obs.mat2 <- obs.mat3<-NULL
exp.mat1 <- exp.mat2 <- exp.mat3<-NULL
V.mat <- list()
set.seed(2022)
for(s in 1:length(test_ids)){
  message(s)
  sid <- test_ids[s]
  data <- dt.pre %>% filter(osler_id == sid) 
  
  data <- data %>% left_join(dt.endpoints, by = "osler_id")
  total.length <- unique(data$endpoint.u)+1  # in status3, patients' endpoint = first-vented-day -> rev_day = 0 is endpoint - 1
  ## fill in missing days
  if(unique(data$status) %in% c(1,2) & total.length != nrow(data)){
    baseline <- dt.pre %>% 
      filter(osler_id == sid) %>% 
      dplyr::select(osler_id,resp_bl, temp_bl, pulse_bl, sao2fio2_bl,
                    crp_bl, alc_bl, dd_bl, gfr_bl, charlson_cat, demo5, sex, bmi_cat,smoking) %>% 
      unique()
    newdt <- do.call("rbind", replicate(total.length, baseline, simplify = F))
    newdt$interval <- c(0:(total.length-1))#newdt %>% mutate(interval = seq(0,j))
    newdt$rev_day <- c(-(total.length-1):0)
    newdt$status <- unique(data$status)
    tmp <- data %>% dplyr::select(osler_id, interval, rev_day, status, sao2_fio2_ratio_o, pulse_o, temp_c_o,
                                  sao2_norm, temp_norm, pulse_norm, endpoint.u)
    
    data <- newdt %>% left_join(tmp)
  }else if(unique(data$status) %in% c(3) & (total.length-1) != nrow(data)){
    baseline <- dt.pre %>% 
      filter(osler_id == sid) %>% 
      dplyr::select(osler_id,resp_bl, temp_bl, pulse_bl, sao2fio2_bl,
                    crp_bl, alc_bl, dd_bl, gfr_bl, charlson_cat, demo5, sex, bmi_cat,smoking) %>% 
      unique()
    newdt <- do.call("rbind", replicate((total.length-1), baseline, simplify = F))
    newdt$interval <- c(0:(total.length-2))#newdt %>% mutate(interval = seq(0,j))
    newdt$rev_day <- c(-(total.length-2):0)
    newdt$status <- unique(data$status)
    tmp <- data %>% dplyr::select(osler_id,interval,rev_day, status, sao2_fio2_ratio_o, pulse_o, temp_c_o,
                                  sao2_norm, temp_norm, pulse_norm, endpoint.u)
    
    data <- newdt %>% left_join(tmp)
  }
  ## fill up to day 20
  if(nrow(data) < (maxT+1)){
    baseline <- dt.pre %>%
      filter(osler_id == sid) %>%
      dplyr::select(osler_id,resp_bl, temp_bl, pulse_bl, sao2fio2_bl,
                    crp_bl, alc_bl, dd_bl, gfr_bl, charlson_cat, demo5, sex, bmi_cat,smoking) %>%
      unique()
    newdt <- do.call("rbind", replicate(maxT+1, baseline, simplify = F))
    newdt$rev_day <- c(-maxT:0)
    newdt$interval <- c(0:maxT)
    newdt$status <- unique(data$status)
    tmp <- data %>% dplyr::select(osler_id, rev_day, status, sao2_fio2_ratio_o, pulse_o, temp_c_o,
                                  sao2_norm, temp_norm, pulse_norm)

    data <- newdt %>% left_join(tmp)
  }else{data <- data %>% filter(rev_day >= -20)}
  obs.tmp1 <- matrix(data$sao2_norm, nrow = 1); obs.tmp1[is.nan(obs.tmp1)] <- NA
  obs.tmp2 <- matrix(data$temp_norm, nrow = 1); obs.tmp2[is.nan(obs.tmp2)] <- NA
  obs.tmp3 <- matrix(data$pulse_norm, nrow = 1); obs.tmp3[is.nan(obs.tmp3)] <- NA
  
  obs.tmp1 <- matrix(obs.tmp1[1:(maxT+1)], nrow = 1)
  obs.tmp2 <- matrix(obs.tmp2[1:(maxT+1)], nrow = 1)
  obs.tmp3 <- matrix(obs.tmp3[1:(maxT+1)], nrow = 1)
  
  rownames(obs.tmp1) <- rownames(obs.tmp2) <- rownames(obs.tmp3) <- sid
  colnames(obs.tmp1) <- colnames(obs.tmp2) <- colnames(obs.tmp3) <- -maxT:0
  obs.mat1 <- rbind(obs.mat1, obs.tmp1)
  obs.mat2 <- rbind(obs.mat2, obs.tmp2)
  obs.mat3 <- rbind(obs.mat3, obs.tmp3)
  
  # get estiedmated trajectory for 20 days (XB + Zb) ---
  mod <- fit.pre
  smy <- summary(mod)
  fm <-mod$Fixed$formula
  fe.spline <- ns(dt.pre$rev_day, 5) 
  re.spline <- ns(dt.pre$rev_day, 2) 
  Bhat <- fixef(mod, use = "mean")
  rn.Bhat <- rownames(Bhat)
  order.Bhat <- c(grep("sao2_norm", rn.Bhat),
                  grep("temp_norm", rn.Bhat),
                  grep("pulse_norm", rn.Bhat))
  Bhat <- matrix(Bhat[order.Bhat], ncol = 1)
  rownames(Bhat) <- rn.Bhat[order.Bhat]
  VCV <- mod$VCV # 600 by 90
  D <- apply(VCV, 2, mean)[1:(nbiom * df.re)^2] # (3*3)^2 # 3 biomarkers, random intercept+2df random slope = 3df for each biom
  rn.D <- names(D)
  D <- matrix(D, ncol = nbiom * df.re)
  colnames(D) <- rbind(lapply(str_split(rn.D[1:(nbiom * df.re)],":trait"), `[[`,1))
  order.D <- c(grep("sao2_norm", colnames(D)),
               grep("temp_norm", colnames(D)),
               grep("pulse_norm", colnames(D)))
  D <- D[order.D,order.D]
  rownames(D) <- colnames(D)
  
  fm.Zi <- as.formula(paste0("~1 + ns(rev_day,knots = c(-3), Boundary.knots = c(-20,0))"))
  Zi1 <- model.matrix(fm.Zi, data)
  Zi2 <- model.matrix(fm.Zi, data)
  Zi3 <- model.matrix(fm.Zi, data)
  Zi <- as.matrix(bdiag(list(Zi1, Zi2, Zi3)))

  
  # check variance --
  Rcov <- smy[["Rcovariances"]][,"post.mean"]
  rn.R <- names(Rcov)
  R = matrix(Rcov, nrow = nbiom)
  colnames(R) <- rownames(R) <-  rbind(lapply(str_split(rn.R[1:(nbiom)],":trait"), `[[`,1))
  vars <- c('sao2_fio2_ratio_m', 'temp_c_m', 'pulse_m')
  days <- c(0:maxT)
  nms.R <- do.call(paste0, expand.grid(vars, days))
  Ri <- as.matrix(bdiag(rep(list(R), maxT+1)))
  colnames(Ri) <- rownames(Ri) <- nms.R
  ord.Ri <- c(grep("sao2_fio2_ratio_m", colnames(Ri)),
              grep("temp_c_m", colnames(Ri)),
              grep("pulse_m", colnames(Ri)))
  Ri <- Ri[ord.Ri, ord.Ri]
  
  Vi = Zi %*% D %*% t(Zi) + Ri
  ## fill in missing patterns for Vi
  na_inx <- c(which(is.na(obs.tmp1)), (maxT+1)+(which(is.na(obs.tmp2))),(2*(maxT+1))+which(is.na(obs.tmp3)))
  Vi[na_inx,na_inx] <- NA
  
  V.mat[[s]] <- Vi
  names(V.mat)[s] <- sid
}

dt.endpoints.test <- dt.endpoints %>% filter(osler_id %in% test_ids)

V.ave <- list()
dt.ave <- NULL
for(e in 1:3){
  tmp.ave <- NULL
  id.select <- dt.endpoints.test$osler_id[dt.endpoints.test$status.u == e]
  tmp.exp1 <- exp.mat1[rownames(exp.mat1) %in% id.select,]
  tmp.exp2 <- exp.mat2[rownames(exp.mat2) %in% id.select,]
  tmp.exp3 <- exp.mat3[rownames(exp.mat3) %in% id.select,]
  
  tmp.obs1 <- obs.mat1[rownames(obs.mat1) %in% id.select,]
  tmp.obs2 <- obs.mat2[rownames(obs.mat2) %in% id.select,]
  tmp.obs3 <- obs.mat3[rownames(obs.mat3) %in% id.select,]
  
  exp.ave1 <- colMeans(tmp.exp1, na.rm = T)
  exp.ave2 <- colMeans(tmp.exp2, na.rm = T)
  exp.ave3 <- colMeans(tmp.exp3, na.rm = T)
  
  obs.ave1 <- colMeans(tmp.obs1, na.rm = T)
  obs.ave2 <- colMeans(tmp.obs2, na.rm = T)
  obs.ave3 <- colMeans(tmp.obs3, na.rm = T)
  
  event = rep(e, length(obs.ave1))
  tmp.ave  <- cbind(tmp.ave,exp.ave1,exp.ave2,exp.ave3,
                    obs.ave1, obs.ave2, obs.ave3,
                    event)
  dt.ave <- rbind(dt.ave,tmp.ave)
}

tmp.V <- V.mat
a <- do.call(cbind,tmp.V)
dim(a) <- c(3*(maxT+1),3*(maxT+1), length(tmp.V))
V.ave <- apply(a, c(1,2), mean, na.rm =T)
V.cor <- diag(1/sqrt(diag(V.ave))) %*% V.ave %*% diag(1/sqrt(diag(V.ave))) 
rownames(V.cor) <- colnames(V.cor) <- rep(c(-15:0), 3)
pdf("paperPrep/figs/modelchecking/model_Vcor.pdf")
corrplot(V.cor, tl.col = "black", tl.cex = 0.6, title = "Model based V-cor")
dev.off()

dt.ave <- as.data.frame(dt.ave) %>% mutate(diff.ave1 = obs.ave1 - exp.ave1,
                                           diff.ave2 = obs.ave2 - exp.ave2,
                                           diff.ave3 = obs.ave3 - exp.ave3,
                                           rev_day = rownames(dt.ave))
dtp <- dt.ave %>% 
  dplyr::select(event, rev_day, diff.ave1, diff.ave2, diff.ave3) %>% 
  dplyr::rename(biom1 = diff.ave1, biom2=diff.ave2, biom3=diff.ave3) %>% 
  gather(biomarker,diff,biom1:biom3)

dtp.obs <- dt.ave %>% 
  dplyr::select(event, rev_day, obs.ave1, obs.ave2, obs.ave3) %>% 
  dplyr::rename(biom1 = obs.ave1, biom2=obs.ave2, biom3=obs.ave3) %>% 
  gather(biomarker,obs,biom1:biom3)

dtp.exp <- dt.ave %>% 
  dplyr::select(event, rev_day, exp.ave1, exp.ave2, exp.ave3) %>% 
  dplyr::rename(biom1 = exp.ave1, biom2=exp.ave2, biom3=exp.ave3) %>% 
  gather(biomarker,exp,biom1:biom3)

pdf("paperPrep/figs/modelchecking/mean_checking.pdf",width = 12, height = 10)
dtp <-dtp %>% mutate(rev_day = as.numeric(rev_day))
dtp <- dtp %>%  mutate(diff = ifelse(biomarker == 'biom2', diff*10, diff))
ggplot()+
  geom_point(data = dtp, aes(x=rev_day, y = diff))+
  geom_line(data = dtp, aes(x=rev_day, y = diff, group = interaction(event, biomarker)))+
  facet_grid(event~biomarker, scales = "free")+
  theme_classic()+
  geom_hline(yintercept = 0)
dev.off()

pdf("paperPrep/figs/modelchecking/mean_checking2.pdf",width = 12, height = 10)
ggplot()+
  geom_point(data = dtp.obs, aes(x=rev_day, y = obs, 
                                 group = interaction(event, biomarker), color = "observed"),size = 1)+
  geom_line(data = dtp.obs, aes(x=rev_day, y = obs, 
                                group = interaction(event, biomarker), color = "observed"),size = 1)+
  geom_point(data = dtp.exp, aes(x=rev_day, y = exp, 
                                 group = interaction(event, biomarker), color = "expected"), size = 1)+
  geom_line(data = dtp.exp, aes(x=rev_day, y = exp, 
                                group = interaction(event, biomarker), color = "expected"), size = 1)+
  facet_grid(event~biomarker)+
  theme_classic()
dev.off()
