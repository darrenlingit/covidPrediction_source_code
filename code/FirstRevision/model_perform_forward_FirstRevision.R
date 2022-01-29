## model perform1 (cv) for forward model 
## check prior and posterior with cross validation
rm(list=ls())
nbiom <- 3
df.re <- 2
projT <- 20
nsim <- 30
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

# for(k in 1:5){
#   message(k)
#   test_ids <- c(test_ids1[[k]],test_ids2[[k]],test_ids3[[k]])
#   test_ids_4d <- test_ids[test_ids %in% dt.endpoints$osler_id]
#   
#   dt.pre.train <- dt.pre %>% filter(!(osler_id %in% test_ids))
#   
#   set.seed(2021)
#   fit.pre = MCMCglmm(cbind(sao2_norm, temp_norm, pulse_norm) ~
#                        -1 + trait + trait:(ns(interval,Boundary.knots = c(0,20), knots = c(1,3,6)) +
#                                              resp_bl + temp_bl + pulse_bl + sao2fio2_bl +
#                                              crp_bl + alc_bl + dd_bl + gfr_bl +
#                                              charlson_cat + demo5 + sex + bmi_cat +
#                                              smoking),
#                      random = ~ us(trait + trait:interval):osler_id,  #b0 + b1day for each outcome/trait
#                      #random = ~ us(trait + trait:ns(rev_day,knots = c(-3), Boundary.knots = c(-20,0))):osler_id,  #b0 + b1day for each outcome/trait
#                      rcov = ~ us(trait):units, #unstructured, coursenote pg67
#                      burnin = nb, nitt = nt, pr = T,
#                      family = rep("gaussian", 3),
#                      data = dt.pre.train %>% filter(interval != 0))
#   save(fit.pre, file = paste0("paperPrep/finalizedCode/model/fit.pre.forw.cv",k,".RData"))
# }

## i = 0 -------------
maxT = 20


dt1c <- dt1c %>% mutate(sao2_fio2_ratio_o = sao2_fio2_ratio_m,
                        temp_c_o = temp_c_m,
                        pulse_o = pulse_m)
dt1c <- dt1c %>%
  filter(complete.cases(temp_c_m, sao2_fio2_ratio_m, pulse_m)) %>% 
  mutate(sao2_norm = orderNorm(sao2_fio2_ratio_o)$x.t,
         temp_norm = orderNorm(temp_c_o)$x.t,
         pulse_norm = orderNorm(pulse_o)$x.t)

pprior_d0 <- list()
pr.cv <- list()
obs.freq.cv <- exp.freq.cv <- list()
pat_prob_save <- list() 
for(k in 1:5){
  message(k)
  test_ids <- c(test_ids1[[k]],test_ids2[[k]],test_ids3[[k]])
  ## split data and fit model using trains
  dt.pre.train <- dt.pre %>% filter(!(osler_id %in% test_ids))
  
  dt.pre.mult.train <- dt1c %>% filter(!(osler_id %in% test_ids)) #length(unique(dt.pre.mult.train$osler_id))
  dt.pre.mult.train <- dt.pre.mult.train %>% filter(!is.na(endpoint)) %>% filter(endpoint != 0)
  dt.pre.mult.train <- dt.pre.mult.train %>% filter(interval != 0)
  
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
  # mult.pre <- multinom(factor(Wt) ~ resp_bl + temp_bl + pulse_bl + sao2fio2_bl +
  #                                   crp_bl + alc_bl + dd_bl + gfr_bl +
  #                                   charlson_cat + demo5 + sex + bmi_cat +
  #                                   smoking + ns(interval, knots = c(4, 6 ,8,12,16)) +
  #                                   prev1_sao2 + prev1_temp + prev1_pulse +
  #                                   slope_sao2 + slope_temp + slope_pulse,
  #                                 na.action = na.omit,
  #                                 data = dt.pre.mult.train, maxit = 1000)
  # save(mult.pre, file = paste0("paperPrep/finalizedCode/model/mult.pre.forw.cv",k,".RData"))
  
  load(paste0("paperPrep/finalizedCode/model/fit.pre.forw.cv",k,".RData"))
  load(paste0("paperPrep/finalizedCode/model/mult.pre.forw.cv",k,".RData"))
  
  test_ids_0d <- test_ids[test_ids %in% dt.endpoints$osler_id]
  priors <- list()
  i=0
  for(s in 1:length(test_ids_0d)){
    if(s %% 10 == 0){
      print(s)
    }
    sid <- test_ids_0d[s]
    data <- dt.pre %>% filter(osler_id == sid) %>% filter(interval == 0)
    # filter(interval <= (i-1))
    ##  calculate joint distribution of [Y(t=0),...,Y(t=19)] ~ G(XB, ZDZ' + R)
    baseline <- dt.pre %>% 
      filter(osler_id == sid) %>% 
      dplyr::select(osler_id,resp_bl, temp_bl, pulse_bl, sao2fio2_bl,
                    crp_bl, alc_bl, dd_bl, gfr_bl, charlson_cat, demo5, sex, bmi_cat,smoking, status) %>% 
      unique()
    compdt <- do.call("rbind", replicate(projT+1, baseline, simplify = F))
    compdt$interval <- c(0:(projT))#newdt %>% mutate(interval = seq(0,j))
    #compdt <- compdt %>% mutate(status = factor(status , levels = c("1", "2", "3")))
    
    fm.Xi <- as.formula(paste0("~ns(interval,knots = c(1,3,6),
                                 Boundary.knots = c(0,20)) + resp_bl + temp_bl +
                                 pulse_bl + sao2fio2_bl + crp_bl + alc_bl + dd_bl + gfr_bl +
                                 charlson_cat + demo5 + sex + bmi_cat + smoking"))
    Xi1 <- model.matrix(fm.Xi, data = compdt)
    Xi2 <- model.matrix(fm.Xi, data = compdt)
    Xi3 <- model.matrix(fm.Xi, data = compdt)
    Xi <- as.matrix(bdiag(list(Xi1, Xi2, Xi3)))
    
    mod <- fit.pre
    Bhat <- fixef(mod, use = "mean")
    rn.Bhat <- rownames(Bhat)
    order.Bhat <- c(grep("sao2_norm", rn.Bhat),
                    grep("temp_norm", rn.Bhat),
                    grep("pulse_norm", rn.Bhat))
    Bhat <- matrix(Bhat[order.Bhat], ncol = 1)
    rownames(Bhat) <- rn.Bhat[order.Bhat]
    Yhati <- (Xi %*% Bhat)
    
    Zi1 =  model.matrix(~1+interval, compdt)
    Zi2 = model.matrix(~1+interval, compdt)
    Zi3 = model.matrix(~1+interval, compdt)
    Zi <- as.matrix(bdiag(list(Zi1, Zi2, Zi3)))
    ni <- length(Yhati)  
    smy <- summary(mod)
    ## process Ri 0120
    Rcov <- smy[["Rcovariances"]][,"post.mean"]
    rn.R <- names(Rcov)
    R = matrix(Rcov, nrow = nbiom)
    colnames(R) <- rownames(R) <-  rbind(lapply(str_split(rn.R[1:(nbiom)],":trait"), `[[`,1))
    vars <- c('sao2_fio2_ratio_m', 'temp_c_m', 'pulse_m')
    days <- compdt$interval
    nms.R <- do.call(paste0, expand.grid(vars, days))
    Ri <- as.matrix(bdiag(rep(list(R), length(days))))
    colnames(Ri) <- rownames(Ri) <- nms.R
    ord.Ri <- c(grep("sao2_fio2_ratio_m", colnames(Ri)),
                grep("temp_c_m", colnames(Ri)),
                grep("pulse_m", colnames(Ri)))
    Ri <- Ri[ord.Ri, ord.Ri]
    
    VCV <- mod$VCV # 600 by 90
    D <- apply(VCV, 2, mean)[1:(nbiom * df.re)^2] # (3*3)^2 # 3 biomarkers, random intercept+2df random slope = 3df for each biom
    rn.D <- names(D)
    D <- matrix(D, ncol = nbiom * df.re)
    colnames(D) <- rbind(lapply(str_split(rn.D[1:(nbiom * df.re)],":trait"), `[[`,1))
    # order.D <- c(grep("sao2_fio2_ratio_m", colnames(D)),
    #              grep("temp_c_m", colnames(D)),
    #              grep("pulse_m", colnames(D)))
    order.D <- c(grep("sao2_norm", colnames(D)),
                 grep("temp_norm", colnames(D)),
                 grep("pulse_norm", colnames(D)))
    D <- D[order.D,order.D]
    rownames(D) <- colnames(D)
    
    Vi = Zi %*% D %*% t(Zi) + Ri
    
    jm.tmp <- Yhati
    jv.tmp <- Vi
    
    nms <- c(paste0('sao2_fio2_ratio_m', '_day',c(0:(length(days)-1))),
             paste0('temp_c_m', '_day',c(0:(length(days)-1))),
             paste0('pulse_m', '_day',c(0:(length(days)-1))))
    rownames(jm.tmp) <- nms
    rownames(jv.tmp) <- colnames(jv.tmp) <-nms
    
    ## given first 8 days of data
    ## calculate [Yi9,...,Yi19|Yi0,...,Yi1]
    ## if only partially observed, conditioning only on observed Y's
    
    toMatch.Yit <- paste0('day', c(0:max((i-1),0)), "$")
    inx.Yit <- grep(paste(toMatch.Yit, collapse = "|"),nms)
    toMatch.Yit1ts <- paste0('day', c(max(i,1):(projT)),  "$")
    inx.Yit1ts <- grep(paste(toMatch.Yit1ts, collapse = "|"),nms)
    mu.Yit = matrix(jm.tmp[inx.Yit], ncol = 1)
    mu.Yit1ts = matrix(jm.tmp[inx.Yit1ts], ncol = 1)
    var.Yit = jv.tmp[inx.Yit,inx.Yit]
    var.Yit1ts = jv.tmp[inx.Yit1ts, inx.Yit1ts]
    cov.Yit1ts.Yit = jv.tmp[inx.Yit1ts,inx.Yit ]
    cov.Yit.Yit1ts =  t(cov.Yit1ts.Yit)
    
    if(nrow(data)==0){
      dtc <- baseline
      dtc$interval = 0
    }else{
      dtc <- data
    }
    
    last.intv <- max(i-1,0)
    if(nrow(dtc) != (last.intv + 1)){
      fk.length <- last.intv + 1
      fk.dtc <- data.frame(osler_id = rep(sid, fk.length),
                           interval = 0:last.intv)
      dtc <- dtc %>% full_join(fk.dtc, by = c('osler_id', 'interval')) %>% arrange(interval)
    }
    if(nrow(data)==0){
      Yit <- matrix(c(NA,NA,NA), ncol = 1)
      nms.Yit <- c(paste0('sao2_fio2_ratio_m', '_day',dtc$interval),
                   paste0('temp_c_m', '_day',dtc$interval),
                   paste0('pulse_m', '_day',dtc$interval))
      rownames(Yit) <- nms.Yit
    }else{
      Yit <- matrix(c(dtc$sao2_norm,
                      dtc$temp_norm,
                      dtc$pulse_norm), ncol = 1)
      nms.Yit <- c(paste0('sao2_fio2_ratio_m', '_day',dtc$interval),
                   paste0('temp_c_m', '_day',dtc$interval),
                   paste0('pulse_m', '_day',dtc$interval))
      rownames(Yit) <- nms.Yit
    }
    if(sum(is.na(Yit)) == 0){ ## if all measures complete
      cM <- mu.Yit1ts + cov.Yit1ts.Yit %*% solve(var.Yit) %*% (Yit - mu.Yit)
      cV <- var.Yit1ts - cov.Yit1ts.Yit %*% solve(var.Yit) %*% cov.Yit.Yit1ts
    }else if(sum(is.na(Yit)) > 0 & sum(!is.na(Yit)) != 0){  ## some days have some missings
      ind.c <- which(is.na(Yit))
      Yit.c <- matrix(Yit[-ind.c],ncol=1)
      mu.Yit.c <- matrix(mu.Yit[-ind.c], ncol = 1)
      var.Yit.c <- var.Yit[-ind.c, -ind.c]
      cov.Yit1ts.Yit.c <- cov.Yit1ts.Yit[, -ind.c]
      cov.Yit.Yit1ts.c =  t(cov.Yit1ts.Yit.c)
      
      mu.Yit1ts.c <- mu.Yit1ts
      var.Yit1ts.c <- var.Yit1ts
      
      cM <- mu.Yit1ts.c + cov.Yit1ts.Yit.c %*% solve(var.Yit.c) %*% (Yit.c - mu.Yit.c)
      cV <- var.Yit1ts.c - cov.Yit1ts.Yit.c %*% solve(var.Yit.c) %*% cov.Yit.Yit1ts.c
    }else if(sum(!is.na(Yit)) == 0){
      cM <- mu.Yit1ts
      cV <- var.Yit1ts
    }
    
    library(mvtnorm)
    set.seed(2021)
    Yfuture.sim <- t(rmvnorm(nsim, cM, cV)) ## interval 8 to interval 19
    
    priors.tmp <- list()
    for(v in 1:nsim){
      newdt <- compdt
      l_na <- length(c(0:max(0,(i-1))))
      Yfut <- Yfuture.sim[,v]
      newdt$sao2_sim <-  c(rep(NA, l_na), Yfut[1:(projT)])
      newdt$temp_sim <- c(rep(NA, l_na), Yfut[(projT+1):(2*(projT))])
      newdt$pulse_sim <- c(rep(NA, l_na),  Yfut[(2*(projT)+1):(3*(projT))])
      
      newdt <- newdt %>% 
        group_by(osler_id) %>% 
        arrange(interval) %>% 
        mutate(prev1_sao2 = lag(sao2_sim),
               prev1_temp = lag(temp_sim),
               prev1_pulse = lag(pulse_sim))%>% 
        mutate(prev2_sao2 = lag(sao2_sim, 2),
               prev2_temp = lag(temp_sim, 2),
               prev2_pulse = lag(pulse_sim, 2))
      newdt <- newdt %>% 
        group_by(osler_id) %>% 
        mutate(prev1_interval = lag(interval),
               prev2_interval = lag(interval, 2)) %>% 
        mutate(slope_sao2 = (prev1_sao2-prev2_sao2)/(prev1_interval-prev2_interval),
               slope_temp = (prev1_temp-prev2_temp)/(prev1_interval-prev2_interval),
               slope_pulse = (prev1_pulse-prev2_pulse)/(prev1_interval-prev2_interval))
      
      newdt <- newdt %>% left_join(data %>% dplyr::select(osler_id, interval, sao2_norm, temp_norm, pulse_norm), by = c("osler_id", "interval"))
      
      ## for day 0, use day1 and day2's simulated value to calculate any slopes
      setDT(newdt)
      newdt[, prev1_sao2 := na.locf(get("prev1_sao2"), na.rm=F, fromLast = T)]
      newdt[, prev1_temp := na.locf(get("prev1_temp"), na.rm=F, fromLast = T)]
      newdt[, prev1_pulse := na.locf(get("prev1_pulse"), na.rm=F, fromLast = T)]
      newdt[, slope_sao2 := na.locf(get("slope_sao2"), na.rm=F, fromLast = T)]
      newdt[, slope_temp := na.locf(get("slope_temp"), na.rm=F, fromLast = T)]
      newdt[, slope_pulse := na.locf(get("slope_pulse"), na.rm=F, fromLast = T)]
      newdt$prev1_sao2[newdt$interval == 0] <- NA
      
      preds <- as.data.frame(predict(mult.pre, type = 'prob', newdata = newdt))
      dtpreds <- cbind(newdt, preds)
      
      cond <- dtpreds %>% ungroup() %>% dplyr::select(`0`,`1`,`2`,`3`)
      tmp.prior <- dtpreds[-1,]
      newdt.prior <- tmp.prior %>% mutate(cumprod = cumprod(`0`), 
                                          cumsum1 = cumsum(`1`), 
                                          cumsum2 = cumsum(`2`), 
                                          cumsum3 = cumsum(`3`)) %>% 
        mutate(temp1 = `1`*`cumprod`/`0`,temp2 = `2`*`cumprod`/`0`, temp3 = `3`*`cumprod`/`0`,
               cir1 = cumsum(temp1), cir2 = cumsum(temp2), cir3 = cumsum(temp3))
      
      rem_p <- 1-sum(newdt.prior$temp1,newdt.prior$temp2,newdt.prior$temp3, na.rm=T)
      newdt.prior <- newdt.prior  %>% 
        dplyr::mutate(prior1 = temp1 + rem_p/(nbiom*(projT)), 
                      prior2 = temp2 + rem_p/(nbiom*(projT)),
                      prior3 = temp3 + rem_p/(nbiom*(projT)))
      newdt.prior <- rbind(dtpreds[which(rowSums(is.na(cond)) == ncol(cond)),], newdt.prior, fill = T)
      priors.tmp[[v]] <- newdt.prior %>% ungroup() %>% dplyr::select(prior1, prior2, prior3)
    }
    tmp <- Reduce(`+`, lapply(priors.tmp, function(x) replace(as.matrix(x), is.na(x),0)))
    tmp <- tmp/nsim
    priors[[s]] <- tmp
  }
  
  ## find observed freq - 1046 for i = 2
  dt.obs <- dt.endpoints %>%
    filter(osler_id %in% test_ids_0d)
  obs.freq <- table(dt.obs$endpoint.u, dt.obs$status.u)
  if(nrow(obs.freq)!= (maxT)+1){
    dayseq <- 0:(maxT)
    day.mis <- dayseq[which(!(dayseq %in% (rownames(obs.freq))))]
    row.mis <- matrix(rep(c(0,0,0), length(day.mis)), nrow = length(day.mis))
    rownames(row.mis) <- day.mis
    obs.freq <- rbind(obs.freq, row.mis)
    obs.freq <- obs.freq[order(as.numeric(rownames(obs.freq))),]
  }
  
  ## find expected freq using forward model
  library(purrr)
  priors.clean <- compact(priors)
  
  dt.forwpreds <- list()
  for(l in 1:length(priors.clean)){
    tmp <- as.data.frame(priors.clean[[l]])
    dt.forwpreds[[l]] <- tmp %>% ungroup() %>% dplyr::select(prior1, prior2, prior3)
  }
  #dt.exp <- Reduce("+", dt.forwpreds)
  
  exp.freq <- Reduce(`+`, lapply(dt.forwpreds, function(x) replace(x, is.na(x),0)))
  pat_prob_save[[k]] <- dt.forwpreds
  exp.freq.cv[[k]]<-exp.freq
  obs.freq.cv[[k]]<-obs.freq
  pr.cv[[k]] <- (obs.freq - exp.freq)/sqrt(exp.freq)
}

dt.exp <- NULL
dt.pr <- NULL
for(k in 1:5){
  day <- 0:20
  tmp.exp <- exp.freq.cv[[k]]; tmp.exp <- cbind(tmp.exp,day)
  tmp.pr <- pr.cv[[k]]; tmp.pr <- cbind(tmp.pr, day)
  dt.exp <- rbind(dt.exp, tmp.exp)
  dt.pr <- rbind(dt.pr, tmp.pr)
}
dt.exp <- as.data.frame(dt.exp) %>% rename(Discharge = prior1, Die = prior2, Vent = prior3)%>% gather(event, exp, 1:3)
dt.pr <- as.data.frame(dt.pr) %>% rename(Discharge = prior1, Die = prior2, Vent = prior3) %>% gather(event, pr, 1:3)
dtp <- full_join(dt.exp, dt.pr)

dt.exp <- NULL
dt.obs <- NULL
for(k in 1:5){
  day <- 0:20
  cv = rep(k, 21)
  tmp.exp <- exp.freq.cv[[k]]; tmp.exp <- cbind(tmp.exp,day); tmp.exp <- cbind(tmp.exp,cv);
  tmp.obs <- obs.freq.cv[[k]]; tmp.obs <- cbind(tmp.obs, day); tmp.obs <- cbind(tmp.obs,cv);
  dt.exp <- rbind(dt.exp, tmp.exp)
  dt.obs <- rbind(dt.obs, tmp.obs)
}
dt.exp <- as.data.frame(dt.exp) %>% rename(Discharge = prior1, Die = prior2, Vent = prior3)%>% gather(event, exp, 1:3)
dt.obs <- as.data.frame(dt.obs) %>% rename(Discharge = `1`, Die = `2`, Vent = `3`) %>% gather(event, obs, 1:3)

dtp <- full_join(dt.exp, dt.obs)
facetlabeller <- c("1" = "Test data 1","2" = "Test data 2","3" = "Test data 3",
                   "4" = "Test data 4","5" = "Test data 5")

p <- dtp %>% 
  dplyr::group_by(day, event) %>% 
  dplyr::summarize(exp.tot = sum(exp, na.rm = T),
                   obs.tot = sum(obs, na.rm = T)) %>% 
  filter(day !=0) %>% 
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

exp.freq <- Reduce("+", exp.freq.cv)[-1,]
obs.freq <- Reduce("+", obs.freq.cv)[-1,]
tmp <- sum((obs.freq-exp.freq)^2/exp.freq)
calib.dt <- tibble(`Chi-square stat` = round(tmp,2))

po <- p+geom_table_npc(data = calib.dt,
                       label = list(calib.dt),
                       npcx = 20, npcy = 1,
                       #table.theme = ttheme_gtlight(core = list(bg_params = list(fill = c("white")))),
                       table.theme = ttheme_gtminimal(),
                       size =4)

forw.save0.S30 <- list(p, po,exp.freq.cv, obs.freq.cv, pat_prob_save)
save(forw.save0.S30, 
     file = paste0("paperPrep/finalizedCode/figs/calibration/forw_di=", '0',"_S30.RData"))

## i = 2,4,8,.... ---------
i <- 8
## filter out patients have an event prior to day 4, keep only event occurred after day 4.
dt.endpoints <- dt.pre.mult.tmp %>% 
  dplyr::group_by(osler_id) %>% 
  dplyr::summarise(endpoint.u = min(endpoint), status.u = min(status))%>% 
  dplyr::select(osler_id, endpoint.u, status.u) %>% 
  unique() %>% 
  filter(endpoint.u > i)

pprob_d4 <- list()
pr.cv <- list()
obs.freq.cv <- exp.freq.cv <- list()
maxT = 20
pat_prob_save <- list()
for(v in 1:5){
  message(v)
  test_ids <- c(test_ids1[[v]],test_ids2[[v]],test_ids3[[v]])
  test_ids_4d <- test_ids[test_ids %in% dt.endpoints$osler_id]
  ## filter out patients missing all biomarkers at every day 
  testids_nonmis <- dt.pre %>%  ## keep patients with at least one marker oberved at any day
    filter(osler_id %in% test_ids_4d) %>%
    filter(complete.cases(temp_c_m, sao2_fio2_ratio_m, pulse_m)) %>% dplyr::select(osler_id) %>% unique()
  test_ids_4d <- test_ids_4d[test_ids_4d %in% testids_nonmis$osler_id]
  
  load(paste0("paperPrep/finalizedCode/model/fit.pre.forw.cv",v,".RData"))
  load(paste0("paperPrep/finalizedCode/model/mult.pre.forw.cv",v,".RData"))
  
  mod <- fit.pre
  VCV <- mod$VCV # 600 by 90
  D <- apply(VCV, 2, mean)[1:(nbiom * df.re)^2] # (3*3)^2 # 3 biomarkers, random intercept+2df random slope = 3df for each biom
  rn.D <- names(D)
  D <- matrix(D, ncol = nbiom * df.re)
  colnames(D) <- rbind(lapply(str_split(rn.D[1:(nbiom * df.re)],":trait"), `[[`,1))
  # order.D <- c(grep("sao2_fio2_ratio_m", colnames(D)),
  #              grep("temp_c_m", colnames(D)),
  #              grep("pulse_m", colnames(D)))
  order.D <- c(grep("sao2_norm", colnames(D)),
               grep("temp_norm", colnames(D)),
               grep("pulse_norm", colnames(D)))
  D <- D[order.D,order.D]
  rownames(D) <- colnames(D)
  priors <- list()
  for(s in 1:length(test_ids_4d)){
    if(s %% 10 == 0){
      print(s)
    }
    sid <- test_ids_4d[s]
    data <- dt.pre %>% 
      filter(osler_id == sid) %>% 
      filter(interval <= (i))
    
    ##  calculate joint distribution of [Y(t=0),...,Y(t=19)] ~ G(XB, ZDZ' + R)
    baseline <- dt.pre %>% 
      filter(osler_id == sid) %>% 
      dplyr::select(osler_id,resp_bl, temp_bl, pulse_bl, sao2fio2_bl,
                    crp_bl, alc_bl, dd_bl, gfr_bl, charlson_cat, demo5, sex, bmi_cat,smoking, status) %>% 
      unique()
    compdt <- do.call("rbind", replicate(projT+1, baseline, simplify = F))
    compdt$interval <- c(0:(projT))#newdt %>% mutate(interval = seq(0,j))
    #compdt <- compdt %>% mutate(status = factor(status , levels = c("1", "2", "3")))
    fm.Xi <- as.formula(paste0("~ns(interval,knots = c(1,3,6),
                                 Boundary.knots = c(0,20)) + resp_bl + temp_bl +
                                 pulse_bl + sao2fio2_bl + crp_bl + alc_bl + dd_bl + gfr_bl +
                                 charlson_cat + demo5 + sex + bmi_cat + smoking"))
    Xi1 <- model.matrix(fm.Xi, data = compdt)
    Xi2 <- model.matrix(fm.Xi, data = compdt)
    Xi3 <- model.matrix(fm.Xi, data = compdt)
    Xi <- as.matrix(bdiag(list(Xi1, Xi2, Xi3)))
    Bhat <- fixef(mod, use = "mean")
    rn.Bhat <- rownames(Bhat)
    # order.Bhat <- c(grep("sao2_fio2_ratio_m", rn.Bhat),
    #                 grep("temp_c_m", rn.Bhat),
    #                 grep("pulse_m", rn.Bhat))
    order.Bhat <- c(grep("sao2_norm", rn.Bhat),
                    grep("temp_norm", rn.Bhat),
                    grep("pulse_norm", rn.Bhat))
    Bhat <- matrix(Bhat[order.Bhat], ncol = 1)
    rownames(Bhat) <- rn.Bhat[order.Bhat]
    Yhati <- (Xi %*% Bhat)
    
    Zi1 =  model.matrix(~1+interval, compdt)
    Zi2 = model.matrix(~1+interval, compdt)
    Zi3 = model.matrix(~1+interval, compdt)
    Zi <- as.matrix(bdiag(list(Zi1, Zi2, Zi3)))
    ni <- length(Yhati)  
    smy <- summary(mod)
    ## process Ri 0120
    Rcov <- smy[["Rcovariances"]][,"post.mean"]
    rn.R <- names(Rcov)
    R = matrix(Rcov, nrow = nbiom)
    colnames(R) <- rownames(R) <-  rbind(lapply(str_split(rn.R[1:(nbiom)],":trait"), `[[`,1))
    vars <- c('sao2_fio2_ratio_m', 'temp_c_m', 'pulse_m')
    days <- compdt$interval
    nms.R <- do.call(paste0, expand.grid(vars, days))
    Ri <- as.matrix(bdiag(rep(list(R), length(days))))
    colnames(Ri) <- rownames(Ri) <- nms.R
    ord.Ri <- c(grep("sao2_fio2_ratio_m", colnames(Ri)),
                grep("temp_c_m", colnames(Ri)),
                grep("pulse_m", colnames(Ri)))
    Ri <- Ri[ord.Ri, ord.Ri]
    
    Vi = Zi %*% D %*% t(Zi) + Ri
    
    jm.tmp <- Yhati
    jv.tmp <- Vi
    
    nms <- c(paste0('sao2_fio2_ratio_m', '_day',c(0:(length(days)-1))),
             paste0('temp_c_m', '_day',c(0:(length(days)-1))),
             paste0('pulse_m', '_day',c(0:(length(days)-1))))
    rownames(jm.tmp) <- nms
    rownames(jv.tmp) <- colnames(jv.tmp) <-nms
    
    ## given first 8 days of data
    ## calculate [Yi9,...,Yi19|Yi0,...,Yi1]
    ## if only partially observed, conditioning only on observed Y's
    
    toMatch.Yit <- paste0('day', c(0:(i)), "$")
    inx.Yit <- grep(paste(toMatch.Yit, collapse = "|"),nms)
    toMatch.Yit1ts <- paste0('day', c((i+1):(projT)),  "$")
    inx.Yit1ts <- grep(paste(toMatch.Yit1ts, collapse = "|"),nms)
    mu.Yit = matrix(jm.tmp[inx.Yit], ncol = 1)
    mu.Yit1ts = matrix(jm.tmp[inx.Yit1ts], ncol = 1)
    var.Yit = jv.tmp[inx.Yit,inx.Yit]
    var.Yit1ts = jv.tmp[inx.Yit1ts, inx.Yit1ts]
    cov.Yit1ts.Yit = jv.tmp[inx.Yit1ts,inx.Yit ]
    cov.Yit.Yit1ts =  t(cov.Yit1ts.Yit)
    
    dtc <- data
    last.intv <- i
    if(nrow(dtc) != (last.intv + 1)){
      fk.length <- last.intv + 1
      fk.dtc <- data.frame(osler_id = rep(sid, fk.length),
                           interval = 0:last.intv)
      dtc <- dtc %>% full_join(fk.dtc, by = c('osler_id', 'interval')) %>% arrange(interval)
    }
    Yit <- matrix(c(dtc$sao2_norm, 
                    dtc$temp_norm, 
                    dtc$pulse_norm), ncol = 1)
    nms.Yit <- c(paste0('sao2_fio2_ratio_m', '_day',dtc$interval),
                 paste0('temp_c_m', '_day',dtc$interval),
                 paste0('pulse_m', '_day',dtc$interval))
    rownames(Yit) <- nms.Yit
    
    
    if(sum(is.na(Yit)) == 0){ ## if all measures complete
      cM <- mu.Yit1ts + cov.Yit1ts.Yit %*% solve(var.Yit) %*% (Yit - mu.Yit)
      cV <- var.Yit1ts - cov.Yit1ts.Yit %*% solve(var.Yit) %*% cov.Yit.Yit1ts
    }else if(sum(is.na(Yit)) > 0 & sum(!is.na(Yit)) != 0){  ## some days have some missings
      ind.c <- which(is.na(Yit))
      Yit.c <- matrix(Yit[-ind.c],ncol=1)
      mu.Yit.c <- matrix(mu.Yit[-ind.c], ncol = 1)
      var.Yit.c <- var.Yit[-ind.c, -ind.c]
      cov.Yit1ts.Yit.c <- cov.Yit1ts.Yit[, -ind.c]
      cov.Yit.Yit1ts.c =  t(cov.Yit1ts.Yit.c)
      
      mu.Yit1ts.c <- mu.Yit1ts
      var.Yit1ts.c <- var.Yit1ts
      
      cM <- mu.Yit1ts.c + cov.Yit1ts.Yit.c %*% solve(var.Yit.c) %*% (Yit.c - mu.Yit.c)
      cV <- var.Yit1ts.c - cov.Yit1ts.Yit.c %*% solve(var.Yit.c) %*% cov.Yit.Yit1ts.c
    }else if(sum(!is.na(Yit)) == 0){
      cM <- mu.Yit1ts
      cV <- var.Yit1ts
    }
    
    library(mvtnorm)
    set.seed(2021)
    Yfuture.sim <- t(rmvnorm(nsim, cM, cV)) ## interval 8 to interval 19
    
    priors.tmp <- list()
    for(k in 1:nsim){
      newdt <- compdt
      l_na <- length(c(0:(i)))
      Yfut <- Yfuture.sim[,k]
      newdt$sao2_sim <-  c(rep(NA, l_na), Yfut[1:(projT-i)])
      newdt$temp_sim <- c(rep(NA, l_na), Yfut[(projT-i+1):(2*(projT-i))])
      newdt$pulse_sim <- c(rep(NA, l_na),  Yfut[(2*(projT-i)+1):(3*(projT-i))])
      
      newdt <- newdt %>% 
        group_by(osler_id) %>% 
        arrange(interval) %>% 
        mutate(prev1_sao2 = lag(sao2_sim),
               prev1_temp = lag(temp_sim),
               prev1_pulse = lag(pulse_sim))%>% 
        mutate(prev2_sao2 = lag(sao2_sim, 2),
               prev2_temp = lag(temp_sim, 2),
               prev2_pulse = lag(pulse_sim, 2))
      newdt <- newdt %>% 
        group_by(osler_id) %>% 
        mutate(prev1_interval = lag(interval),
               prev2_interval = lag(interval, 2)) %>% 
        mutate(slope_sao2 = (prev1_sao2-prev2_sao2)/(prev1_interval-prev2_interval),
               slope_temp = (prev1_temp-prev2_temp)/(prev1_interval-prev2_interval),
               slope_pulse = (prev1_pulse-prev2_pulse)/(prev1_interval-prev2_interval))
      
      ## day 8 day 9's previous biomarker value and trend calculated from the observed data
      ## if observed missing, use simulated
      newdt$j.sao2.sim <- Yhati[1:(projT+1)]
      newdt$j.temp.sim <- Yhati[(projT+1+1):(2*(projT+1))]
      newdt$j.pulse.sim <- Yhati[(2*(projT+1)+1):(3*(projT+1))]
      newdt <- newdt %>% left_join(data %>% dplyr::select(osler_id, interval, sao2_norm, temp_norm, pulse_norm), by = c("osler_id", "interval"))
      
      j <- i+1
      crit_sao2 <- is.na(newdt$sao2_norm[newdt$interval %in% c((j-2), (j-1))])
      crit_temp <- is.na(newdt$temp_norm[newdt$interval %in% c((j-2), (j-1))])
      crit_pulse <- is.na(newdt$pulse_norm[newdt$interval %in% c((j-2), (j-1))])
      
      if(sum(crit_sao2) == 0){ ## none is missing
        newdt$prev1_sao2[newdt$interval == j] <- newdt$sao2_norm[newdt$interval == (j-1)]
        newdt$prev2_sao2[newdt$interval == j] <- newdt$sao2_norm[newdt$interval == (j-2)]
        newdt$prev2_sao2[newdt$interval == (j+1)] <- newdt$prev1_sao2[newdt$interval == j]
      }else if(sum(crit_sao2) == 1){
        if(which(crit_sao2 == TRUE) == 1){ ## missing i-2 observed value (i=4, missing day 2)
          newdt$prev1_sao2[newdt$interval == j] <- newdt$sao2_norm[newdt$interval == (j-1)]
          newdt$prev2_sao2[newdt$interval == j] <- newdt$j.sao2.sim[newdt$interval == (j-2)]
          newdt$prev2_sao2[newdt$interval == (j+1)] <- newdt$prev1_sao2[newdt$interval == j]
        }else if(which(crit_sao2 == TRUE) == 2){ ## missing i-1 observed value (i = 4, missing day 3)
          newdt$prev1_sao2[newdt$interval == j] <- newdt$j.sao2.sim[newdt$interval == (j-1)]
          newdt$prev2_sao2[newdt$interval == (j+1)] <- newdt$prev1_sao2[newdt$interval == j]
          newdt$prev2_sao2[newdt$interval == j] <- newdt$sao2_norm[newdt$interval == (j-2)] #because i-2 is not missing, supply prev2_sao2 with observed i-2 value
        }
      }else if(sum(crit_sao2) == 2){ ## both i-1 and i-2 missing, supply both with simulated value
        newdt$prev1_sao2[newdt$interval == j] <- newdt$j.sao2.sim[newdt$interval == (j-1)]
        newdt$prev2_sao2[newdt$interval == j] <- newdt$j.sao2.sim[newdt$interval == (j-2)]
        newdt$prev2_sao2[newdt$interval == (j+1)] <- newdt$prev1_sao2[newdt$interval == j]
        newdt$slope_sao2[newdt$interval == j] <- (newdt$prev1_sao2[newdt$interval == j]-newdt$prev2_sao2[newdt$interval == j])/(newdt$prev1_interval[newdt$interval == j]-newdt$prev2_interval[newdt$interval == j])
        newdt$slope_sao2[newdt$interval == (j+1)] <- (newdt$prev1_sao2[newdt$interval == (j+1)]-newdt$prev2_sao2[newdt$interval == (j+1)])/(newdt$prev1_interval[newdt$interval == (j+1)]-newdt$prev2_interval[newdt$interval == (j+1)])
      }
      
      
      if(sum(crit_temp) == 0){ ## none is missing
        newdt$prev1_temp[newdt$interval == j] <- newdt$temp_norm[newdt$interval == (j-1)]
        newdt$prev2_temp[newdt$interval == j] <- newdt$temp_norm[newdt$interval == (j-2)]
        newdt$prev2_temp[newdt$interval == (j+1)] <- newdt$prev1_temp[newdt$interval == j]
        
      }else if(sum(crit_temp) == 1){
        if(which(crit_temp == TRUE) == 1){ ## missing i-2 observed value (i=4, missing day 2)
          newdt$prev1_temp[newdt$interval == j] <- newdt$temp_norm[newdt$interval == (j-1)]
          newdt$prev2_temp[newdt$interval == j] <- newdt$j.temp.sim[newdt$interval == (j-2)]
          newdt$prev2_temp[newdt$interval == (j+1)] <- newdt$prev1_temp[newdt$interval == j]
        }else if(which(crit_temp == TRUE) == 2){ ## missing i-1 observed value (i = 4, missing day 3)
          newdt$prev1_temp[newdt$interval == j] <- newdt$j.temp.sim[newdt$interval == (j-1)]
          newdt$prev2_temp[newdt$interval == (j+1)] <- newdt$prev1_temp[newdt$interval == j]
          newdt$prev2_temp[newdt$interval == j] <- newdt$temp_norm[newdt$interval == (j-2)] #because i-2 is not missing, supply prev2_temp with observed i-2 value
        }
      }else if(sum(crit_temp) == 2){ ## both i-1 and i-2 missing, supply both with simulated value
        newdt$prev1_temp[newdt$interval == j] <- newdt$j.temp.sim[newdt$interval == (j-1)]
        newdt$prev2_temp[newdt$interval == j] <- newdt$j.temp.sim[newdt$interval == (j-2)]
        newdt$prev2_temp[newdt$interval == (j+1)] <- newdt$prev1_temp[newdt$interval == j]
        newdt$slope_temp[newdt$interval == j] <- (newdt$prev1_temp[newdt$interval == j]-newdt$prev2_temp[newdt$interval == j])/(newdt$prev1_interval[newdt$interval == j]-newdt$prev2_interval[newdt$interval == j])
        newdt$slope_temp[newdt$interval == (j+1)] <- (newdt$prev1_temp[newdt$interval == (j+1)]-newdt$prev2_temp[newdt$interval == (j+1)])/(newdt$prev1_interval[newdt$interval == (j+1)]-newdt$prev2_interval[newdt$interval == (j+1)])
      }  
      
      if(sum(crit_pulse) == 0){ ## none is missing
        newdt$prev1_pulse[newdt$interval == j] <- newdt$pulse_norm[newdt$interval == (j-1)]
        newdt$prev2_pulse[newdt$interval == j] <- newdt$pulse_norm[newdt$interval == (j-2)]
        newdt$prev2_pulse[newdt$interval == (j+1)] <- newdt$prev1_pulse[newdt$interval == j]
      }else if(sum(crit_pulse) == 1){
        if(which(crit_pulse == TRUE) == 1){ ## missing i-2 observed value (i=4, missing day 2)
          newdt$prev1_pulse[newdt$interval == j] <- newdt$pulse_norm[newdt$interval == (j-1)]
          newdt$prev2_pulse[newdt$interval == j] <- newdt$j.pulse.sim[newdt$interval == (j-2)]
          newdt$prev2_pulse[newdt$interval == (j+1)] <- newdt$prev1_pulse[newdt$interval == j]
        }else if(which(crit_pulse == TRUE) == 2){ ## missing i-1 observed value (i = 4, missing day 3)
          newdt$prev1_pulse[newdt$interval == j] <- newdt$j.pulse.sim[newdt$interval == (j-1)]
          newdt$prev2_pulse[newdt$interval == (j+1)] <- newdt$prev1_pulse[newdt$interval == j]
          newdt$prev2_pulse[newdt$interval == j] <- newdt$pulse_norm[newdt$interval == (j-2)] #because i-2 is not missing, supply prev2_pulse with observed i-2 value
        }
      }else if(sum(crit_pulse) == 2){ ## both i-1 and i-2 missing, supply both with simulated value
        newdt$prev1_pulse[newdt$interval == j] <- newdt$j.pulse.sim[newdt$interval == (j-1)]
        newdt$prev2_pulse[newdt$interval == j] <- newdt$j.pulse.sim[newdt$interval == (j-2)]
        newdt$prev2_pulse[newdt$interval == (j+1)] <- newdt$prev1_pulse[newdt$interval == j]
        newdt$slope_pulse[newdt$interval == j] <- (newdt$prev1_pulse[newdt$interval == j]-newdt$prev2_pulse[newdt$interval == j])/(newdt$prev1_interval[newdt$interval == j]-newdt$prev2_interval[newdt$interval == j])
        newdt$slope_pulse[newdt$interval == (j+1)] <- (newdt$prev1_pulse[newdt$interval == (j+1)]-newdt$prev2_pulse[newdt$interval == (j+1)])/(newdt$prev1_interval[newdt$interval == (j+1)]-newdt$prev2_interval[newdt$interval == (j+1)])
      }
      
      newdt$slope_sao2[newdt$interval == j] <- (newdt$prev1_sao2[newdt$interval == j]-newdt$prev2_sao2[newdt$interval == j])/(newdt$prev1_interval[newdt$interval == j]-newdt$prev2_interval[newdt$interval == j])
      newdt$slope_sao2[newdt$interval == (j+1)] <- (newdt$prev1_sao2[newdt$interval == (j+1)]-newdt$prev2_sao2[newdt$interval == (j+1)])/(newdt$prev1_interval[newdt$interval == (j+1)]-newdt$prev2_interval[newdt$interval == (j+1)])
      newdt$slope_temp[newdt$interval == j] <- (newdt$prev1_temp[newdt$interval == j]-newdt$prev2_temp[newdt$interval == j])/(newdt$prev1_interval[newdt$interval == j]-newdt$prev2_interval[newdt$interval == j])
      newdt$slope_temp[newdt$interval == (j+1)] <- (newdt$prev1_temp[newdt$interval == (j+1)]-newdt$prev2_temp[newdt$interval == (j+1)])/(newdt$prev1_interval[newdt$interval == (j+1)]-newdt$prev2_interval[newdt$interval == (j+1)])
      newdt$slope_pulse[newdt$interval == j] <- (newdt$prev1_pulse[newdt$interval == j]-newdt$prev2_pulse[newdt$interval == j])/(newdt$prev1_interval[newdt$interval == j]-newdt$prev2_interval[newdt$interval == j])
      newdt$slope_pulse[newdt$interval == (j+1)] <- (newdt$prev1_pulse[newdt$interval == (j+1)]-newdt$prev2_pulse[newdt$interval == (j+1)])/(newdt$prev1_interval[newdt$interval == (j+1)]-newdt$prev2_interval[newdt$interval == (j+1)])
      
      
      preds <- as.data.frame(predict(mult.pre, type = 'prob', newdata = newdt))
      dtpreds <- cbind(newdt, preds)
      
      cond <- dtpreds %>% ungroup() %>% dplyr::select(`0`,`1`,`2`,`3`)
      tmp.prior <- dtpreds[-(1:(i+1)),]
      newdt.prior <- tmp.prior %>% mutate(cumprod = cumprod(`0`), 
                                          cumsum1 = cumsum(`1`), 
                                          cumsum2 = cumsum(`2`), 
                                          cumsum3 = cumsum(`3`)) %>% 
        mutate(temp1 = `1`*`cumprod`/`0`,temp2 = `2`*`cumprod`/`0`, temp3 = `3`*`cumprod`/`0`,
               cir1 = cumsum(temp1), cir2 = cumsum(temp2), cir3 = cumsum(temp3))
      
      rem_p <- 1-sum(newdt.prior$temp1,newdt.prior$temp2,newdt.prior$temp3, na.rm=T)
      newdt.prior <- newdt.prior  %>% 
        dplyr::mutate(prior1 = temp1 + rem_p/(nbiom*(projT-i)), 
                      prior2 = temp2 + rem_p/(nbiom*(projT-i)),
                      prior3 = temp3 + rem_p/(nbiom*(projT-i)))
      newdt.prior <- rbind(dtpreds[which(rowSums(is.na(cond)) == ncol(cond)),], newdt.prior)
      priors.tmp[[k]] <- newdt.prior %>% ungroup() %>% dplyr::select(prior1, prior2, prior3)
    }
    tmp <- Reduce("+", priors.tmp)/nsim
    priors[[s]] <- tmp
  }
  
  ## find observed freq - 1046 for i = 2
  dt.obs <- dt.endpoints %>%
    filter(osler_id %in% test_ids_4d)
  obs.freq <- table(dt.obs$endpoint.u, dt.obs$status.u)
  if(ncol(obs.freq)!=3){
    eventseq <- 1:3
    event.mis <- eventseq[which(!(eventseq %in% (colnames(obs.freq))))]
    col.mis <- matrix(rep(0, nrow(obs.freq)), ncol = 1)
    colnames(col.mis) <- event.mis
    obs.freq <- cbind(obs.freq, col.mis)
    obs.freq <- obs.freq[,order(as.numeric(colnames(obs.freq)))]
  }
  if(nrow(obs.freq)!= (maxT+1)){
    dayseq <- 0:(maxT)
    day.mis <- dayseq[which(!(dayseq %in% (rownames(obs.freq))))]
    row.mis <- matrix(rep(c(0,0,0), length(day.mis)), nrow = length(day.mis))
    rownames(row.mis) <- day.mis
    obs.freq <- rbind(obs.freq, row.mis)
    obs.freq <- obs.freq[order(as.numeric(rownames(obs.freq))),]
  }
  
  ## find expected freq using forward model
  library(purrr)
  priors.clean <- compact(priors)
  
  dt.forwpreds <- list()
  for(l in 1:length(priors.clean)){
    tmp <- priors.clean[[l]]
    dt.forwpreds[[l]] <- tmp %>% ungroup() %>% dplyr::select(prior1, prior2, prior3)
  }
  #dt.exp <- Reduce("+", dt.forwpreds)
  
  exp.freq <- Reduce(`+`, lapply(dt.forwpreds, function(x) replace(x, is.na(x),0)))
  
  pat_prob_save[[v]] <- priors.clean
  exp.freq.cv[[v]]<-exp.freq
  obs.freq.cv[[v]]<-obs.freq
  pr.cv[[v]] <- (obs.freq - exp.freq)/sqrt(exp.freq)
  
}

dt.exp <- NULL
dt.pr <- NULL
dt.obs <- NULL
for(k in 1:5){
  day <- 0:20
  cv = rep(k, length(day))
  tmp.exp <- exp.freq.cv[[k]]; tmp.exp <- cbind(tmp.exp,day); tmp.exp <- cbind(tmp.exp,cv);
  tmp.obs <- obs.freq.cv[[k]]; tmp.obs <- cbind(tmp.obs,day); tmp.obs <- cbind(tmp.obs,cv);
  tmp.pr <- pr.cv[[k]]; tmp.pr <- cbind(tmp.pr, day); tmp.pr <- cbind(tmp.pr,cv);
  dt.exp <- rbind(dt.exp, tmp.exp)
  dt.obs <- rbind(dt.obs, tmp.obs)
  dt.pr <- rbind(dt.pr, tmp.pr)
}
dt.exp <- as.data.frame(dt.exp) %>% rename(Discharge = prior1, Die = prior2, Vent = prior3)%>% gather(event, exp, 1:3)
dt.pr <- as.data.frame(dt.pr) %>% rename(Discharge = prior1, Die = prior2, Vent = prior3) %>% gather(event, pr, 1:3)
dt.obs <- as.data.frame(dt.obs) %>% rename(Discharge = `1`, Die = `2`, Vent = `3`) %>% gather(event, obs, 1:3)

dtp <- dt.exp %>% full_join(dt.pr) %>% full_join(dt.obs) 

p <- dtp %>% 
  dplyr::group_by(day, event) %>% 
  dplyr::summarize(exp.tot = sum(exp, na.rm = T),
                   obs.tot = sum(obs, na.rm = T)) %>% 
  filter(day > i)%>%
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

exp.freq <- Reduce("+", exp.freq.cv)[-(1:(i+1)),]
obs.freq <- Reduce("+", obs.freq.cv)[-(1:(i+1)),]
tmp <- sum((obs.freq-exp.freq)^2/exp.freq)
calib.dt <- tibble(`Chi-square stat` = round(tmp,2))

po <- p+geom_table_npc(data = calib.dt,
                       label = list(calib.dt),
                       npcx = 20, npcy = 1,
                       #table.theme = ttheme_gtlight(core = list(bg_params = list(fill = c("white")))),
                       table.theme = ttheme_gtminimal(),
                       size =4)

forw.save8.S30 <- list(p, po,exp.freq.cv, obs.freq.cv, pat_prob_save)

save(forw.save8.S30, 
     file = paste0("paperPrep/finalizedCode/figs/calibration/forw_di=", i,"_S30.RData"))
