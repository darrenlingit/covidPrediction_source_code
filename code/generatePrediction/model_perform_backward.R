## model performance check script 1
## check prior and posterior probablity (calibration)
## with cross validation for backward model
rm(list=ls())
nbiom <- 3
df.re <- 3
nb = 2000; nt = 8000
source("00_loadlibs.R")
setwd(wd)

##load data
dt1c <- data

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
# cv_ids <- list(test_ids1=test_ids1,test_ids2=test_ids2, test_ids3=test_ids3)
# save(cv_ids, 
#      file = "paperPrep/finalizedCode/figs/calibration/cv_ids.RData")
## check prior model -------------
maxT = 20
pprior_d0 <- list()
pr.cv <- list()
obs.freq.cv <- exp.freq.cv <- list()
pat_prob_save <- list()
for(k in 1:5){
  message(k)
  test_ids <- c(test_ids1[[k]],test_ids2[[k]],test_ids3[[k]])
  ## split data and fit model using trains
  #dt.pre.train <- dt.pre %>% filter(!(osler_id %in% test_ids))
  
  dt.pre.mult.train <- dt1c %>% filter(!(osler_id %in% test_ids)) #length(unique(dt.pre.mult.train$osler_id))
  dt.pre.mult.train <- dt.pre.mult.train %>% filter(!is.na(endpoint)) %>% filter(endpoint != 0)
  dt.pre.mult.train <- dt.pre.mult.train %>% filter(interval != 0)
  mult.pre <- multinom(factor(Wt) ~ resp_bl + temp_bl + pulse_bl + sao2fio2_bl +
                         crp_bl + alc_bl + dd_bl + gfr_bl +
                         charlson_cat + demo5 + sex + bmi_cat +
                         smoking + ns(interval, knots = c(2,4, 6 ,12)),
                       data = dt.pre.mult.train, maxit = 1000)
  save(mult.pre, file = paste0("model/mult.pre.cv",k,".RData"))

  load(paste0("model/mult.pre.cv",k,".RData"))
  
  test_ids_0d <- test_ids[test_ids %in% dt.endpoints$osler_id]
  prior.save <- list()
  for(s in 1:length(test_ids_0d)){
    sid <- test_ids_0d[s]
    data <- dt.pre %>% filter(osler_id == sid)
    data <- data %>% left_join(dt.endpoints, by = "osler_id")
    data.prior <- data %>% dplyr::select(osler_id,resp_bl, temp_bl, pulse_bl, sao2fio2_bl,
                                         crp_bl, alc_bl, dd_bl, gfr_bl,
                                         charlson_cat, demo5, sex, bmi_cat,smoking) %>% unique()
    
    j=maxT
    newdt <- data.prior
    newdt <- do.call("rbind", replicate(j, newdt, simplify = F))
    newdt$interval <- c(1:(j))#newdt %>% mutate(interval = seq(0,j))
    newdt$rev_day <- c((-j):-1)
    dtl<-newdt
    
    dtl <- cbind(dtl, predict(mult.pre, type = 'prob', newdata = dtl))
    cir <- dtl %>% mutate(cumprod = cumprod(`0`), 
                          cumsum1 = cumsum(`1`), 
                          cumsum2 = cumsum(`2`), 
                          cumsum3 = cumsum(`3`)) %>% 
      mutate(temp1 = `1`*`cumprod`/`0`,temp2 = `2`*`cumprod`/`0`, temp3 = `3`*`cumprod`/`0`,
             cir1 = cumsum(temp1), cir2 = cumsum(temp2), cir3 = cumsum(temp3)) %>% 
      dplyr::rename(prior1 = temp1, 
                    prior2 = temp2,
                    prior3 = temp3)
    
    xi <- cir %>% dplyr::select(osler_id, interval, rev_day, prior1, prior2, prior3)
    if(sum(xi$prior1, xi$prior2,xi$prior3) != 1){
      prem <- (1- sum(xi$prior1, xi$prior2,xi$prior3))/(3*nrow(xi))
      xi <- xi %>% mutate(prior1 = prior1 + prem,
                          prior2 = prior2 + prem,
                          prior3 = prior3 + prem)
    }
    prior.save[[s]] <- xi 
  }
  pprior.save <- list()
  for(s in 1:length(prior.save)){
    tmp <- prior.save[[s]]
    tmp <- tmp %>% dplyr::select(prior1:prior3)
    pprior.save[[s]] <- as.matrix(tmp)
  }
  pprior_d0[[k]] <- Reduce("+", pprior.save)
  
  #observed table
  dt.endpoints.test <- dt.endpoints %>% filter(osler_id %in% test_ids_0d)
  obs.freq <- table(dt.endpoints.test$endpoint.u, dt.endpoints.test$status.u)
  if(nrow(obs.freq)!= maxT){
    day.mis <- which(!(1:maxT %in% (rownames(obs.freq))))
    row.mis <- matrix(rep(c(0,0,0), length(day.mis)), nrow = length(day.mis))
    rownames(row.mis) <- day.mis
    obs.freq <- rbind(obs.freq, row.mis)
    obs.freq <- obs.freq[order(as.numeric(rownames(obs.freq))),]
  }
  exp.freq <- pprior_d0[[k]]
  
  pat_prob_save[[k]] <- pprior.save
  exp.freq.cv[[k]]<-exp.freq
  obs.freq.cv[[k]]<-obs.freq
  pr.cv[[k]] <- (obs.freq - exp.freq)/sqrt(exp.freq)
}

dt.exp <- NULL
dt.pr <- NULL
for(k in 1:5){
  day <- 1:20
  tmp.exp <- exp.freq.cv[[k]]; tmp.exp <- cbind(tmp.exp,day)
  tmp.pr <- pr.cv[[k]]; tmp.pr <- cbind(tmp.pr, day)
  dt.exp <- rbind(dt.exp, tmp.exp)
  dt.pr <- rbind(dt.pr, tmp.pr)
}
dt.exp <- as.data.frame(dt.exp) %>% rename(Discharge = prior1, Die = prior2, Vent = prior3)%>% gather(event, exp, 1:3)
dt.pr <- as.data.frame(dt.pr) %>% rename(Discharge = `1`, Die = `2`, Vent = `3`) %>% gather(event, pr, 1:3)
dtp <- full_join(dt.exp, dt.pr)

dt.exp <- NULL
dt.obs <- NULL
for(k in 1:5){
  day <- 1:20
  cv = rep(k, 20)
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
  ggplot()+
  geom_line(aes(x = day, y = obs.tot, color = event, linetype = 'obs')) + 
  geom_line(aes(x = day, y = exp.tot, color = event, linetype = 'exp')) +
  # facet_grid(cv ~. ,labeller = as_labeller(facetlabeller))+
  theme_classic()+
  scale_color_manual(breaks = c("Discharge", "Die", "Vent"),
                     values = c("#2FB3CA","#F1564F","#90C13E"),
                     name=" ")+
  scale_linetype_manual(name = "", 
                        breaks = c("obs", "exp"),
                        values = c(1,2) ,
                        label = c("Observed", "Expected"))+
  xlab("Day of Hospitalization")+
  ylab("Number of events")

exp.freq <- Reduce("+", exp.freq.cv)
obs.freq <- Reduce("+", obs.freq.cv)
tmp <- sum((obs.freq-exp.freq)^2/exp.freq)
calib.dt <- tibble(`Chi-square stat` = round(tmp,2))

po <- p+geom_table_npc(data = calib.dt,
                       label = list(calib.dt),
                       npcx = 20, npcy = 1,
                       #table.theme = ttheme_gtlight(core = list(bg_params = list(fill = c("white")))),
                       table.theme = ttheme_gtminimal(),
                       size =4)

bkwr.save0 <- list(p, po,exp.freq.cv, obs.freq.cv, pat_prob_save)
save(bkwr.save0, 
     file = paste0("model/bkwr_di=", "0",".RData"))




## check posterior -----------------
attr(ns(dt.pre$rev_day,5), "Boundary.knots"); attr(ns(dt.pre$rev_day,5), "knots")
attr(ns(dt.pre$rev_day,2), "Boundary.knots");attr(ns(dt.pre$rev_day,2), "knots")
for(k in 1:5){
  message(k)
  test_ids <- c(test_ids1[[k]],test_ids2[[k]],test_ids3[[k]])
  test_ids_4d <- test_ids[test_ids %in% dt.endpoints$osler_id]
  
  dt.pre.train <- dt.pre %>% filter(!(osler_id %in% test_ids)) 
  
  set.seed(2021)
  fit.pre = MCMCglmm(cbind(sao2_norm, temp_norm, pulse_norm) ~
                       -1 + trait + trait:(status*ns(rev_day,knots = c(-2, -5, -10,-15), Boundary.knots = c(-20,0)) +
                                             resp_bl + temp_bl + pulse_bl + sao2fio2_bl +
                                             crp_bl + alc_bl + dd_bl + gfr_bl +
                                             charlson_cat + demo5 + sex + bmi_cat +
                                             smoking),
                     random = ~ us(trait + trait:ns(rev_day,knots = c(-3), Boundary.knots = c(-20,0))):osler_id,  #b0 + b1day for each outcome/trait
                     rcov = ~ us(trait):units, #unstructured, coursenote pg67
                     burnin = nb, nitt = nt, pr = T,
                     family = rep("gaussian", 3),
                     data = dt.pre.train%>% filter(interval != 0))
  save(fit.pre, file = paste0("model/fit.pre.cv",k,"_diffknots.RData"))
}

## day 4, i=2, 4, 8, ...)
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
v.df <- 6
for(v in 1:5){
  message(v)
  test_ids <- c(test_ids1[[v]],test_ids2[[v]],test_ids3[[v]])
  test_ids_4d <- test_ids[test_ids %in% dt.endpoints$osler_id]
  ## filter out patients missing all biomarkers at every day 
  testids_nonmis <- dt.pre %>%  ## keep patients with at least one marker oberved at any day
    filter(osler_id %in% test_ids_4d) %>%
    filter(complete.cases(temp_c_m, sao2_fio2_ratio_m, pulse_m)) %>% dplyr::select(osler_id) %>% unique()
  test_ids_4d <- test_ids_4d[test_ids_4d %in% testids_nonmis$osler_id]
  
  dt.pre.train <- dt.pre %>% filter(!(osler_id %in% test_ids_4d))
  
  load(paste0("model/fit.pre.cv",v,".RData"))
  load(paste0("model/mult.pre.cv",v,".RData"))
  mod <- fit.pre
  smy <- summary(mod)
  fm <-mod$Fixed$formula
  fe.spline <- ns(dt.pre$rev_day, 5) 
  re.spline <- ns(dt.pre$rev_day, 2) 
  Bhat <- fixef(mod, use = "mean")
  rn.Bhat <- rownames(Bhat)
  # order.Bhat <- c(grep("sao2_fio2_ratio_m", rn.Bhat),
  #                 grep("temp_c_m", rn.Bhat),
  #                 grep("pulse_m", rn.Bhat))
  order.Bhat <- c(grep("sao2_norm", rn.Bhat),
                  grep("temp_norm", rn.Bhat),
                  grep("pulse_norm",rn.Bhat))
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
  
  dtprob.save <- list()
  for(s in 1:length(test_ids_4d)){
    sid <- test_ids_4d[s]
    data <- dt.pre %>% filter(osler_id == sid) 
    
    data <- data %>% left_join(dt.endpoints, by = "osler_id")
    total.length <- unique(data$endpoint.u)+1  # in status3, patients' endpoint = first-vented-day -> rev_day = 0 is endpoint - 1
    if(unique(data$status) %in% c(1,2) & total.length != nrow(data)){
      baseline <- dt.pre %>% 
        filter(osler_id == sid) %>% 
        dplyr::select(osler_id,resp_bl, temp_bl, pulse_bl, sao2fio2_bl,
                      crp_bl, alc_bl, dd_bl, gfr_bl, charlson_cat, demo5, sex, bmi_cat,smoking) %>% 
        unique()
      newdt <- do.call("rbind", replicate(total.length, baseline, simplify = F))
      newdt$interval <- c(0:(total.length-1))#newdt %>% mutate(interval = seq(0,j))
      newdt$rev_day <- c(-(total.length-1):0)
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
      tmp <- data %>% dplyr::select(osler_id,interval,rev_day, status, sao2_fio2_ratio_o, pulse_o, temp_c_o,
                                    sao2_norm, temp_norm, pulse_norm, endpoint.u)
      
      data <- newdt %>% left_join(tmp)
    }
    
    ## likelihood -----------
    ## on day 4
    Lt <-  Lt.wt  <- NULL
    Lg <- Lg.wt <- NULL
    projT <- maxT-i
    for(j in 1:projT){
      lkh <- lkh.wt <-  NULL
      lkh.g <- lkh.g.wt <- NULL
      select =  data$interval <= i # i=1, j=0->if i'm on day0, I should be just using day 0 's data to project day 1,...,14
      dtl <- data[select,]
      dtl <- dtl %>% mutate(endpoint = i+j, rev_day = interval - endpoint)
      Zi1 =  model.matrix(~1+ns(rev_day, Boundary.knots = attr(re.spline,"Boundary.knots"),
                                knots = attr(re.spline,"knots")), dtl)
      Zi2 = model.matrix(~1+ns(rev_day, Boundary.knots = attr(re.spline,"Boundary.knots"),
                               knots = attr(re.spline,"knots")), dtl)
      Zi3 = model.matrix(~1+ns(rev_day, Boundary.knots = attr(re.spline,"Boundary.knots"),
                               knots = attr(re.spline,"knots")), dtl)
      Zi <- as.matrix(bdiag(list(Zi1, Zi2, Zi3)))
      #Zi <- as.matrix(bdiag(list(Zi1, Zi2)))
      ni <- (i+1)*nbiom  
      ## process Ri 0120
      Rcov <- smy[["Rcovariances"]][,"post.mean"]
      rn.R <- names(Rcov)
      R = matrix(Rcov, nrow = nbiom)
      colnames(R) <- rownames(R) <-  rbind(lapply(str_split(rn.R[1:(nbiom)],":trait"), `[[`,1))
      vars <- c('sao2_fio2_ratio_m', "temp_c_m", 'pulse_m')
      days <- c(0:(i))
      nms.R <- do.call(paste0, expand.grid(vars, days))
      Ri <- as.matrix(bdiag(rep(list(R), i+1)))
      colnames(Ri) <- rownames(Ri) <- nms.R
      ord.Ri <- c(grep("sao2_fio2_ratio_m", colnames(Ri)),
                  grep("temp_c_m", colnames(Ri)),
                  grep("pulse_m", colnames(Ri)))
      Ri <- Ri[ord.Ri, ord.Ri]
      
      Yi.o <- matrix(c(dtl$sao2_norm, 
                       dtl$temp_norm, 
                       dtl$pulse_norm), ncol = 1)
      Vi.o = Zi %*% D %*% t(Zi) + Ri
      Vi.invsqrtm.o = solve(sqrtm(Vi.o))

      row1 <- c(rep(0,(i+1)-1), 1, rep(0, (nbiom-1)*(i+1)))
      row2 <- c(rep(0,(nbiom-1)*(i+1)-1), 1, rep(0, (i+1)))
      row3 <- c(rep(0,nbiom*(i+1)-1),1)
      row4 <- c(rep(0, (i+1)-3), -1/2, 0, 1/2, rep(0,  (nbiom-1)*(i+1)))
      row5 <- c(rep(0,(i+1)), rep(0, (i+1)-3), -1/2, 0, 1/2, rep(0,(i+1)))
      row6 <- c(rep(0,(nbiom-1)*(i+1)), rep(0, (i+1)-3), -1/2, 0, 1/2)
      
      M <- rbind(row1, row2, row3, row4, row5, row6)
      
      for(k in 1:3){
        dtl2 <- dtl %>% mutate(status = factor(k , levels = c("1", "2", "3")))
        fm.Xi <- as.formula(paste0("~status * ns(rev_day,knots = c(-2, -5, -10,-15),
                                 Boundary.knots = c(-20,0)) + resp_bl + temp_bl +
                                 pulse_bl + sao2fio2_bl + crp_bl + alc_bl + dd_bl + gfr_bl +
                                 charlson_cat + demo5 + sex + bmi_cat + smoking"))

        Xi1 <- model.matrix(fm.Xi, data = dtl2) #;colnames(Xi1) <- rownames(Bhat)[grep("sao2_fio2_ratio_m", rownames(Bhat))]
        Xi2 <- model.matrix(fm.Xi, data = dtl2)
        Xi3 <- model.matrix(fm.Xi, data = dtl2)
        Xi <- as.matrix(bdiag(list(Xi1, Xi2, Xi3)))
        
        Yhat.o <-  Xi %*% Bhat 
        Yhati.o = matrix(Yhat.o, ncol = 1)
        
        Yi = Yi.o
        Yhati = Yhati.o
        Vi = Vi.o
        Vi.invsqrtm = Vi.invsqrtm.o #Vi^{-1/2}
        # Y ~ G(Yhat, Vi) 
        
        Yi.wt = matrix(M %*% Yi.o, ncol = 1)
        Yhati.wt = matrix(M %*% Yhati.o, ncol = 1)
        Vi.wt  = M %*% Vi.o %*% t(M)
        Vi.invsqrtm.wt = solve(sqrtm(Vi.wt))
        ni.wt <- length(Yi.wt)
        
        if(length(Yhati.o) >= nbiom & sum(is.na(Yi)) == 0){  ## at future days , all complete
          tmp.g <- (2*pi)^(-ni/2) * det(Vi)^(-1/2) * 
            exp((-1/2) * t(Yi - Yhati) %*% solve(Vi) %*% (Yi - Yhati)) ## decorR = t(Yi - Yhati) %*% invsqrtm(Vi) %*% (Yi - Yhati); kernel = decorR %*% decorR
          
          tmp.g.wt <- (2*pi)^(-ni.wt/2) * det(Vi.wt)^(-1/2) * 
            exp((-1/2) * t(Yi.wt - Yhati.wt) %*% solve(Vi.wt) %*% (Yi.wt - Yhati.wt)) ## decorR = t(Yi - Yhati) %*% invsqrtm(Vi) %*% (Yi - Yhati); kernel = decorR %*% decorR
          
          ## change to t-likelihood  
          # decorRes <- matrix(t(Yi - Yhati)  %*% Vi.invsqrtm,ncol = 1)
          decorRes <- matrix(t(Yi - Yhati)  %*% solve(Vi) %*% (Yi-Yhati),ncol = 1)
          const <- gamma((v.df+ni)/2)/((sqrt(v.df^ni * pi^ni)*gamma(v.df/2))*sqrt(det(Vi)))
          t.loglik <- log(const) - ((v.df+ni)/2)*log(1 + 1/v.df * decorRes)
          tmp <- exp(t.loglik)
          
          decorRes.wt <- matrix(t(Yi.wt - Yhati.wt)  %*% solve(Vi.wt) %*% (Yi.wt-Yhati.wt),ncol = 1)
          const.wt <- gamma((v+ni)/2)/((sqrt(v^ni * pi^ni)*gamma(v/2))*sqrt(det(Vi.wt)))
          t.loglik.wt <- log(const.wt) - ((v+ni)/2)*log(1 + 1/v * decorRes.wt)
          tmp.wt <- exp(t.loglik.wt)
          
        }else if(length(Yhati.o) >= nbiom & sum(is.na(Yi))!= 0){ ## at future days, some measure missing
          position <- which(!is.na(Yi.o))
          Vi.pos <- Vi.o[position, position]
          Yhati.pos <- Yhati.o[position]
          Yi.pos <- Yi.o[position]
          ni.pos <- length(position)
          Vi.invsqrtm.pos = Vi.invsqrtm.o[position, position]
          
          vals <- c(Yi.o[c((i-1),(i+1))],
                    Yi.o[c((2*(i+1)-2),(2*(i+1)))],
                    Yi.o[c((3*(i+1)-2),(3*(i+1)))])
          position.wt <- which(!is.na(vals))
          if(length(position.wt) == 2*nbiom){
            Yi.n <- Yi.o;Yi.n[is.na(Yi.n)] <- 0
            
            Yi.wt = matrix(M %*% Yi.n, ncol = 1)
            Yhati.wt = matrix(M %*% Yhati.o, ncol = 1)
            Vi.wt  = M %*% Vi.o %*% t(M)
            Vi.invsqrtm.wt = solve(sqrtm(Vi.wt))
            ni.wt <- length(Yi.wt)
            
            Vi.pos.wt <- Vi.wt
            Yhati.pos.wt <- Yhati.wt
            Yi.pos.wt <- Yi.wt
            ni.pos.wt <- length(Yi.pos.wt)
            Vi.invsqrtm.pos.wt = Vi.invsqrtm.wt
          }else{
            Vi.pos.wt <- Vi.wt
            Yhati.pos.wt <- Yhati.wt
            Yi.pos.wt <- Yi.wt
            ni.pos.wt <- length(Yi.pos.wt)
            message(s," NA in rencent 3 observaionts")
          }
          
          
          ## gaussian likelihood
          tmp.g <- (2*pi)^(-ni.pos/2) * det(Vi.pos)^(-1/2) * exp((-1/2) * t(Yi.pos - Yhati.pos) %*% solve(Vi.pos) %*% (Yi.pos - Yhati.pos))
          tmp.g.wt <- (2*pi)^(-ni.pos.wt/2) * det(Vi.pos.wt)^(-1/2) * exp((-1/2) * t(Yi.pos.wt - Yhati.pos.wt) %*% solve(Vi.pos.wt) %*% (Yi.pos.wt - Yhati.pos.wt)) ## decorR = t(Yi - Yhati) %*% invsqrtm(Vi) %*% (Yi - Yhati); kernel = decorR %*% decorR
          
          ## t likelihood
          decorRes <- matrix(t(Yi.pos - Yhati.pos)  %*% solve(Vi.pos) %*% (Yi.pos-Yhati.pos),ncol = 1)
          const <- gamma((v.df+ni)/2)/((sqrt(v.df^ni * pi^ni)*gamma(v.df/2))*sqrt(det(Vi.pos)))
          t.loglik <- log(const) - ((v.df+ni)/2)*log(1 + 1/v.df * decorRes)
          tmp <- exp(t.loglik)
          
          decorRes.wt <- matrix(t(Yi.pos.wt - Yhati.pos.wt)  %*% solve(Vi.pos.wt) %*% (Yi.pos.wt-Yhati.pos.wt),ncol = 1)
          const.wt <- gamma((v+ni)/2)/((sqrt(v^ni * pi^ni)*gamma(v/2))*sqrt(det(Vi.pos.wt)))
          t.loglik.wt <- log(const.wt) - ((v+ni)/2)*log(1 + 1/v * decorRes.wt)
          tmp.wt <- exp(t.loglik.wt)
          
        }
        lkh <- cbind(lkh, tmp)
        lkh.wt <- cbind(lkh.wt, tmp.wt)
        lkh.g <- cbind(lkh.g, tmp.g)
        lkh.g.wt <- cbind(lkh.g.wt, tmp.g.wt)
      }
      Lt <- rbind(Lt, lkh)
      colnames(Lt) <- c("L1", "L2", "L3")
      Lt.wt <- rbind(Lt.wt, lkh.wt)
      colnames(Lt.wt) <- c("L1", "L2", "L3")
      
      Lg <- rbind(Lg, lkh.g)
      colnames(Lg) <- c("L1", "L2", "L3")
      Lg.wt <- rbind(Lg.wt, lkh.g.wt)
      colnames(Lg.wt) <- c("L1", "L2", "L3")
    }
    L <- Lg.wt
    ## prior -----------------
    j <- (i+1)+projT
    # data <- dt.pre.mult %>% filter(osler_id == sid) 
    data.prior <- data %>% dplyr::select(osler_id,resp_bl, temp_bl, pulse_bl, sao2fio2_bl,
                                         crp_bl, alc_bl, dd_bl, gfr_bl,
                                         charlson_cat, demo5, sex, bmi_cat,smoking) %>% unique()
    
    newdt <- data.prior
    newdt <- do.call("rbind", replicate(projT, newdt, simplify = F))
    newdt$interval <- c((i+1):(j-1))#newdt %>% mutate(interval = seq(0,j))
    #newdt$rev_day <- c(-(j-1):-9)
    dtl<-newdt
    
    dtl <- cbind(dtl, predict(mult.pre, type = 'prob', newdata = dtl))
    cir <- dtl %>% mutate(cumprod = cumprod(`0`), 
                          cumsum1 = cumsum(`1`), 
                          cumsum2 = cumsum(`2`), 
                          cumsum3 = cumsum(`3`)) %>% 
      mutate(temp1 = `1`*`cumprod`/`0`,temp2 = `2`*`cumprod`/`0`, temp3 = `3`*`cumprod`/`0`,
             cir1 = cumsum(temp1), cir2 = cumsum(temp2), cir3 = cumsum(temp3)) 
    
    rem_p <- 1-sum(cir$temp1,cir$temp2,cir$temp3, na.rm=T)
    cir <- cir  %>% 
      dplyr::mutate(prior1 = temp1 + rem_p/(nbiom*(projT-1)), 
                    prior2 = temp2 + rem_p/(nbiom*(projT-1)),
                    prior3 = temp3 + rem_p/(nbiom*(projT-1)))
    xi <- cir %>% dplyr::select(osler_id, interval, prior1, prior2, prior3)
    
    ## bind Likelihood with prior --------------
    lkh <- L
    dt.comp <- cbind(xi, lkh)
    dt.comp <- dt.comp %>% mutate(tmp1 = prior1*L1,
                                  tmp2 = prior2*L2,
                                  tmp3 = prior3*L3) 
    
    denom <- sum(dt.comp$tmp1, dt.comp$tmp2,dt.comp$tmp3, na.rm = T)
    dt.comp <- dt.comp %>% 
      mutate(pprob1 = tmp1/denom, pprob2 = tmp2/denom, pprob3=  tmp3/denom)
    
    dtprob.save[[s]] <- dt.comp 
    ## end --------------
  }
  
  pprob.save <- list()
  for(s in 1:length(dtprob.save)){
    tmp <- dtprob.save[[s]]
    tmp <- tmp %>% dplyr::select(pprob1:pprob3)
    pprob.save[[s]] <- as.matrix(tmp)
  }
  pat_prob_save[[v]] <- pprob.save
  # pprob_d4[[v]] <- Reduce("+", pprob.save)[(i+2):21,]
  pprob_d4[[v]] <- Reduce("+", pprob.save)
  #observed table
  dt.endpoints.test <- dt.endpoints %>% 
    filter(osler_id %in% test_ids)
  obs.freq <- table(dt.endpoints.test$endpoint.u, dt.endpoints.test$status.u)
  
  if(ncol(obs.freq)!= 3){
    event.mis <- which(!(c(1,2,3) %in% colnames(obs.freq)))
    col.mis <- matrix(rep(0,nrow(obs.freq)), ncol = 1)
    colnames(col.mis) <- event.mis
    obs.freq <- cbind(obs.freq, col.mis)
    obs.freq <- obs.freq[order(as.numeric(rownames(obs.freq))),]
  }
  
  if(nrow(obs.freq)!= (maxT-i)){
    dayseq <- (i+1):maxT
    day.mis <- dayseq[which(!((i+1):maxT %in% (rownames(obs.freq))))]
    row.mis <- matrix(rep(c(0,0,0), length(day.mis)), nrow = length(day.mis))
    rownames(row.mis) <- day.mis
    obs.freq <- rbind(obs.freq, row.mis)
    obs.freq <- obs.freq[order(as.numeric(rownames(obs.freq))),]
  }
  
  
  
  exp.freq <- pprob_d4[[v]]
  
  exp.freq.cv[[v]]<-exp.freq
  obs.freq.cv[[v]]<-obs.freq
  pr.cv[[v]] <- (obs.freq - exp.freq)/sqrt(exp.freq)
}

exp.freq.save <- exp.freq.cv
obs.freq.save <- obs.freq.cv


pprob_d4_new <- obs.freq.new <- list()
id.remove <- NULL
## if no recent observation available
for(v in 1:5){
  test_ids <- c(test_ids1[[v]],test_ids2[[v]],test_ids3[[v]])
  test_ids_4d <- test_ids[test_ids %in% dt.endpoints$osler_id]
  exp.freq <-  exp.freq.cv[[v]]
  obs.freq <- obs.freq.cv[[v]]
  pprob.save <- pat_prob_save[[v]]
  if(sum(!is.na(exp.freq)) == 0){#if contains NA
    tmp <- lapply(pprob.save, is.na)
    tmp <- do.call(rbind, lapply(tmp, sum))
    sNA <- which(tmp != 0)
    pprob.save.new <- pat_prob_save[[v]]
    for(s in 1:length(sNA)){
      nid <- sNA[s]
      pprob.save.new[[nid]] <- 0

      ## remove observed values
      sid <- test_ids_4d[nid]
      sid.end <- dt.endpoints %>% filter(osler_id == sid)
      endpoint.u <- sid.end$endpoint.u
      status.u <- sid.end$status.u
      obs.freq[rownames(obs.freq) == endpoint.u, colnames(obs.freq) == status.u] <- 
        obs.freq[rownames(obs.freq) == endpoint.u, colnames(obs.freq) == status.u] - 1
      obs.freq.new[[v]] <- obs.freq
      
      id.remove <- rbind(id.remove, c(v, nid, sid))
    }
    pprob_d4_new[[v]] <- Reduce("+", pprob.save.new)
  }else{
    pprob_d4_new[[v]] <- exp.freq.cv[[v]]
    obs.freq.new[[v]] <- obs.freq.cv[[v]]
  }
}
dt.endpoints %>% filter(osler_id %in% id.remove[,3])
exp.freq.cv <- pprob_d4_new
obs.freq.cv <- obs.freq.new

dt.exp <- NULL
dt.pr <- NULL
dt.obs <- NULL
for(k in 1:5){
  day <- (i+1):20
  cv = rep(k, length(day))
  tmp.exp <- exp.freq.cv[[k]]; tmp.exp <- cbind(tmp.exp,day); tmp.exp <- cbind(tmp.exp,cv);
  tmp.obs <- obs.freq.cv[[k]]; tmp.obs <- cbind(tmp.obs,day); tmp.obs <- cbind(tmp.obs,cv);
  tmp.pr <- pr.cv[[k]]; tmp.pr <- cbind(tmp.pr, day); tmp.pr <- cbind(tmp.pr,cv);
  dt.exp <- rbind(dt.exp, tmp.exp)
  dt.obs <- rbind(dt.obs, tmp.obs)
  dt.pr <- rbind(dt.pr, tmp.pr)
}
dt.exp <- as.data.frame(dt.exp) %>% rename(Discharge = pprob1, Die = pprob2, Vent = pprob3)%>% gather(event, exp, 1:3)
dt.pr <- as.data.frame(dt.pr) %>% rename(Discharge = `1`, Die = `2`, Vent = `3`) %>% gather(event, pr, 1:3)
dt.obs <- as.data.frame(dt.obs) %>% rename(Discharge = `1`, Die = `2`, Vent = `3`) %>% gather(event, obs, 1:3)

dtp <- dt.exp %>% full_join(dt.pr) %>% full_join(dt.obs) 

p <- dtp %>% 
  dplyr::group_by(day, event) %>% 
  dplyr::summarize(exp.tot = sum(exp, na.rm = T),
                   obs.tot = sum(obs, na.rm = T)) %>% 
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

exp.freq <- Reduce("+", exp.freq.cv)
obs.freq <- Reduce("+", obs.freq.cv)
tmp <- sum((obs.freq-exp.freq)^2/exp.freq)
calib.dt <- tibble(`Chi-square stat` = round(tmp,2))

po <- p+geom_table_npc(data = calib.dt,
                       label = list(calib.dt),
                       npcx = 20, npcy = 1,
                       #table.theme = ttheme_gtlight(core = list(bg_params = list(fill = c("white")))),
                       table.theme = ttheme_gtminimal(),
                       size =4)

bkwr.save2 <- list(p, po,exp.freq.cv, obs.freq.cv, pat_prob_save)
save(bkwr.save2, 
     file = paste0("model/bkwr_di=", i,"_M.RData"))

save(id.remove, 
     file = paste0("model/id_remove_", i,".RData"))


