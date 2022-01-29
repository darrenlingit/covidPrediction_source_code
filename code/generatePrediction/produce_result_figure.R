## result figures updated
# 00 MCMCglmm for pre and post  ----------
rm(list=ls())
nbiom <- 3
df.re <- 3
setwd("S:/COVID_LDA/code/zwang_mgbowring/prediction_model")
source("00_loadlibs.R")
dt_joint <- readRDS("data/data_for_joint.rds")

# dtjimp <- dt_joint%>% dplyr::select(osler_id, resp_bl, temp_bl, pulse_bl, sao2fio2_bl, crp_bl, alc_bl,
#                                     dd_bl, gfr_bl,charlson_cat, demo5, sex, 
#                                     bmi_cat, smoking) %>% unique()
# dtjimp$smoking <- factor(dtjimp$smoking, levels = c("Never Smoker", "Former Smoker", "Current Smoker"))
# dtjimp$sex <- factor(dtjimp$sex, levels = c("Female", "Male"))
# 
# imp1 <- mice(dtjimp,m=10,set.seed = 2020)
# imp1$loggedEvents
# saveRDS(imp1, "paperPrep/imp_dtjoint.rds")
imp1 <- readRDS("paperPrep/imp_dtjoint.rds")
dtj <- complete(imp1, action=1) #take the first imputed dataset

dt.joint <- dt_joint %>% dplyr::select(osler_id, interval, rev_day, status, sao2_fio2_ratio_m, pulse_m, temp_c_m, outcome)%>% 
  left_join(dtj) %>% 
  filter(complete.cases(dd_bl, crp_bl, alc_bl, gfr_bl, 
                        resp_bl, sao2fio2_bl, pulse_bl, temp_bl, charlson_cat, demo5, 
                        sex, bmi_cat, smoking, rev_day, status))

nb = 2000; nt = 8000
## filter for pre-ventilation
dt.pre <- dt.joint %>% filter(status %in% c(1,2,3)) #length(unique(dt.pre$osler_id)) #1762 patients
dt.pre$long <- ifelse(dt.pre$rev_day < -20,'long','short')
id.check.pre <- dt.pre %>% dplyr::select(osler_id, long, status) %>% filter(long == "long") %>% unique()#28 patients
dt.pre <- dt.pre %>% dplyr::filter(!(osler_id %in% id.check.pre$osler_id)) ## 1734 patients

dt.pre <- dt.pre %>% 
  dplyr::rename(sao2_fio2_ratio_o = sao2_fio2_ratio_m,
                temp_c_o = temp_c_m,
                pulse_o = pulse_m)

dt.pre <- dt.pre  %>% mutate(sex = factor(sex, levels = c("Female", "Male")),
                             smoking = factor(smoking, levels = c("Never Smoker", "Former Smoker", "Current Smoker")),
                             charlson_cat = factor(charlson_cat, levels = c("0", "1-2" , "3-4", ">5")),
                             demo5 = factor(demo5, levels =   c("Young_White", "Old_White" , "Young_Black", "Old_Black", "Latinx/Other")),
                             bmi_cat = factor(bmi_cat, levels = c("<30", ">=30")))
library(bestNormalize)
dt.pre <- dt.pre %>% 
  filter(complete.cases(sao2_fio2_ratio_o, temp_c_o, pulse_o)) %>% 
  mutate(sao2_norm = orderNorm(sao2_fio2_ratio_o)$x.t,
         temp_norm = orderNorm(temp_c_o)$x.t,
         pulse_norm = orderNorm(pulse_o)$x.t)

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

## likelihood --------------
maxT <- 20
#sid <- "00f8b5b9-7462-4738-a626-588d697345bb"  ## not use
#sid <- "d05319ea-c7a4-4979-b399-27bc4f843b4d"  ## day10
sid <- "5a9cb1e1-d2f3-4bbb-b25e-7acb970727e7"  ## day 12
#sid <- "7899a2c7-47ef-4cdb-b4d4-62ac484901cb"   ## day 4 or 5

data <- dt.pre %>% filter(osler_id == sid) 
load("paperPrep/finalizedCode/model/alldata_jtgaussian_model.RData")
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

L <- list()
for(i in 0:max(data$interval)){ #for(i in 0:unique(data$endpoint))
  Lg <- NULL
  Lg.wt <- NULL
  projT <- maxT-i
  for(j in 1:projT){
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
  
    ## create weight
    if(i<2){
      M <- diag(1, nrow = length(Yi.o))
    }else{
      row1 <- c(rep(0,(i+1)-1), 1, rep(0, (nbiom-1)*(i+1)))
      row2 <- c(rep(0,(nbiom-1)*(i+1)-1), 1, rep(0, (i+1)))
      row3 <- c(rep(0,nbiom*(i+1)-1),1)
      row4 <- c(rep(0, (i+1)-3), -1/2, 0, 1/2, rep(0,  (nbiom-1)*(i+1)))
      row5 <- c(rep(0,(i+1)), rep(0, (i+1)-3), -1/2, 0, 1/2, rep(0,(i+1)))
      row6 <- c(rep(0,(nbiom-1)*(i+1)), rep(0, (i+1)-3), -1/2, 0, 1/2)
      
      M <- rbind(row1, row2, row3, row4, row5, row6)
    }
    
    for(k in 1:3){
      # inner knots: -2, -5, -10,-15(wherevern ventilation ends)
      dtl2 <- dtl %>% mutate(status = factor(k , levels = c("1", "2", "3")))
      # fm.Xi <- as.formula(paste0("~status * ns(rev_day,knots = c(-7,-4,-2,-1),
      #                          Boundary.knots = c(-20,0)) + resp_bl + temp_bl +
      #                          pulse_bl + sao2fio2_bl + crp_bl + alc_bl + dd_bl + gfr_bl +
      #                          charlson_cat + demo5 + sex + bmi_cat + smoking"))
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
      
      ## U = M'Y ~ G(M'mu_i, M'V_iM)
      ## M n_i x 6, intercept and linear slope over past 3 days for 3 biomarkers
      
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
      }
      lkh.g <- cbind(lkh.g, tmp.g)
      lkh.g.wt <- cbind(lkh.g.wt, tmp.g.wt)
    }
    Lg <- rbind(Lg, lkh.g)
    colnames(Lg) <- c("L1", "L2", "L3")
    Lg.wt <- rbind(Lg.wt, lkh.g.wt)
    colnames(Lg.wt) <- c("L1", "L2", "L3")
  }
  L[[i+1]] <- Lg.wt
}

# prior and posterior probabilities ------
load("paperPrep/finalizedCode/model/alldata_multinom_model.RData")
dt1c <- readRDS("data/dt_pre.rds")
dt1c <- dt1c %>% filter(!is.na(endpoint)) %>% filter(endpoint != 0) ## drop patients AC or had an event on day 0

imp1 <- readRDS("model/mice_imp_pre_1029.rds") ## imputation 
dt1imp <- complete(imp1, action=1) #take the first imputed dataset

dt1c <- dt1c %>% dplyr::select(osler_id, interval, outcome2, status, first_vented, last_interval, group2, 
                               rev_day, daily_who, first_vented_interval, endpoint, Wt) %>% left_join(dt1imp, by = 'osler_id')
dt1c.sid <- dt1c %>% filter(osler_id == sid)
data <- data %>% left_join(dt1c.sid %>% dplyr::select(osler_id, interval, status, first_vented, last_interval,
                                                      group2, first_vented_interval, endpoint, Wt, maxwho))
prior <- list()
for(i in 0:unique(data$endpoint)){
  r <- i+1
  prt <- NULL
  projT <- maxT-i
  j <- r+projT

    baseline <- data %>% 
      dplyr::select(osler_id,resp_bl, temp_bl, pulse_bl, sao2fio2_bl,
                    crp_bl, alc_bl, dd_bl, gfr_bl, charlson_cat, demo5, sex, bmi_cat,smoking) %>% 
      unique()
    total.length <- maxT  - i
    newdt <- do.call("rbind", replicate(total.length, baseline, simplify = F))
    newdt$interval <- c((i+1):(maxT))#newdt %>% mutate(interval = seq(0,j))
    dtl <- newdt

  
  select =  dtl$interval < j
  dtl <- dtl[select,]
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
  
  prior[[r]] <- cir
}

dt.full1 <- list()
for(i in 1:length(L)){
  d <- i-1
  dt.lkh <- NULL
  lkh <- L[[i]]
  xi <- prior[[i]]
  projT <- maxT-d
  dt.comp <- cbind(xi, lkh) 
  dt.comp <- dt.comp %>% mutate(tmp1 = prior1*L1,
                                tmp2 = prior2*L2,
                                tmp3 = prior3*L3) 
  
  denom <- sum(dt.comp$tmp1, dt.comp$tmp2,dt.comp$tmp3, na.rm = T)
  dt.comp <- dt.comp %>% 
    mutate(pprob1 = tmp1/denom, pprob2 = tmp2/denom, pprob3=  tmp3/denom)
  dt.full1[[i]] <- dt.comp
}

### plot prior, likelihood and posterior ---------
library(viridis)
library(RColorBrewer)
library(scales)
fill.cols <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
int_breaks <- function(n=5,...){
  fxn <- function(x) {
    breaks <- floor(pretty(x,n,...))
    names(breaks) <- attr(breaks ,"labels")
    breaks
  }
  return(fxn)
}

dt.pprior <- list()
ran <- NULL
for(i in 1:length(prior)){
  dt.p <- prior[[i]]
  dt.p <- dt.p 
  dt.p <- dt.p%>% mutate(prior1.norm = log(prior1) - max(log(dt.p$prior1),log(dt.p$prior2),log(dt.p$prior3)),
                         prior2.norm = log(prior2) - max(log(dt.p$prior1),log(dt.p$prior2),log(dt.p$prior3)),
                         prior3.norm = log(prior3) - max(log(dt.p$prior1),log(dt.p$prior2),log(dt.p$prior3)))
  dt.p <- dt.p %>% gather(Ltype, prior, c(prior1.norm:prior3.norm))
  ran <- rbind(ran, range(dt.p$prior, na.rm = T))
  dt.pprior[[i]] <- dt.p
}
ran.prior <- range(ran)

dt.ppost <- list()
ran <- NULL
options(warn = 1)
for(i in 1:length(dt.full1)){
  dt.p <- dt.full1[[i]]
  dt.p <- dt.p%>% mutate(post1.norm = log(pprob1) - max(log(dt.p$pprob1),log(dt.p$pprob2),log(dt.p$pprob3), na.rm = T),
                         post2.norm = log(pprob2) - max(log(dt.p$pprob1),log(dt.p$pprob2),log(dt.p$pprob3), na.rm = T),
                         post3.norm = log(pprob3) - max(log(dt.p$pprob1),log(dt.p$pprob2),log(dt.p$pprob3), na.rm = T))
  dt.p <- dt.p %>% gather(Ltype, posterior, c(post1.norm:post3.norm))
  ran <- rbind(ran, range(dt.p$posterior, na.rm = T))
  dt.ppost[[i]] <- dt.p
}
ran.post <- range(ran)

ran <- NULL
dt.plik <- list()
for(i in 1:length(L)){
  dt.p <- dt.full1[[i]]
  dt.p <- dt.p %>% dplyr::select(interval, L1, L2, L3)
  dt.p <- dt.p%>% mutate(ll1.CB = log(L1) - max(log(dt.p$L1),log(dt.p$L2),log(dt.p$L3)),
                         ll2.CB = log(L2) - max(log(dt.p$L1),log(dt.p$L2),log(dt.p$L3)),
                         ll3.CB = log(L3) - max(log(dt.p$L1),log(dt.p$L2),log(dt.p$L3)))
  dt.p <- dt.p %>% gather(Ltype, llt, c(ll1.CB:ll3.CB))
  ran <- rbind(ran, range(dt.p$llt, na.rm = T))
  dt.plik[[i]] <- dt.p
}
ran.lik <- range(ran)
ran <- range(ran.post, ran.lik, ran.prior)

pr <- list()
lim.x <- c(0,maxT)
for(i in 1:length(prior)){
  dt.p <- dt.pprior[[i]]
  pr[[i]] <-  ggplot(data = dt.p,aes(y = Ltype,  x = interval, fill = prior))+
    geom_tile()+
    theme_minimal()+
    scale_fill_viridis(name = "",
                       limits = ran, direction = -1)+
    scale_x_continuous(name = "",limits =c(0, maxT))+
    ggtitle(paste0("Day",(i)))+
    ylab("")+
    scale_y_discrete(limits = c("prior3.norm", "prior2.norm", "prior1.norm"),
                     labels = c("prior1.norm" = "Discharge", "prior2.norm" = "Die", "prior3.norm" = "Vent"))+
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8))
}


pp <- list()
lim.x <- c(0,maxT)
for(i in 1:length(dt.full1)){
  dt.p <- dt.ppost[[i]]
  pp[[i]] <-  ggplot(data = dt.p,aes(y = Ltype,  x = interval, fill = posterior))+
    geom_tile()+
    theme_minimal()+
    scale_fill_viridis(name = "",
                       limits = ran, direction = -1)+
    scale_x_continuous(name = "",limits =c(0, maxT))+
    ggtitle(paste0("Day",(i)))+
    ylab("")+
    scale_y_discrete(limits = c("post3.norm", "post2.norm", "post1.norm"),
                     labels = c("post1.norm" = "Discharge", "post2.norm" = "Die", "post3.norm" = "Vent"))+
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8))
}

pl <-list()
for(i in 1:length(L)){
  dt.p <- dt.plik[[i]]
  pl[[i]] <-  ggplot(data = dt.p,aes(y = Ltype,  x = interval, fill = llt))+
    geom_tile()+
    theme_minimal()+
    scale_fill_viridis(name = "",
                       limits = ran, direction = -1)+
    scale_x_continuous( name = "",limits =c(0, maxT))+
    ggtitle(paste0("Day",(i)))+
    ylab("")+
    scale_y_discrete(limits = c("ll3.CB", "ll2.CB", "ll1.CB"),
                     labels = c("ll1.CB" = "Discharge", "ll2.CB" = "Die", "ll3.CB" = "Vent"))+
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8))
  
}

## conditional mean and result ---------
joint.mean <- list()
joint.var <- list()
for(i in 0:unique(data$endpoint)){
  dtc <- data
  last.intv <- unique(dtc$last_interval)
  info.pat <- dtc %>% dplyr::select(-interval, -rev_day, -Wt, -sao2_norm, -temp_norm, -pulse_norm,
                                    -sao2_fio2_ratio_o, -pulse_o, -temp_c_o, -temp_norm, -pulse_norm) %>% unique()
  if(nrow(dtc) != (last.intv + 1)){
    fk.length <- last.intv + 1
    fk.dtc <- data.frame(osler_id = rep(unique(dtc$osler_id), fk.length),
                         interval = 0:last.intv)
    fk.dtc <- fk.dtc %>% left_join(info.pat, by = c('osler_id'))
    dtc <- dtc %>% full_join(fk.dtc) %>% arrange(interval)
  }
  cM <- list()
  cV <- list()
  projT <- maxT - i
  for(j in 1:projT){
    cM.tmp <- list()
    select =   dtc$interval <= i+j
    dtl <- dtc[select,]
    dtl <- dtl %>% mutate(endpoint = i+j, rev_day = interval - endpoint)
    if(max(dtl$rev_day)!=0){
      dtl2 <- rbind(dtl, dtl[1:(0-max(dtl$rev_day)),])
      dtl2<- dtl2 %>%  
        mutate(rev_day = c(-(i+j):0))
      dtl2$gfr_bl <- na.locf(dtl2$gfr_bl)
      dtl2$alc_bl <- na.locf(dtl2$alc_bl)
      dtl2$crp_bl <- na.locf(dtl2$crp_bl)
      dtl2$dd_bl <- na.locf(dtl2$dd_bl)
      dtl2$pulse_bl <- na.locf(dtl2$pulse_bl)
      dtl2$temp_bl <- na.locf(dtl2$temp_bl)
      dtl2$resp_bl <- na.locf(dtl2$resp_bl)
      dtl2$sao2fio2_bl <- na.locf(dtl2$sao2fio2_bl)
      
      dtl2$charlson_cat <- na.locf(dtl2$charlson_cat)
      dtl2$sex <- na.locf(dtl2$sex)
      dtl2$demo5 <- na.locf(dtl2$demo5)
      dtl2$bmi_cat <- na.locf(dtl2$bmi_cat)
      dtl2$smoking <- na.locf(dtl2$smoking)
      
      Zi1 =  model.matrix(~1+ns(rev_day, Boundary.knots = attr(re.spline,"Boundary.knots"),
                                knots = attr(re.spline,"knots")), dtl2)
      Zi2 = model.matrix(~1+ns(rev_day, Boundary.knots = attr(re.spline,"Boundary.knots"),
                               knots = attr(re.spline,"knots")), dtl2)
      Zi3 = model.matrix(~1+ns(rev_day, Boundary.knots = attr(re.spline,"Boundary.knots"),
                               knots = attr(re.spline,"knots")), dtl2)
      Zi <- as.matrix(bdiag(list(Zi1, Zi2, Zi3)))
      ni <- (i+1)*nbiom    #ni = n or nk?
      
      ## process Ri 0120
      Rcov <- smy[["Rcovariances"]][,"post.mean"]
      rn.R <- names(Rcov)
      R = matrix(Rcov, nrow = nbiom)
      colnames(R) <- rownames(R) <-  rbind(lapply(str_split(rn.R[1:(nbiom)],":trait"), `[[`,1))
      vars <- c('sao2_fio2_ratio_m', 'temp_c_m', 'pulse_m')
      days <- c(0:(i+j))
      nms.R <- do.call(paste0, expand.grid(vars, days))
      Ri <- as.matrix(bdiag(rep(list(R), i+j+1)))
      colnames(Ri) <- rownames(Ri) <- nms.R
      ord.Ri <- c(grep("sao2_fio2_ratio_m", colnames(Ri)),
                  grep("temp_c_m", colnames(Ri)),
                  grep("pulse_m", colnames(Ri)))
      Ri <- Ri[ord.Ri, ord.Ri]
      
      Yi <- matrix(c(dtl$sao2_norm, 
                     dtl$temp_norm, 
                     dtl$pulse_norm), ncol = 1)
      Vi = Zi %*% D %*% t(Zi) + Ri   ## future improvement: D variance matrix could be different based on status 0106
    }else{
      Zi1 =  model.matrix(~1+ns(rev_day, Boundary.knots = attr(re.spline,"Boundary.knots"),
                                knots = attr(re.spline,"knots")), dtl)
      Zi2 = model.matrix(~1+ns(rev_day, Boundary.knots = attr(re.spline,"Boundary.knots"),
                               knots = attr(re.spline,"knots")), dtl)
      Zi3 = model.matrix(~1+ns(rev_day, Boundary.knots = attr(re.spline,"Boundary.knots"),
                               knots = attr(re.spline,"knots")), dtl)
      Zi <- as.matrix(bdiag(list(Zi1, Zi2, Zi3)))
      ni <- (i+1)*nbiom  #ni = n or nk?
      
      ## process Ri 0120
      Rcov <- smy[["Rcovariances"]][,"post.mean"]
      rn.R <- names(Rcov)
      R = matrix(Rcov, nrow = nbiom)
      colnames(R) <- rownames(R) <-  rbind(lapply(str_split(rn.R[1:(nbiom)],":trait"), `[[`,1))
      vars <- c('sao2_fio2_ratio_m', 'temp_c_m', 'pulse_m')
      days <- c(0:(i+j))
      nms.R <- do.call(paste0, expand.grid(vars, days))
      Ri <- as.matrix(bdiag(rep(list(R), i+j+1)))
      colnames(Ri) <- rownames(Ri) <- nms.R
      ord.Ri <- c(grep("sao2_fio2_ratio_m", colnames(Ri)),
                  grep("temp_c_m", colnames(Ri)),
                  grep("pulse_m", colnames(Ri)))
      Ri <- Ri[ord.Ri, ord.Ri]
      
      Yi <- matrix(c(dtl$sao2_norm, 
                     dtl$temp_norm, 
                     dtl$pulse_norm), ncol = 1)
      Vi = Zi %*% D %*% t(Zi) + Ri 
    }
    for(k in 1:3){
      fm.Xi <- as.formula(paste0("~status * ns(rev_day,knots = c(-7,-4,-2,-1),
                                 Boundary.knots = c(-20,0)) + resp_bl + temp_bl +
                                 pulse_bl + sao2fio2_bl + crp_bl + alc_bl + dd_bl + gfr_bl +
                                 charlson_cat + demo5 + sex + bmi_cat + smoking"))
      if(max(dtl$rev_day)!=0){
        dtl2<- dtl2 %>% 
          mutate(status = factor(k , levels = c("1", "2", "3")))
        Xi1 <- model.matrix(fm.Xi, data = dtl2) #;colnames(Xi1) <- rownames(Bhat)[grep("sao2_fio2_ratio_m", rownames(Bhat))]
        Xi2 <- model.matrix(fm.Xi, data = dtl2) #;colnames(Xi2) <- rownames(Bhat)[grep("resp_rate_m", rownames(Bhat))]
        Xi3 <- model.matrix(fm.Xi, data = dtl2)
        X <- as.matrix(bdiag(list(Xi1, Xi2, Xi3)))
        
      }else{
        Xi1 <- model.matrix(fm.Xi, data = dtl %>% mutate(status = factor(k , levels = c("1", "2", "3")))) #;colnames(Xi1) <- rownames(Bhat)[grep("sao2_fio2_ratio_m", rownames(Bhat))]
        Xi2 <- model.matrix(fm.Xi, data = dtl %>% mutate(status = factor(k , levels = c("1", "2", "3")))) #;colnames(Xi2) <- rownames(Bhat)[grep("resp_rate_m", rownames(Bhat))]
        Xi3 <- model.matrix(fm.Xi, data = dtl %>% mutate(status = factor(k , levels = c("1", "2", "3")))) #;colnames(Xi2) <- rownames(Bhat)[grep("resp_rate_m", rownames(Bhat))]
        X <- as.matrix(bdiag(list(Xi1, Xi2, Xi3)))
      }
      
      Yhat <-  matrix(X %*% Bhat , ncol = 1)
      cM.tmp[[k]] <- Yhat
    }
    cV[[j]] <- Vi
    cM[[j]] <- cM.tmp
  }
  joint.mean[[i+1]] <- cM
  joint.var[[i+1]] <- cV
}
cond.mean <- list()
cond.var <-list()
for(i in 0:unique(data$endpoint)){
  r <- i+1
  dtc <- data
  last.intv <- unique(dtc$last_interval)
  if(nrow(dtc) != (last.intv + 1)){
    fk.length <- last.intv + 1
    fk.dtc <- data.frame(osler_id = rep(unique(dtc$osler_id), fk.length),
                         interval = 0:last.intv)
    dtc <- dtc %>% full_join(fk.dtc, by = c('osler_id', 'interval')) %>% arrange(interval)
  }
  cV <- cM <- list()
  projT <- maxT - i
  for(j in 1:projT){
    cM.tmp <- cV.tmp <- list()
    for(k in 1:3){
      Yit <- matrix(c(dtc$sao2_norm[1:r], 
                      dtc$temp_norm[1:r],
                      dtc$pulse_norm[1:r]), ncol = 1)
      j.mean <- joint.mean[[r]][[j]][[k]]
      j.var <- joint.var[[r]][[j]]
      
      nms <- c(paste0('sao2_fio2_ratio_m', '_day',c(0:(i+j))),
               paste0('temp_c_m', '_day',c(0:(i+j))),
               paste0('pulse_m', '_day',c(0:(i+j))))
      rownames(j.mean) <- nms
      colnames(j.var) <-nms
      
      toMatch.Yit <- paste0('day', c(0:i), "$")
      inx.Yit <- grep(paste(toMatch.Yit, collapse = "|"),nms)
      toMatch.Yit1ts <- paste0('day', c((i+1):(i+j)),  "$")
      inx.Yit1ts <- grep(paste(toMatch.Yit1ts, collapse = "|"),nms)
      mu.Yit = matrix(j.mean[inx.Yit], ncol = 1)
      mu.Yit1ts = matrix(j.mean[inx.Yit1ts], ncol = 1)
      var.Yit = j.var[inx.Yit,inx.Yit]
      var.Yit1ts = j.var[inx.Yit1ts, inx.Yit1ts]
      cov.Yit1ts.Yit = j.var[inx.Yit1ts,inx.Yit ]
      cov.Yit.Yit1ts =  t(cov.Yit1ts.Yit)
      
      if(sum(is.na(Yit)) == 0){  ## measures complete at all days
        cM.tmp[[k]] <- mu.Yit1ts + cov.Yit1ts.Yit %*% solve(var.Yit) %*% (Yit - mu.Yit)
        cV.tmp[[k]] <- var.Yit1ts - cov.Yit1ts.Yit %*% solve(var.Yit) %*% cov.Yit.Yit1ts
      }else if(sum(is.na(Yit)) > 0){  ## some days have some missings
        ind.c <- which(is.na(Yit))
        Yit.c <- matrix(Yit[-ind.c],ncol=1)
        mu.Yit.c <- matrix(mu.Yit[-ind.c], ncol = 1)
        var.Yit.c <- var.Yit[-ind.c, -ind.c]
        cov.Yit1ts.Yit.c <- matrix(cov.Yit1ts.Yit[-ind.c, -ind.c], nrow = j) ## this line is wrong
        cov.Yit.Yit1ts.c =  t(cov.Yit1ts.Yit.c)
        
        mu.Yit1ts.c <- mu.Yit1ts
        var.Yit1ts.c <- var.Yit1ts
        
        cM.tmp[[k]] <- mu.Yit1ts.c + cov.Yit1ts.Yit.c %*% solve(var.Yit.c) %*% (Yit.c - mu.Yit.c)
        cV.tmp[[k]] <- var.Yit1ts.c - cov.Yit1ts.Yit.c %*% solve(var.Yit.c) %*% cov.Yit.Yit1ts.c
      }
      
      rownames(cM.tmp[[k]]) <- nms[inx.Yit1ts]
    }
    cV[[j]] <- cV.tmp
    cM[[j]] <- cM.tmp
  }
  cond.mean[[i+1]] <- cM
  cond.var[[i+1]] <- cV
}

## plot patient result figure --------
fontsize <- 18
pred_figure <- pred_table <- list()
library(plyr) 
dt.ave<-list()
p.all <- list()
d0 <- data
options(warn =1)

load("paperPrep/finalizedCode/model/normobj.RData")
normobj.sao2 <- normobj[[1]]
normobj.temp <- normobj[[2]]
normobj.pulse <- normobj[[3]]
library(grid) 
library(cowplot)
library(formattable)
for(s in 1:length(cond.mean)){
  tmp.ave <-list()
  dtc <- d0
  cmean1 <- cond.mean[[s]] 
  dtfull <- dt.full1[[s]]
  last.intv <- unique(dtc$last_interval)
  if(nrow(dtc) != (last.intv + 1)){
    fk.length <- last.intv + 1
    fk.dtc <- data.frame(osler_id = rep(unique(dtc$osler_id), fk.length),
                         interval = 0:last.intv)
    dtc <- dtc %>% full_join(fk.dtc, by = c('osler_id', 'interval')) %>% arrange(interval)
  }
  
  nms <- rownames(cmean1[[1]][[1]])
  inx.biom1 <- grep('sao2_fio2_ratio_m', nms)
  inx.biom3 <- grep('temp_c_m', nms)
  inx.biom4 <- grep('pulse_m', nms)
  if(nrow(cmean1[[1]][[1]])/nbiom+s <= nrow(dtc)){
    leng.t <-length(dtc$sao2_norm[1:s])
    leng.b <- nrow(dtc) - length(cmean1[[1]][[1]][inx.biom1])-length(dtc$sao2_norm[1:s]) 
    dtc$cmeank1.biom1 <- c(rep(NA,leng.t), cmean1[[1]][[1]][inx.biom1], rep(NA, leng.b))
    dtc$cmeank2.biom1 <- c(rep(NA,leng.t), cmean1[[1]][[2]][inx.biom1], rep(NA, leng.b))
    dtc$cmeank3.biom1 <-  c(rep(NA,leng.t), cmean1[[1]][[3]][inx.biom1], rep(NA, leng.b))
    
    dtc$cmeank1.biom3 <- c(rep(NA,leng.t), cmean1[[1]][[1]][inx.biom3], rep(NA, leng.b))
    dtc$cmeank2.biom3 <- c(rep(NA,leng.t), cmean1[[1]][[2]][inx.biom3], rep(NA, leng.b))
    dtc$cmeank3.biom3 <-  c(rep(NA,leng.t), cmean1[[1]][[3]][inx.biom3], rep(NA, leng.b))
    
    dtc$cmeank1.biom4 <- c(rep(NA,leng.t), cmean1[[1]][[1]][inx.biom4], rep(NA, leng.b))
    dtc$cmeank2.biom4 <- c(rep(NA,leng.t), cmean1[[1]][[2]][inx.biom4], rep(NA, leng.b))
    dtc$cmeank3.biom4 <-  c(rep(NA,leng.t), cmean1[[1]][[3]][inx.biom4], rep(NA, leng.b))
    
    full1 <- dtfull %>% dplyr::select(osler_id, interval, pprob1, pprob2, pprob3)
    full1 <- full1[full1$interval == s,]
    dtc <- dtc %>% left_join(full1 %>% dplyr::select(osler_id, pprob1, pprob2, pprob3), by = c('osler_id'))
  }else{
    leng.t <-length(dtc$sao2_norm[1:s])
    cmeank1.biom1 <- c(rep(NA,leng.t), cmean1[[1]][[1]][inx.biom1])
    cmeank2.biom1 <- c(rep(NA,leng.t), cmean1[[1]][[2]][inx.biom1])
    cmeank3.biom1 <- c(rep(NA,leng.t), cmean1[[1]][[3]][inx.biom1])
    
    cmeank1.biom3 <- c(rep(NA,leng.t), cmean1[[1]][[1]][inx.biom3])
    cmeank2.biom3 <- c(rep(NA,leng.t), cmean1[[1]][[2]][inx.biom3])
    cmeank3.biom3 <- c(rep(NA,leng.t), cmean1[[1]][[3]][inx.biom3])
    
    cmeank1.biom4 <- c(rep(NA,leng.t), cmean1[[1]][[1]][inx.biom4])
    cmeank2.biom4 <- c(rep(NA,leng.t), cmean1[[1]][[2]][inx.biom4])
    cmeank3.biom4 <- c(rep(NA,leng.t), cmean1[[1]][[3]][inx.biom4])
    
    dt.ext <- data.frame(cmeank1.biom1, cmeank2.biom1, cmeank3.biom1, 
                         cmeank1.biom3,cmeank2.biom3,cmeank3.biom3,
                         cmeank1.biom4,cmeank2.biom4,cmeank3.biom4)
    dt.ext$interval <- 0:(1+s-1)
    dt.ext$osler_id <- sid
    dtc <- data %>% full_join(dt.ext, by =c('osler_id', "interval"))
    full1 <- dtfull %>% dplyr::select(osler_id, interval, pprob1, pprob2, pprob3)
    full1 <- full1[full1$interval == s-1+1,]
    dtc <- dtc %>% left_join(full1 %>% dplyr::select(osler_id, pprob1, pprob2, pprob3), by = 'osler_id')
    
  }
  tmp.ave[[1]] <- dtc
  
  dtc.tmp <- dtc%>% filter(interval <= s) %>% dplyr::rename(obs.biom1 = sao2_norm,
                                                            obs.biom3 = temp_norm,
                                                            obs.biom4 = pulse_norm)
  dtc.tmp$obs.biom1.inv = predict(normobj.sao2, newdata = dtc.tmp$obs.biom1, inverse = T)
  dtc.tmp$obs.biom3.inv = predict(normobj.temp, newdata = dtc.tmp$obs.biom3, inverse = T)
  dtc.tmp$obs.biom4.inv = predict(normobj.pulse, newdata = dtc.tmp$obs.biom4, inverse = T)
  
  dtc.tmp$cmeank1.biom1.inv = predict(normobj.sao2, newdata = dtc.tmp$cmeank1.biom1, inverse = T)
  dtc.tmp$cmeank2.biom1.inv = predict(normobj.sao2, newdata = dtc.tmp$cmeank2.biom1, inverse = T)
  dtc.tmp$cmeank3.biom1.inv = predict(normobj.sao2, newdata = dtc.tmp$cmeank3.biom1, inverse = T)
  
  dtc.tmp$cmeank1.biom3.inv = predict(normobj.temp, newdata = dtc.tmp$cmeank1.biom3, inverse = T)
  dtc.tmp$cmeank2.biom3.inv = predict(normobj.temp, newdata = dtc.tmp$cmeank2.biom3, inverse = T)
  dtc.tmp$cmeank3.biom3.inv = predict(normobj.temp, newdata = dtc.tmp$cmeank3.biom3, inverse = T)
  
  dtc.tmp$cmeank1.biom4.inv = predict(normobj.pulse, newdata = dtc.tmp$cmeank1.biom4, inverse = T)
  dtc.tmp$cmeank2.biom4.inv = predict(normobj.pulse, newdata = dtc.tmp$cmeank2.biom4, inverse = T)
  dtc.tmp$cmeank3.biom4.inv = predict(normobj.pulse, newdata = dtc.tmp$cmeank3.biom4, inverse = T)
  
  dtc <- dtc.tmp %>% 
    gather(key, cmean, c(obs.biom1.inv:obs.biom4.inv, cmeank1.biom1.inv:cmeank3.biom4.inv))%>% 
    separate(key, c('type', 'biom','inv'), "\\.") %>% 
    spread(type, cmean)
  
  biom.labeller <- c("SpO2/FiO2", "Temperature", "Pulse")
  names(biom.labeller) <- c("biom1", "biom3", "biom4")
  p0 <- ggplot()+
    geom_point(data = dtc[dtc$interval <= s-1,], aes(x = interval, y = obs), color = 'black')+
    geom_line(data = dtc[dtc$interval <= s-1,], aes(x = interval,y = obs), color = 'black')+
    geom_point(data = dtc, aes(x = interval,y = cmeank1, color = "Discharge", size = pprob1), alpha = 0.6)+
    geom_point(data = dtc, aes(x = interval,y = cmeank2, color = "Died", size = pprob2), alpha = 0.6)+
    geom_point(data = dtc, aes(x = interval,y = cmeank3, color = "Vent", size = pprob3), alpha = 0.6)+
    facet_grid(biom~., scales = "free_y", labeller = labeller(biom = biom.labeller))
  
  p <- p0
  
  for(i in 2:length(cmean1)){ #s=2: currently at day s-1 = 1; i=13: project next 13 days
    nms <- rownames(cmean1[[i]][[1]])
    inx.biom1 <- grep('sao2_fio2_ratio_m', nms)
    inx.biom3 <- grep('temp_c_m', nms)
    inx.biom4 <- grep('pulse_m', nms)
    if(nrow(cmean1[[i]][[1]])/nbiom+s <= nrow(data)){
      dtc <- data
      leng.t <-length(dtc$sao2_norm[1:s])
      leng.b <- nrow(dtc) - length(cmean1[[i]][[1]][inx.biom1])-length(dtc$sao2_norm[1:s])
      dtc$cmeank1.biom1 <- c(rep(NA,leng.t), cmean1[[i]][[1]][inx.biom1], rep(NA, leng.b))
      dtc$cmeank2.biom1 <- c(rep(NA,leng.t), cmean1[[i]][[2]][inx.biom1], rep(NA, leng.b))
      dtc$cmeank3.biom1 <-  c(rep(NA,leng.t), cmean1[[i]][[3]][inx.biom1], rep(NA, leng.b))
      
      dtc$cmeank1.biom3 <- c(rep(NA,leng.t), cmean1[[i]][[1]][inx.biom3], rep(NA, leng.b))
      dtc$cmeank2.biom3 <- c(rep(NA,leng.t), cmean1[[i]][[2]][inx.biom3], rep(NA, leng.b))
      dtc$cmeank3.biom3 <-  c(rep(NA,leng.t), cmean1[[i]][[3]][inx.biom3], rep(NA, leng.b))
      
      dtc$cmeank1.biom4 <- c(rep(NA,leng.t), cmean1[[i]][[1]][inx.biom4], rep(NA, leng.b))
      dtc$cmeank2.biom4 <- c(rep(NA,leng.t), cmean1[[i]][[2]][inx.biom4], rep(NA, leng.b))
      dtc$cmeank3.biom4 <-  c(rep(NA,leng.t), cmean1[[i]][[3]][inx.biom4], rep(NA, leng.b))
      
      
      full1 <- dtfull %>% dplyr::select(osler_id, interval, pprob1, pprob2, pprob3)
      full1 <- full1[full1$interval == s-1+i,]
      dtc <- dtc %>% left_join(full1 %>% dplyr::select(osler_id, pprob1, pprob2, pprob3), by = 'osler_id')
    }else{
      leng.t <- length(1:s)
      cmeank1.biom1 <- c(rep(NA,leng.t), cmean1[[i]][[1]][inx.biom1] )
      cmeank2.biom1 <- c(rep(NA,leng.t), cmean1[[i]][[2]][inx.biom1] )
      cmeank3.biom1 <- c(rep(NA,leng.t), cmean1[[i]][[3]][inx.biom1] )
      
      cmeank1.biom3 <- c(rep(NA,leng.t), cmean1[[i]][[1]][inx.biom3] )
      cmeank2.biom3 <- c(rep(NA,leng.t), cmean1[[i]][[2]][inx.biom3] )
      cmeank3.biom3 <- c(rep(NA,leng.t), cmean1[[i]][[3]][inx.biom3] )
      
      cmeank1.biom4 <- c(rep(NA,leng.t), cmean1[[i]][[1]][inx.biom4] )
      cmeank2.biom4 <- c(rep(NA,leng.t), cmean1[[i]][[2]][inx.biom4] )
      cmeank3.biom4 <- c(rep(NA,leng.t), cmean1[[i]][[3]][inx.biom4] )
      
      dt.ext <- data.frame(cmeank1.biom1, cmeank2.biom1, cmeank3.biom1,
                           cmeank1.biom3, cmeank2.biom3, cmeank3.biom3,
                           cmeank1.biom4, cmeank2.biom4, cmeank3.biom4)
      dt.ext$interval <- 0:(i+s-1)
      dt.ext$osler_id <- sid
      dtc <- data %>% full_join(dt.ext, by =c('osler_id', "interval"))
      full1 <- dtfull %>% dplyr::select(osler_id, interval, pprob1, pprob2, pprob3)
      full1 <- full1[full1$interval == s-1+i,]
      dtc <- dtc %>% left_join(full1 %>% dplyr::select(osler_id, pprob1, pprob2, pprob3), by = 'osler_id')
    }
    tmp.ave[[i]] <- dtc
    
    dtc.tmp <- dtc%>% filter(interval <= s-1+i) %>% dplyr::rename(obs.biom1 = sao2_norm,
                                                                  obs.biom3 = temp_norm,
                                                                  obs.biom4 = pulse_norm)
    
    dtc.tmp$obs.biom1.inv = predict(normobj.sao2, newdata = dtc.tmp$obs.biom1, inverse = T)
    dtc.tmp$obs.biom3.inv = predict(normobj.temp, newdata = dtc.tmp$obs.biom3, inverse = T)
    dtc.tmp$obs.biom4.inv = predict(normobj.pulse, newdata = dtc.tmp$obs.biom4, inverse = T)
    
    dtc.tmp$cmeank1.biom1.inv = predict(normobj.sao2, newdata = dtc.tmp$cmeank1.biom1, inverse = T)
    dtc.tmp$cmeank2.biom1.inv = predict(normobj.sao2, newdata = dtc.tmp$cmeank2.biom1, inverse = T)
    dtc.tmp$cmeank3.biom1.inv = predict(normobj.sao2, newdata = dtc.tmp$cmeank3.biom1, inverse = T)
    
    dtc.tmp$cmeank1.biom3.inv = predict(normobj.temp, newdata = dtc.tmp$cmeank1.biom3, inverse = T)
    dtc.tmp$cmeank2.biom3.inv = predict(normobj.temp, newdata = dtc.tmp$cmeank2.biom3, inverse = T)
    dtc.tmp$cmeank3.biom3.inv = predict(normobj.temp, newdata = dtc.tmp$cmeank3.biom3, inverse = T)
    
    dtc.tmp$cmeank1.biom4.inv = predict(normobj.pulse, newdata = dtc.tmp$cmeank1.biom4, inverse = T)
    dtc.tmp$cmeank2.biom4.inv = predict(normobj.pulse, newdata = dtc.tmp$cmeank2.biom4, inverse = T)
    dtc.tmp$cmeank3.biom4.inv = predict(normobj.pulse, newdata = dtc.tmp$cmeank3.biom4, inverse = T)
    
    dtc <- dtc.tmp %>% 
      gather(key, cmean, c(obs.biom1.inv:obs.biom4.inv, cmeank1.biom1.inv:cmeank3.biom4.inv))%>% 
      separate(key, c('type', 'biom','inv'), "\\.") %>% 
      spread(type, cmean)
    p <- p + 
      geom_line(data  = dtc, aes(x = interval, y = cmeank1, color = "Discharge"), alpha = 0.6)+
      geom_point(data  = dtc %>% filter(interval == s-1+i), aes(x = interval, y = cmeank1, color = "Discharge", size = pprob1), alpha = 0.6)+
      geom_line(data  = dtc, aes(x = interval, y = cmeank2, color = "Died"), alpha = 0.6)+
      geom_point(data  = dtc %>% filter(interval ==  s-1+i), aes(x = interval, y = cmeank2, color = "Died", size = pprob2), alpha = 0.6)+
      geom_line(data  = dtc, aes(x = interval, y = cmeank3, color = "Vent"), alpha = 0.6)+
      geom_point(data  = dtc %>% filter(interval ==  s-1+i), aes(x = interval, y = cmeank3, color = "Vent", size = pprob3), alpha = 0.6) +
      facet_grid(biom~., scales = "free_y", labeller = labeller(biom = biom.labeller))
    #    p
  }
  dt.ave[[s]] <-tmp.ave
  summ <- data.frame(type = c("Discharge", "Die", "Vent"),
                     sum.P = c(percent(sum(dtfull$pprob1, na.rm = T)),
                               percent(sum(dtfull$pprob2, na.rm = T)),
                               percent(sum(dtfull$pprob3, na.rm = T))))
  colnames(summ) <- c("Event", "Probability")
  
  facet_bounds <- data.frame(ymin = c(90, 36, 60),
                             ymax = c(410, 37.2, 120),
                             breaks = c(50,0.2,10),
                             biom = c("biom1", "biom3", "biom4")) 
  
  pf <- p+
    theme_bw()+
    scale_size_continuous(name=" ", range = c(0,10),breaks = seq(0,1,length.out=5), limits = c(0,1))+ 
    scale_color_manual(
      breaks = c("Discharge", "Died", "Vent"),
      values = c("#2FB3CA","#F1564F","#90C13E"),
      name=" ")+ylab("")+xlab("Day of Hospitalization") + 
    geom_blank(data = facet_bounds, aes(ymin = ymin)) + 
    geom_blank(data = facet_bounds, aes(ymax = ymax))
    #guides(size = F, color = F)

  summ.theme <- ttheme_minimal(core = list(bg_params = list(fill = c("#2FB3CA","#F1564F","#90C13E")),
                                           fg_params=list(cex=2.0)))
  pt <- grid.arrange(tableGrob(summ, theme = summ.theme, rows = NULL ,cols = NULL))
  # plist <- list(pt, pp)
  # p.all[[s]] <- plot_grid(plotlist = plist, align="v", axis = "l", ncol = 1, nrow = 2,
  #                         rel_heights = c(1,3))
  
  
  ## original graph
  my_hist <-pf + guides(color = F) 
  
  ## get the legend 
  tmp <- ggplot_gtable(ggplot_build(my_hist))
  leg <- which(sapply(tmp$grobs, function(x) x$name) ==  "guide-box")
  legend <- tmp$grobs[[leg]]
  
  ## create inset table
  my_table <- pt
  
  plist <- list(my_hist + guides(size = F), arrangeGrob(legend, my_table))
  
  
  p.all[[s]] <- plot_grid(plotlist = plist, align="h", axis = "l", ncol = 2, nrow = 1,
                          rel_widths = c(4,1))
  pred_figure[[s]] <- my_hist
  pred_table[[s]] <- my_table
  
}

## process heatmaps ----- 
## get the legend for heatmap
p.output <- list()
pleft.list <- list()
pright.list <- list()
output.list <- list()
newlegrange <- pr[[1]]
tmp <- ggplot_gtable(ggplot_build(newlegrange))
leg <- which(sapply(tmp$grobs, function(x) x$name) ==  "guide-box")
legend.heat <- tmp$grobs[[leg]]
library(ggpubr)
as_ggplot(legend.heat)
for(i in 1:length(p.all)){
  plist.heat <- list(pr[[i]] + guides(fill = F) + ggtitle(""), 
                     NULL,
                     pl[[i]] + guides(fill = F) + ggtitle(""),
                     NULL,
                     pp[[i]] + guides(fill = F) + ggtitle(""), 
                     legend.heat)
  pleft.list[[i]] <- plot_grid(plotlist = plist.heat, align="hv", axis = "b", ncol = 2, nrow = 3,
                                labels = c("A", "", "B", "", "C"),
                                rel_widths = c(5,1))
  

  pright.list[[i]] <- p.all[[i]]

  output.list[[i]] <- list(pleft.list[[i]], pright.list[[i]])
  p.output[[i]] <- plot_grid(plotlist = output.list[[i]], 
                             align="h", axis = "l", ncol = 2, nrow = 1,
            rel_widths = c(1,1.5),labels = c("", "D"))
}


pdf(paste0("paperPrep/finalizedCode/figs/result_figure_",sid,".pdf"), width = 16, height = 8)
for(i in 1:length(p.output)){
  print(p.output[[i]])
}
dev.off()


## for day 12 specifically
i <- 12
newlegrange <- pp[[i]]+ scale_fill_viridis(limits=c(-6,0))
tmp <- ggplot_gtable(ggplot_build(newlegrange))
leg <- which(sapply(tmp$grobs, function(x) x$name) ==  "guide-box")
legend.heat <- tmp$grobs[[leg]]
for(i in 12){
  plist.heat <- list(pr[[i]] + guides(fill = F) + ggtitle("")+ scale_fill_viridis(limits=c(-6,0)), 
                     NULL,
                     pl[[i]] + guides(fill = F) + ggtitle("")+ scale_fill_viridis(limits=c(-6,0)),
                     NULL,
                     pp[[i]] + guides(fill = F) + ggtitle("")+ scale_fill_viridis(limits=c(-6,0)), 
                     legend.heat)
  pleft.list[[i]] <- plot_grid(plotlist = plist.heat, align="hv", axis = "b", ncol = 2, nrow = 3,
                               labels = c("A", "", "B", "", "C"),
                               rel_widths = c(5,1))
  
  
  pright.list[[i]] <- p.all[[i]]
  
  output.list[[i]] <- list(pleft.list[[i]], pright.list[[i]])
  p.output[[i]] <- plot_grid(plotlist = output.list[[i]], 
                             align="h", axis = "l", ncol = 2, nrow = 1,
                             rel_widths = c(1,1.5),labels = c("", "D"))
}

pdf(paste0("paperPrep/finalizedCode/figs/result_figure_d12_",sid,".pdf"), width = 16, height = 8)
p.output[[i]] 
dev.off()


pdf(paste0("paperPrep/finalizedCode/figs/result_figure_d2&12_",sid,".pdf"), width = 16, height = 18)
grid.arrange(p.output[[2]] , p.output[[12]] , ncol = 1)
dev.off()

