## backward prediction AUC
## given data on day 0, day 2, 4, and 8 
rm(list=ls())
nbiom <- 3
df.re <- 3
nb = 2000; nt = 8000
source("00_loadlibs.R")
setwd(wd)
#load data
dt.pre
dt.pre.mult.tmp

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


i <- 2 #given data until day 2
## filter out patients have an event prior to day 4, keep only event occurred after day 4.
dt.endpoints <- dt.pre.mult.tmp %>% 
  dplyr::group_by(osler_id) %>% 
  dplyr::summarise(endpoint.u = min(endpoint), status.u = min(status))%>% 
  dplyr::select(osler_id, endpoint.u, status.u) %>% 
  unique() %>% 
  filter(endpoint.u > i)

load(paste0(path, "bkwr_di=", i ,"_M.RData"))
load(paste0(path, "forw_di=", i ,"_M.RData"))

f.prob.save <- forw.save2[[5]]
b.prob.save <- bkwr.save2[[5]]

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

#test <- data.roc %>% dplyr::select(osler_id, cv) %>% unique();table(test$cv)

b_auc_uet<- f_auc_uet <- NULL
for(t in (i+1):20){
  # backward
  b.tmp.auc <- NULL
  datasurv <- b.data.roc %>% filter(day == t)
  auc1 <- auc(datasurv$cumO1,datasurv$cumE1)
  auc2 <- auc(datasurv$cumO2,datasurv$cumE2)
  auc3 <- auc(datasurv$cumO3,datasurv$cumE3)

  b.tmp.auc <- data.frame(event = c("disch", "death", "vent"), 
                          day = rep(t, 3), 
                          mean = c(auc1, auc2, auc3))

  b_auc_uet <- rbind(b_auc_uet, b.tmp.auc)
  
  # forward
  f.tmp.auc <- NULL
  datasurv <- f.data.roc %>% filter(day == t)
  auc1 <- auc(datasurv$cumO1,datasurv$cumE1)
  auc2 <- auc(datasurv$cumO2,datasurv$cumE2)
  auc3 <- auc(datasurv$cumO3,datasurv$cumE3)

  f.tmp.auc <- data.frame(event = c("disch", "death", "vent"), 
                          day = rep(t, 3), 
                          mean = c(auc1, auc2, auc3))
  f_auc_uet <- rbind(f_auc_uet, f.tmp.auc)
}

m.diff <- sum(b_auc_uet$mean - f_auc_uet$mean)/((max(f_auc_uet$day)-min(f_auc_uet$day)))
lwr.diff <- sum(b_auc_uet$lwr - f_auc_uet$lwr)/((max(f_auc_uet$day)-min(f_auc_uet$day)))
upr.diff <- sum(b_auc_uet$upr - f_auc_uet$upr)/((max(f_auc_uet$day)-min(f_auc_uet$day)))
diff <- paste0(round(m.diff,3), "(", round(lwr.diff,3), ",", round(upr.diff,3), ")")

pauc<-ggplot()+
  geom_line(data = b_auc_uet,aes(x= day, y = mean, group = event, color = event, linetype = "backward"))+
  geom_line(data = f_auc_uet,aes(x= day, y = mean, group = event, color = event, linetype = "forward"))+
  theme_classic()+
  ylim(0.5,1)+
  ylab("AUC")+
  xlab("Days projected")+
  scale_color_manual(breaks = c("disch", "death", "vent"),
                     label = c("Discharge", "Die", "Vent"),
                     values = c("#2FB3CA","#F1564F","#90C13E"),name=" ")+
  scale_fill_manual(breaks = c("disch", "death", "vent"),
                    values = c("#2FB3CA","#F1564F","#90C13E"),
                    label = c("Discharge", "Die", "Vent"), name=" ")+
  scale_linetype(name = "")+
  ggtitle(paste0("AUC(u=",i,", e, t)"))+
  xlim(0,20)

auc2 <- list(pauc, diff, b_auc_uet, f_auc_uet)
save(auc2,
     file = paste0("pmodel/auc_di=", i,"_M.RData"))







