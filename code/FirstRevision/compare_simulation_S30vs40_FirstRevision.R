
path<-"paperPrep/finalizedCode/figs/calibration/"
load(paste0(path,"id_remove_2.RData")); id_remove2 <- id.remove
load(paste0(path,"id_remove_4.RData")); id_remove4 <- id.remove
load(paste0(path,"id_remove_8.RData")); id_remove8 <- id.remove
id_remove2  <- as.data.frame(id_remove2);id_remove4  <- as.data.frame(id_remove4);id_remove8  <- as.data.frame(id_remove8)
colnames(id_remove2) <- colnames(id_remove4) <- colnames(id_remove8) <- c("cv", "position", "osler_id")


i <- 0 #given data until day 2
## filter out patients have an event prior to day 4, keep only event occurred after day 4.
dt.endpoints <- dt.pre.mult.tmp %>% 
  dplyr::group_by(osler_id) %>% 
  dplyr::summarise(endpoint.u = min(endpoint), status.u = min(status))%>% 
  dplyr::select(osler_id, endpoint.u, status.u) %>% 
  unique() %>% 
  filter(endpoint.u > i)

load(paste0(path, "forw_di=", i ,"_S30.RData"))

load(paste0(path, "forw_di=", i ,"_M.RData"))
load(paste0(path, "forw_di=",  "0.RData"))

b.prob.save <- forw.save8[[5]]
f.prob.save <- forw.save8.S30[[5]]
id_remove <- id_remove8

b.prob.save <- forw.save0[[5]]
f.prob.save <- forw.save0.S30[[5]]
id_remove <- NULL


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
      # if(i == 0){data <- cbind(data, day = 1:20)}else{data <- cbind(data, day = (i+1):20)}
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

#test <- data.roc %>% dplyr::select(osler_id, cv) %>% unique();table(test$cv)

library(pROC)
b_auc_uet<- f_auc_uet <- NULL
for(t in (i+1):20){
  # backward
  b.tmp.auc <- NULL
  datasurv <- b.data.roc %>% filter(day == t)
  auc1 <- auc(datasurv$cumO1,datasurv$cumE1)
  #  lwr1 <- ci.auc(auc1, method = "bootstrap")[1]; upr1 <- ci.auc(auc1, method = "bootstrap")[3]
  auc2 <- auc(datasurv$cumO2,datasurv$cumE2)
  #  lwr2 <- ci.auc(auc2, method = "bootstrap")[1]; upr2 <- ci.auc(auc2, method = "bootstrap")[3]
  auc3 <- auc(datasurv$cumO3,datasurv$cumE3)
  #  lwr3 <- ci.auc(auc3, method = "bootstrap")[1]; upr3 <- ci.auc(auc3, method = "bootstrap")[3]
  
  b.tmp.auc <- data.frame(event = c("disch", "death", "vent"), 
                          day = rep(t, 3), 
                          mean = c(auc1, auc2, auc3))
  #lwr = c(lwr1, lwr2, lwr3),
  #upr = c(upr1, upr2, upr3))
  b_auc_uet <- rbind(b_auc_uet, b.tmp.auc)
  
  # forward
  f.tmp.auc <- NULL
  datasurv <- f.data.roc %>% filter(day == t)
  auc1 <- auc(datasurv$cumO1,datasurv$cumE1)
  # lwr1 <- ci.auc(auc1, method = "bootstrap")[1]; upr1 <- ci.auc(auc1, method = "bootstrap")[3]
  auc2 <- auc(datasurv$cumO2,datasurv$cumE2)
  # lwr2 <- ci.auc(auc2, method = "bootstrap")[1]; upr2 <- ci.auc(auc2, method = "bootstrap")[3]
  auc3 <- auc(datasurv$cumO3,datasurv$cumE3)
  # lwr3 <- ci.auc(auc3, method = "bootstrap")[1]; upr3 <- ci.auc(auc3, method = "bootstrap")[3]
  
  f.tmp.auc <- data.frame(event = c("disch", "death", "vent"), 
                          day = rep(t, 3), 
                          mean = c(auc1, auc2, auc3))
  #lwr = c(lwr1, lwr2, lwr3),
  #upr = c(upr1, upr2, upr3))
  f_auc_uet <- rbind(f_auc_uet, f.tmp.auc)
}

pauc<-ggplot()+
  geom_line(data = b_auc_uet,aes(x= day, y = mean, group = event, color = event, linetype = "S=40"))+
  geom_line(data = f_auc_uet,aes(x= day, y = mean, group = event, color = event, linetype = "S=30"))+
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

dtf0 <- forw.save0[[1]]$data; dtf0$type = "Prosp"; dtf0$U = 0; dtf0$S <- 40
dtf0S30 <- forw.save0.S30[[1]]$data; dtf0S30$type = "Prosp"; dtf0S30$U = 0; dtf0S30$S <- 30
data.all <- rbind(dtf0, dtf0S30)
pcal <- ggplot(data.all) +
  geom_point(data = data.all, aes(x = day, y = obs.tot, color = event, shape = "Observed")) +
  geom_line(data = data.all, aes(x = day, y = exp.tot, color = event, linetype = factor(S))) +
  # geom_vline(data = data.U, aes(xintercept = nstay+1), color = "grey")+
  # facet_rep_grid(U~.)+
  theme_classic()+
  scale_color_manual(label = c("Discharge", "Died", "Vent"),
                     breaks = c("Discharge", "Die", "Vent"),
                     values = c("#2FB3CA","#F1564F","#90C13E"),
                     name="")+
  xlab("Day of Hospitalization")+
  ylab("Number of Events")+
  theme(strip.text = element_blank())+
  scale_y_continuous(trans = "sqrt", position = "right",breaks = c(0,5,20,50,100,150,200))+
  scale_x_continuous(breaks = c(0,5,10,15,20), limits = c(0,20))

pcal
pdf("paperPrep/finalizedCode/figs/compare_simulationS30vs40_d0.pdf", 
    width = 10, height = 6)
grid.arrange(pcal + guides(linetype = F, shape = F) +ggtitle(""),
             pauc + ggtitle(""), nrow = 1)
dev.off()
# auc2 <- list(pauc, diff, b_auc_uet, f_auc_uet)
# save(auc2,
#      file = paste0("paperPrep/finalizedCode/figs/calibration/auc_di=", i,"_M.RData"))







