rm(list=ls())
path <- "model/"

load(paste0(path, "bkwr_di=0.RData"))
load(paste0(path, "forw_di=0.RData"))
load(paste0(path, "bkwr_di=2_M.RData"))
load(paste0(path, "forw_di=2_M.RData"))
load(paste0(path, "bkwr_di=4_M.RData"))
load(paste0(path, "forw_di=4_M.RData"))
load(paste0(path, "bkwr_di=8_M.RData"))
load(paste0(path, "forw_di=8_M.RData"))
load(paste0(path, "auc_di=0.RData"))
load(paste0(path, "auc_di=2_M.RData"))
load(paste0(path, "auc_di=4_M.RData"))
load(paste0(path, "auc_di=8_M.RData"))

dtb0 <- bkwr.save0[[1]]$data; dtb0$type = "Retro"; dtb0$U = 0
dtb2 <- bkwr.save2[[1]]$data; dtb2$type = "Retro"; dtb2$U = 2
dtb4 <- bkwr.save4[[1]]$data; dtb4$type = "Retro"; dtb4$U = 4
dtb8 <- bkwr.save8[[1]]$data; dtb8$type = "Retro"; dtb8$U = 8

dtf0 <- forw.save0[[1]]$data; dtf0$type = "Prosp"; dtf0$U = 0
dtf2 <- forw.save2[[1]]$data; dtf2$type = "Prosp"; dtf2$U = 2
dtf4 <- forw.save4[[1]]$data; dtf4$type = "Prosp"; dtf4$U = 4
dtf8 <- forw.save8[[1]]$data; dtf8$type = "Prosp"; dtf8$U = 8
data.all <- rbind(dtb0, dtb2, dtb4, dtb8, dtf0, dtf2, dtf4, dtf8)
data.U <- data.all %>% group_by(U) %>% summarize(nstay = U) %>% unique()
library(lemon)
pj<-ggplot(data.all) +
  geom_point(data = data.all, aes(x = day, y = obs.tot, color = event), shape=18) +
  geom_line(data = data.all, aes(x = day, y = exp.tot, color = event, linetype = type)) +
  geom_vline(data = data.U, aes(xintercept = nstay+1), color = "grey")+
  facet_rep_grid(U~.)+
  theme_classic()+
  scale_color_manual(breaks = c("Discharge", "Die", "Vent"),
                     values = c("#2FB3CA","#F1564F","#90C13E"),
                     name="")+
  xlab("Day of Hospitalization")+
  ylab("Number of events")+
  theme(strip.text = element_blank())+
  scale_y_continuous(trans = "sqrt")+
  scale_x_continuous(breaks = c(0,5,10,15,20), limits = c(0,20))+
  scale_linetype_manual(name = "",
                        breaks = c("Prosp", "Retro"),
                        labels = c("Prospective","Retrospective"),
                        values = c(2,1))

auc.b0 <- auc0[[3]]; auc.b0$type = "Retro"; auc.b0$U = 0
auc.b2 <- auc2[[3]]; auc.b2$type = "Retro"; auc.b2$U = 2
auc.b4 <- auc4[[3]]; auc.b4$type = "Retro"; auc.b4$U = 4
auc.b8 <- auc8[[3]]; auc.b8$type = "Retro"; auc.b8$U = 8

auc.f0 <- auc0[[4]]; auc.f0$type = "Prosp"; auc.f0$U = 0
auc.f2 <- auc2[[4]]; auc.f2$type = "Prosp"; auc.f2$U = 2
auc.f4 <- auc4[[4]]; auc.f4$type = "Prosp"; auc.f4$U = 4
auc.f8 <- auc8[[4]]; auc.f8$type = "Prosp"; auc.f8$U = 8

auc.all <- rbind(auc.b0,auc.b2,auc.b4,auc.b8,
                 auc.f0, auc.f2, auc.f4, auc.f8)
fkdt <- data.frame(day = 100,mean = 0)
pa <- ggplot(auc.all)+
  geom_point(data=fkdt, aes(x= day, y = mean,shape = "Observed"))+
  geom_line(data = auc.all, aes(x= day, y = mean, group = interaction(type, event), color = event, linetype = type))+
  facet_rep_grid(U~.)+
  theme_classic()+
  ylab("AUC")+
  xlab("Days projected")+
  scale_color_manual(breaks = c("disch", "death", "vent"),
                     label = c("Discharge", "Die", "Vent"),
                     values = c("#2FB3CA","#F1564F","#90C13E"),name=" ")+
  theme(strip.background = element_blank(),
        strip.text = element_blank())+
  scale_linetype_manual(name = "",
                        breaks = c("Prosp", "Retro"),
                        values = c(2,1))+
  scale_y_continuous(breaks = c(0.5,0.6,0.7,0.8,0.9), limits = c(0.5,1))+
  scale_x_continuous(breaks = c(0,5,10,15,20), limits = c(0,20))+
  scale_shape_manual(name = "",
                     breaks = c("Observed"), values = 18)

grid.arrange(pj + guides(linetype = F) + 
               theme(legend.position = c(0.8,1),
                     legend.background = element_rect(fill = NA,colour = NA)),
             pa + guides(color = F)+ 
               theme(legend.position = c(0.8,1),
                     legend.background = element_rect(fill = NA,colour = NA)),ncol = 2)

fontsize <- 18

grid.arrange(pj + guides(linetype = F, shape = F) + 
               theme(legend.position = c(0.8,0.95),
                     legend.background = element_rect(fill = NA,colour = NA),
                     legend.text = element_text(size = fontsize),
                     legend.key.size = unit(0.8,'cm'),
                     axis.title = element_text(size = fontsize),
                     axis.text = element_text(size = fontsize)),
             pa + guides(color = F)+ 
               theme(legend.position = c(0.8,0.85),
                     legend.background = element_rect(fill = NA,colour = NA),
                     legend.text = element_text(size = fontsize),
                     legend.key.size = unit(0.8,'cm'),
                     axis.title = element_text(size = fontsize),
                     axis.text = element_text(size = fontsize),
                     legend.spacing = unit(-1, "cm")),ncol = 2)


