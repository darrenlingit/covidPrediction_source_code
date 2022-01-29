plot(sqrt(data.all$exp.tot[data.all$type == "Retro"]),
     sqrt(data.all$exp.tot[data.all$type == "Prosp"]))
abline(0,1)
#average prospecitve and retrospective
ave_pro_retro <- data.all %>% 
  ungroup() %>% 
  group_by(day, event, U) %>% 
  summarize(ave_preds = mean(exp.tot)) %>% 
  arrange(U, day, event)
data.all <- data.all %>% arrange(U, day, event)

obs <- sqrt(data.all$obs.tot[data.all$type == "Retro"])
pros <- sqrt(data.all$exp.tot[data.all$type == "Prosp"])
retrosp <- sqrt(data.all$exp.tot[data.all$type == "Retro"])
average <-sqrt(ave_pro_retro$ave_preds)
dt_pairs <- cbind(obs, pros, retrosp, average)
pairs(dt_pairs, upper.panel = my_line)
my_line <- function(x,y,...){
  points(x,y,... )
  abline(a= 0,b=1,...)
}

cor(dt_pairs)
lm(obs ~ pros + retrosp,data = as.data.frame(dt_pairs))

cor(data.all$exp.tot[data.all$type == "Retro"],
     data.all$exp.tot[data.all$type == "Prosp"])

unique.event <- c("Discharge", "Die", "Vent")
par(mfrow = c(2,2))
for(i in 1:length(unique.event)){
  Retrosp <- data.all$exp.tot[data.all$type == "Retro" & data.all$event == unique.event[i]]
  Prosp <- data.all$exp.tot[data.all$type == "Prosp" & data.all$event == unique.event[i]]
  plot(Retrosp,
       Prosp)
  cori <- round(cor(Retrosp, Prosp, method = "spearman" ),2)
  text(max(Retrosp)/4, 3*max(Prosp)/4, paste(unique.event[i], cori))
}

unique.U <- c(0,2,4,8)
par(mfrow = c(2,2))
for(i in 1:length(unique.U)){
  Retrosp <- data.all$exp.tot[data.all$type == "Retro" & data.all$U == unique.U[i]]
  Prosp <- data.all$exp.tot[data.all$type == "Prosp" & data.all$U == unique.U[i]]
  plot(Retrosp,
       Prosp)
  cori <- round(cor(Retrosp, Prosp, method = "spearman" ),2)
  text(max(Retrosp)/4, 3*max(Prosp)/4, paste(unique.U[i], cori))
}
