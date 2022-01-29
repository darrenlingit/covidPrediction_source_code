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



## backward ----
dt.pre.mult <- dt1c

dt.pre.mult.train <- dt1c %>% filter((osler_id %in% unique(dt1c$osler_id))) #length(unique(dt.pre.mult.train$osler_id))
dt.pre.mult.train <- dt.pre.mult.train %>% filter(!is.na(endpoint)) %>% filter(endpoint != 0)
#dt.pre.mult.train <- dt.pre.mult.train %>% filter(interval != 0)


dt.pre.mult.train <- dt.pre.mult.train  %>% mutate(sex = factor(sex, levels = c("Female", "Male")),
                                                   smoking = factor(smoking, levels = c("Never Smoker", "Former Smoker", "Current Smoker")),
                                                   charlson_cat = factor(charlson_cat, levels = c("0", "1-2" , "3-4", ">5")),
                                                   demo5 = factor(demo5, levels =   c("Young_White", "Old_White" , "Young_Black", "Old_Black", "Latinx/Other")),
                                                   bmi_cat = factor(bmi_cat, levels = c("<30", ">=30")))
# mult.pre <- multinom(factor(Wt) ~ resp_bl + temp_bl + pulse_bl + sao2fio2_bl +
#                        crp_bl + alc_bl + dd_bl + gfr_bl +
#                        charlson_cat + demo5 + sex + bmi_cat +
#                        smoking + ns(interval, knots = c(2,4, 6 ,8,12,16)),
#                      data = dt.pre.mult.train, maxit = 1000)
# 
# 
# save(mult.pre, file = "paperPrep/finalizedCode/model/alldata_multinom_model.RData")

load("paperPrep/finalizedCode/model/alldata_multinom_model.RData")
t(summary(mult.pre)$coefficients)
t(summary(mult.pre)$standard.errors)


## forward------
dt1c <- dt1c %>% mutate(sao2_fio2_ratio_o = sao2_fio2_ratio_m,
                        temp_c_o = temp_c_m,
                        pulse_o = pulse_m)
dt1c <- dt1c %>%
  # filter(complete.cases(temp_c_m, sao2_fio2_ratio_m, pulse_m)) %>% 
  mutate(sao2_norm = orderNorm(sao2_fio2_ratio_o)$x.t,
         temp_norm = orderNorm(temp_c_o)$x.t,
         pulse_norm = orderNorm(pulse_o)$x.t)

dt.pre.mult <- dt1c

# 1734 unique patients for dt.pre and 1602 for mult.pre
dt.pre.mult.train <- dt1c %>% filter((osler_id %in% unique(dt1c$osler_id))) #length(unique(dt.pre.mult.train$osler_id))
dt.pre.mult.train <- dt.pre.mult.train %>% filter(!is.na(endpoint)) %>% filter(endpoint != 0)
#dt.pre.mult.train <- dt.pre.mult.train %>% filter(interval != 0)



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
dt.pre.mult.train <- dt.pre.mult.train  %>% mutate(sex = factor(sex, levels = c("Female", "Male")),
                                                   smoking = factor(smoking, levels = c("Never Smoker", "Former Smoker", "Current Smoker")),
                                                   charlson_cat = factor(charlson_cat, levels = c("0", "1-2" , "3-4", ">5")),
                                                   demo5 = factor(demo5, levels =   c("Young_White", "Old_White" , "Young_Black", "Old_Black", "Latinx/Other")),
                                                   bmi_cat = factor(bmi_cat, levels = c("<30", ">=30")))

# mult.pre.forw <- multinom(factor(Wt) ~ resp_bl + temp_bl + pulse_bl + sao2fio2_bl +
#                                   crp_bl + alc_bl + dd_bl + gfr_bl +
#                                   charlson_cat + demo5 + sex + bmi_cat +
#                                   smoking + ns(interval, knots = c(4, 6 ,8,12,16)) +
#                                   prev1_sao2 + prev1_temp + prev1_pulse +
#                                   slope_sao2 + slope_temp + slope_pulse,
#                                 na.action = na.omit,
#                                 data = dt.pre.mult.train, maxit = 1000)
# save(mult.pre.forw, file = "paperPrep/finalizedCode/model/alldata_multinom_model_forw.RData")


load("paperPrep/finalizedCode/model/alldata_multinom_model.RData")
coef_bkwr <- t(summary(mult.pre)$coefficients)
se_bkwr <- t(summary(mult.pre)$standard.errors)

load("paperPrep/finalizedCode/model/alldata_multinom_model_forw.RData")
coef_fowr <- t(summary(mult.pre.forw)$coefficients)
se_fowr <- t(summary(mult.pre.forw)$standard.errors)

bkwr_mult_result <- list(coef_bkwr=coef_bkwr, se_bkwr = se_bkwr)
save(bkwr_mult_result, file="paperPrep/finalizedCode/model/bkwr_mult_coef_se.RData")
fowr_mult_result <- list(coef_fowr=coef_fowr, se_fowr = se_fowr)
save(fowr_mult_result, file="paperPrep/finalizedCode/model/fowr_mult_coef_se.RData")

## start here -------- 
load("paperPrep/finalizedCode/model/fowr_mult_coef_se.RData")
load("paperPrep/finalizedCode/model/bkwr_mult_coef_se.RData")
coef_fowr <- fowr_mult_result[[1]]; se_fowr <- fowr_mult_result[[2]]
coef_bkwr <- bkwr_mult_result[[1]]; se_bkwr <- bkwr_mult_result[[2]]

colnames(coef_bkwr) <- colnames(coef_fowr) <-  c("disch", "death", "vent")
colnames(se_bkwr) <- colnames(se_fowr) <- c("se_disch", "se_death", "se_vent")

rows_mult10_bkwr <- c("resp_bl", "temp_bl", "alc_bl")
rows_mult100_bkwr <- c("sao2fio2_bl", "pulse_bl", "dd_bl", "gfr_bl","crp_bl")

rows_mult10_fowr <- c("resp_bl", "temp_bl", "alc_bl", 
                      "prev1_temp" , "slope_temp")
rows_mult100_fowr <- c("sao2fio2_bl", "pulse_bl", "dd_bl", "gfr_bl","crp_bl")

reformat_table <- function(mat, rows_mult10, rows_mult100){
  mat <- as.data.frame(mat)
  mat[rownames(mat) %in% rows_mult10,] <- as.data.frame(mat)[rownames(mat) %in% rows_mult10,]*10
  mat[rownames(mat) %in% rows_mult100,] <- as.data.frame(mat)[rownames(mat) %in% rows_mult100,]*100
  return(mat)
}
coef_bkwr_ref <- reformat_table(coef_bkwr, rows_mult10_bkwr, rows_mult100_bkwr)
se_bkwr_ref <- reformat_table(se_bkwr, rows_mult10_bkwr, rows_mult100_bkwr)

coef_fowr_ref <- reformat_table(coef_fowr, rows_mult10_fowr, rows_mult100_fowr)
se_fowr_ref <- reformat_table(se_fowr, rows_mult10_fowr, rows_mult100_fowr)

bkwr_ref <- cbind(coef_bkwr_ref, se_bkwr_ref)
fowr_ref <- cbind(coef_fowr_ref, se_fowr_ref)

bkwr_ref <- bkwr_ref[!(grepl("ns",rownames(bkwr_ref))),]
cci <- bkwr_ref[grepl("charlson_cat",rownames(bkwr_ref)),]
bkwr_ref <- bkwr_ref[!(grepl("charlson_cat",rownames(bkwr_ref))),]
bkwr_ref <- rbind(bkwr_ref, cci)
bkwr_ref$var <- c("Intercept","Respiratory rate/10, BL", "Temperature/10, BL", 
                  "Pulse/100, BL", "SpO2-FiO2/100, BL",
                  "CRP/100, BL", "ALC/10, BL", "D-dimer/100, BL", "eGFR/100, BL",
                  "White, Age>74", "Black, Age<=74", "Black, Age>74", "Latinx/Other, All ages", 
                  "Male (vs. Female)", "BMI>=30 (vs <30)", "Former smoker", "Current smoker",
                  "CCI 1-2", "CCI 3-4", "CCI >=5")

fowr_ref <- fowr_ref[!(grepl("ns",rownames(fowr_ref))),]
cci <- fowr_ref[grepl("charlson_cat",rownames(fowr_ref)),]
prev_n_slope <- fowr_ref[grepl("prev1|slope",rownames(fowr_ref)),]

fowr_ref <- fowr_ref[!(grepl("charlson_cat|prev1|slope",rownames(fowr_ref))),]
fowr_ref <- rbind(fowr_ref, cci, prev_n_slope)
fowr_ref$var <- c("Intercept","Respiratory rate/10, BL", "Temperature/10, BL", 
                  "Pulse/100, BL", "SpO2-FiO2/100, BL",
                  "CRP/100, BL", "ALC/10, BL", "D-dimer/100, BL", "eGFR/100, BL",
                  "White, Age>74", "Black, Age<=74", "Black, Age>74", "Latinx/Other, All ages", 
                  "Male (vs. Female)", "BMI>=30 (vs <30)", "Former smoker", "Current smoker",
                  "CCI 1-2", "CCI 3-4", "CCI >=5", 
                  "SpO2-FiO2, PV", "Temperature/10, PV", "Pulse, PV",
                  "SpO2-FiO2, PS", "Temperature/10, PS", "Pulse, PS")


bkwr_ref <- bkwr_ref %>% mutate(lwr_disch = disch - qnorm(0.975) * se_disch,
                                upr_disch = disch + qnorm(0.975) * se_disch,
                                lwr_death = death - qnorm(0.975) * se_death,
                                upr_death = death + qnorm(0.975) * se_death,
                                lwr_vent = vent - qnorm(0.975) * se_vent,
                                upr_vent = vent + qnorm(0.975) * se_vent) 
bkwr_ref$var <- factor(bkwr_ref$var, levels = bkwr_ref$var )
bkwr_ref <- bkwr_ref %>% filter(var != "Intercept")
p.bkwr <- ggplot() +
  geom_pointrange(data = bkwr_ref, aes(x = var, y=disch, ymin = lwr_disch, ymax = upr_disch, color = "Discharge"))+
  geom_pointrange(data = bkwr_ref, aes(x = var, y=death, ymin = lwr_death, ymax = upr_death, color = "Die"), position = position_jitter(height = 0.1))+
  geom_pointrange(data = bkwr_ref, aes(x = var, y=vent, ymin = lwr_vent, ymax = upr_vent, color = "Vent"), position = position_jitter(height = 0.1))+
  scale_x_discrete(limit= rev)+
  scale_color_manual(
    breaks = c("Discharge", "Die", "Vent"),
    values = c("#2FB3CA","#F1564F","#90C13E"),
    name=" ")+ylab("")+xlab("Day of Hospitalization") + 
  coord_flip()+
  geom_hline(yintercept = 0, color = "grey")+
  theme_classic()+
  guides(color  = "none")

fowr_ref <- fowr_ref %>% mutate(lwr_disch = disch - qnorm(0.975) * se_disch,
                                upr_disch = disch + qnorm(0.975) * se_disch,
                                lwr_death = death - qnorm(0.975) * se_death,
                                upr_death = death + qnorm(0.975) * se_death,
                                lwr_vent = vent - qnorm(0.975) * se_vent,
                                upr_vent = vent + qnorm(0.975) * se_vent) 
fowr_ref$var <- factor(fowr_ref$var, levels = fowr_ref$var )
fowr_ref <- fowr_ref %>% filter(var != "Intercept")
p.fowr <- ggplot() +
  geom_pointrange(data = fowr_ref, aes(x = var, y=disch, ymin = lwr_disch, ymax = upr_disch, color = "Discharge"))+
  geom_pointrange(data = fowr_ref, aes(x = var, y=death, ymin = lwr_death, ymax = upr_death, color = "Die"), position = position_jitter(height = 0.1))+
  geom_pointrange(data = fowr_ref, aes(x = var, y=vent, ymin = lwr_vent, ymax = upr_vent, color = "Vent"), position = position_jitter(height = 0.1))+
  scale_x_discrete(limit= rev)+
  scale_color_manual(
    breaks = c("Discharge", "Die", "Vent"),
    values = c("#2FB3CA","#F1564F","#90C13E"),
    name=" ")+ylab("")+xlab("Day of Hospitalization") + 
  coord_flip()+
  geom_hline(yintercept = 0, color = "grey")+
  theme_classic()+
  guides(color  = "none")

#grid.arrange(p.fowr, p.bkwr, nrow =1)

fowr_ref$type = "Prospective"
bkwr_ref$type = "Retrospective"
plot_ref <- rbind(fowr_ref, bkwr_ref)
#plot_ref$fac <- factor(unlist(sapply(as.vector(table(plot_ref$var)), seq_len)))


col_disch = plot_ref[,c("disch", "lwr_disch", "upr_disch", "type", "var")]; col_disch$event = "Discharge"
col_death = plot_ref[,c("death", "lwr_death", "upr_death", "type", "var")]; col_death$event = "Die"
col_vent = plot_ref[,c("vent", "lwr_vent", "upr_vent", "type", "var")]; col_vent$event = "Vent"
colnames(col_disch)<-colnames(col_death) <-colnames(col_vent) <- c("mean", "lwr", "upr", "type", "var", "event")

plot_ref2 <- rbind(col_disch, col_death, col_vent)
plot_ref2 <- plot_ref2 %>% mutate(mean_trim = ifelse(mean < -5, -5, 
                                                     ifelse(mean  > 5, 5, mean)),
                                  lwr_trim = ifelse(lwr < -5, -5, 
                                                     ifelse(lwr  > 5, 5, lwr)),
                                  upr_trim = ifelse(upr < -5, -5, 
                                                     ifelse(upr  > 5, 5, upr)))
plot_ref2$event <- factor(plot_ref2$event, levels = c("Vent", "Die", "Discharge"))
fontsize = 26
nm.grey <- unique(plot_ref2$var)[seq(1, length(unique(plot_ref2$var)), 2)]
rectangles <- data.frame(x = nm.grey, y=1)
pmult <- ggplot() +
  # geom_pointrange(data = plot_ref2, aes(x = var,y = mean, ymin = lwr, ymax = upr,  color = event, shape = type),
  #                 position = position_dodge2(width = 1),
  #                 size = 0.8,
  #                 fatten = 1.5)+
  geom_point(data = plot_ref2, aes(x = var,y = mean,   color = event, shape = type),
             position = position_dodge2(width = 1), size = 6)+
  geom_linerange(data = plot_ref2, aes(x = var,y = mean, ymin = lwr, ymax = upr,  color = event),
                 position = position_dodge2(width = 1), size = 1)+
  coord_flip(ylim  = c(-8,8))+
  scale_x_discrete(limit= rev, name="")+
  scale_color_manual(
    labels  = c("Discharge", "Died", "Vent"),
    limits=c("Discharge", "Die", "Vent"),
    values = c("Discharge"="#2FB3CA",
               "Die" = "#F1564F",
               "Vent" = "#90C13E"),
    name="")+
  ylab("Coefficients of Multinomial Regression")+
  geom_hline(yintercept = 0, color = "grey")+
  theme_minimal()+
  theme(axis.title = element_text(size = fontsize),
        axis.text = element_text(size = fontsize),
        strip.text = element_text(size = fontsize),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = c(0.15,0.93),
        legend.background = element_rect(fill = NA,colour = NA),
        legend.text = element_text(size = fontsize-2),
        legend.spacing = unit(-1, "cm"),
        legend.key.size=unit(1.5,"cm"))+
  scale_shape_manual(name = "",
                       breaks = c("Retrospective", "Prospective"),
                       values = c("\u25C4","\u25BA"))+
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))

cairo_pdf("paperPrep/finalizedCode/figs/multinom_results_redim.pdf", width = 12, height = 26)
pmult
dev.off()


###
plot_ref2$var <- factor(plot_ref2$var, levels = unique(plot_ref2$var))
g <- ggplot(data = plot_ref2, aes(x = mean,y = var))+
  geom_vline(xintercept = 0, linetype = "solid", colour = "black", size = 0.5, alpha = 0.6)+
  
  geom_effect(
    ggplot2::aes(xmin = lwr, xmax=upr,
                 colour = event, shape = type),
    fatten = 4, size = 1,
    position = ggstance::position_dodgev(height =1)
  )
g2 <- g +
  theme_minimal()+
  geom_stripes(odd = '#11111111', even = '#00000000')+
  scale_color_manual(
    labels = c("Discharge", "Died", "Vent"),
    breaks = c("Discharge", "Die", "Vent"),
    values = c("#2FB3CA","#F1564F","#90C13E"),
    name=" ", 
    limits = (levels(plot_ref2$event)))+xlab("Coefficients of Multinomial Regression")+ylab("") + 
  theme(axis.title = element_text(size = fontsize),
        axis.text = element_text(size = fontsize),
        strip.text = element_text(size = fontsize),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = c(0.2,0.95),
        legend.background = element_rect(fill = NA,colour = NA),
        legend.spacing = unit(-0.8, "cm"),
        legend.key =  element_blank(),
        legend.text = element_text(size = fontsize - 1))+
       # legend.margin = margin(0,0,0,0,"cm"))+
  scale_shape_manual(name = "",
                     breaks = c("Retrospective", "Prospective"),
                     values = c(24,25))+
  scale_y_discrete(limits = rev(levels(plot_ref2$var)))+
  guides(color = guide_legend(order = 1,override.aes=list(size = 5), keywidth = 1, keyheight = 1, default.unit = "cm"), 
         shape = guide_legend(order = 2, override.aes=list(size = 5), keywidth = 1, keyheight = 1, default.unit = "cm"))
g2

pdf("paperPrep/finalizedCode/figs/multinom_results_stripe.pdf", width = 12, height = 26)
g2
dev.off()

