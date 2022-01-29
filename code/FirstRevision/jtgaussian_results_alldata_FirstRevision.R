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
  dplyr::select(osler_id, interval, rev_day, status, sao2_fio2_ratio_m,temp_c_m, pulse_m, outcome)%>% 
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
                pulse_o = pulse_m) 
library(bestNormalize)
dt.pre <- dt.pre %>% 
  filter(complete.cases(sao2_fio2_ratio_o, temp_c_o, pulse_o)) %>% 
  mutate(sao2_norm = orderNorm(sao2_fio2_ratio_o)$x.t,
         temp_norm = orderNorm(temp_c_o)$x.t,
         pulse_norm = orderNorm(pulse_o)$x.t)

dt.pre <- dt.pre  %>% mutate(sex = factor(sex, levels = c("Female", "Male")),
                             smoking = factor(smoking, levels = c("Never Smoker", "Former Smoker", "Current Smoker")),
                             charlson_cat = factor(charlson_cat, levels = c("0", "1-2" , "3-4", ">5")),
                             demo5 = factor(demo5, levels =   c("Young_White", "Old_White" , "Young_Black", "Old_Black", "Latinx/Other")),
                             bmi_cat = factor(bmi_cat, levels = c("<30", ">=30")))
# fit.pre = MCMCglmm(cbind(sao2_norm, temp_norm, pulse_norm) ~
#                      -1 + trait + trait:(status*ns(rev_day,knots = c(-7,-4,-2,-1), Boundary.knots = c(-20,0)) +
#                                            resp_bl + temp_bl + pulse_bl + sao2fio2_bl +
#                                            crp_bl + alc_bl + dd_bl + gfr_bl +
#                                            charlson_cat + demo5 + sex + bmi_cat +
#                                            smoking),
#                    random = ~ us(trait + trait:ns(rev_day,knots = c(-3), Boundary.knots = c(-20,0))):osler_id,  #b0 + b1day for each outcome/trait
#                    rcov = ~ us(trait):units, #unstructured, coursenote pg67
#                    burnin = nb, nitt = nt, pr = T,pl=T,
#                    family = rep("gaussian", nbiom),
#                    data = dt.pre %>% filter(interval!=0))
# save(fit.pre, file = "paperPrep/finalizedCode/model/alldata_jtgaussian_model.RData")

# set.seed(2021)
# fit.pre.forw = MCMCglmm(cbind(sao2_norm, temp_norm, pulse_norm) ~
#                      -1 + trait + trait:(ns(interval,Boundary.knots = c(0,20), knots = c(1,3,6)) +
#                                            resp_bl + temp_bl + pulse_bl + sao2fio2_bl +
#                                            crp_bl + alc_bl + dd_bl + gfr_bl +
#                                            charlson_cat + demo5 + sex + bmi_cat +
#                                            smoking),
#                    random = ~ us(trait + trait:interval):osler_id,  #b0 + b1day for each outcome/trait
#                    rcov = ~ us(trait):units, #unstructured, coursenote pg67
#                    burnin = nb, nitt = nt, pr = T,
#                    family = rep("gaussian", 3),
#                    data = dt.pre %>% filter(interval != 0))
# save(fit.pre.forw, file = "paperPrep/finalizedCode/model/alldata_jtgaussian_model_forw.RData")

## backward in time-----
load("paperPrep/finalizedCode/model/alldata_jtgaussian_model.RData")
normobj.sao2 <- orderNorm(dt.pre$sao2_fio2_ratio_o)
normobj.temp <- orderNorm(dt.pre$temp_c_o)
normobj.pulse <- orderNorm(dt.pre$pulse_o)

# normobj <- list(normobj.sao2, normobj.temp, normobj.pulse)
# save(normobj, file = "paperPrep/finalizedCode/model/normobj.RData")

fit.pre <- fit.pre
Bhat <- fixef(fit.pre, use = "mean")
rn.Bhat <- rownames(Bhat)
order.Bhat <- c(grep("sao2_norm", rn.Bhat),
                grep("temp_norm", rn.Bhat),
                grep("pulse_norm", rn.Bhat))
Bhat <- matrix(Bhat[order.Bhat], ncol = 1)
rownames(Bhat) <- rn.Bhat[order.Bhat]

smy <- summary(fit.pre)
smy.sol <- smy$solutions
nms <- rownames(smy.sol)
order.sol <- c(grep("sao2_norm", nms),
               grep("temp_norm", nms),
               grep("pulse_norm", nms))
sol <- (smy.sol[order.sol,])

grp <- 3
dt.baseline <- dt.pre %>% dplyr::select(resp_bl, temp_bl, pulse_bl, sao2fio2_bl, crp_bl, alc_bl,
                                        dd_bl, gfr_bl, charlson_cat, demo5, sex, bmi_cat, smoking) %>% unique()

newdata <- data.frame(osler_id = "Patient",
                      interval = rep(0:20, grp),
                      rev_day = rep(-20:0, grp),
                      status = rep(1:3, each  = 21),
                      resp_bl = mean(dt.baseline$resp_bl),
                      temp_bl = mean(dt.baseline$temp_bl),
                      pulse_bl = mean(dt.baseline$pulse_bl),
                      crp_bl = mean(dt.baseline$crp_bl),
                      alc_bl = mean(dt.baseline$alc_bl),
                      sao2fio2_bl = mean(dt.baseline$sao2fio2_bl),
                      dd_bl = mean(dt.baseline$dd_bl),
                      gfr_bl = 2,
                      charlson_cat = "0",
                      demo5 = "Young_White",
                      sex = "Female",
                      bmi_cat = "<30",
                      smoking ="Never Smoker",
                      sao2_fio2_ratio_m = NA,
                      temp_c_m = NA,
                      pulse_m = NA)
newdata <- newdata %>%  mutate(sex = factor(sex, levels = c("Female", "Male")),
                               smoking = factor(smoking, levels = c("Never Smoker", "Former Smoker", "Current Smoker")),
                               charlson_cat = factor(charlson_cat, levels = c("0", "1-2" , "3-4", ">5")),
                               demo5 = factor(demo5, levels =   c("Young_White", "Old_White" , "Young_Black", "Old_Black", "Latinx/Other")),
                               bmi_cat = factor(bmi_cat, levels = c("<30", ">=30")),
                               status = factor(status, levels = c("1", "2", "3")))


fm.Xi <- as.formula(paste0("~status*ns(rev_day,knots = c(-7,-4,-2,-1), Boundary.knots = c(-20,0)) +
                                           resp_bl + temp_bl + pulse_bl + sao2fio2_bl +
                                           crp_bl + alc_bl + dd_bl + gfr_bl +
                                           charlson_cat + demo5 + sex + bmi_cat +
                                           smoking"))
Xi1 <- model.matrix(fm.Xi, data = newdata) #;colnames(Xi1) <- rownames(Bhat)[grep("sao2_fio2_ratio_m", rownames(Bhat))]
Xi2 <- model.matrix(fm.Xi, data = newdata)
Xi3 <- model.matrix(fm.Xi, data = newdata)
Xi <- as.matrix(bdiag(list(Xi1, Xi2, Xi3)))#Xi <- as.matrix(bdiag(list(Xi1, Xi3)))#

Yhati <-  Xi %*% Bhat
newdata$sao_fitted = Yhati[1:nrow(Xi1)]
newdata$temp_fitted = Yhati[(nrow(Xi1)+1):(2*nrow(Xi1))]
newdata$pulse_fitted = Yhati[(2*nrow(Xi1)+1):(3*nrow(Xi1))]

v.samp <- fit.pre$Sol[,1:length(Bhat)] # 600 by 111
nms <- colnames(v.samp)
order.v <- c(grep("sao2_norm", nms),
             grep("temp_norm", nms),
             grep("pulse_norm", nms))
v.samp <- v.samp[order.v,order.v]
varFix <- var(v.samp)


predvar <- diag(Xi %*% varFix %*% t(Xi))
predse <- sqrt(predvar)
newdata$sao_se = predse[1:nrow(Xi1)]
newdata$temp_se = predse[(nrow(Xi1)+1):(2*nrow(Xi1))]
newdata$pulse_se = predse[(2*nrow(Xi1)+1):(3*nrow(Xi1))]


## get observations backwards in time
dt.joint <- readRDS("data/data_for_joint.rds")
tmp <- dt.joint %>% group_by(osler_id) %>% dplyr::summarize(maxstatus = max(status))
dt <- dt.joint %>% left_join(tmp, by = 'osler_id') %>% filter(complete.cases(status))
load("S:/COVID_LDA/code/zwang_mgbowring/data/vitals_data0809.RData")
dt.bls <- vitals_merge %>% dplyr::select(osler_id, status,first_vented, outcome) %>% unique()
dt <- dt %>% left_join(dt.bls)

dtp <- dt %>% filter(status %in% c(1,2,3))
tmp <- dtp %>% dplyr::group_by(osler_id) %>% dplyr::summarize(endpoint = max(rev_day)) %>% mutate(rev_day = endpoint)
dtp <- dtp %>% left_join(tmp)

dtp <- dtp %>% dplyr::filter(rev_day >= -15) %>% 
  mutate(max_interval_trim = ifelse(endpoint < -15, NA, endpoint))

## sample 50 persons
idlist <- unique(dtp$osler_id[!is.na(dtp$max_interval_trim)])
set.seed(2022)
ids <- dtp %>% 
  dplyr::select(osler_id, status) %>% 
  unique() %>% 
  group_by(status) %>% sample_n(100)
dtps <- dtp %>% filter(osler_id %in% ids$osler_id)

dtps <- dtps %>% filter(complete.cases(sao2_fio2_ratio_m, pulse_m, temp_c_m))
dtps$sao2_norm <- predict(normobj.sao2, dtps$sao2_fio2_ratio_m)
dtps$temp_norm <- predict(normobj.temp, dtps$temp_c_m)
dtps$pulse_norm <- predict(normobj.pulse, dtps$pulse_m)

facet.lab <- c("Discharge", "Die", "Ventilation")
names(facet.lab) <- c("1", "2", "3")

newdata.p <- newdata %>% filter(rev_day >= -15)

newdata.p$sao2_inv <- predict(normobj.sao2, newdata = newdata.p$sao_fitted, inverse = T)
newdata.p$temp_inv <- predict(normobj.temp, newdata = newdata.p$temp_fitted, inverse = T)
newdata.p$pulse_inv <- predict(normobj.pulse, newdata = newdata.p$pulse_fitted, inverse = T)

newdata.p <- newdata.p %>% mutate(sao2_lwr = sao_fitted - 1.96*sao_se,
                                  sao2_upr = sao_fitted + 1.96*sao_se,
                                  temp_lwr = temp_fitted - 1.96*temp_se,
                                  temp_upr = temp_fitted + 1.96*temp_se,
                                  pulse_lwr = pulse_fitted - 1.96*pulse_se,
                                  pulse_upr = pulse_fitted + 1.96*pulse_se)
newdata.p$sao2_upr_inv <- predict(normobj.sao2, newdata = newdata.p$sao2_upr, inverse = T)
newdata.p$sao2_lwr_inv <- predict(normobj.sao2, newdata = newdata.p$sao2_lwr, inverse = T)
newdata.p$pulse_upr_inv <- predict(normobj.pulse, newdata = newdata.p$pulse_upr, inverse = T)
newdata.p$pulse_lwr_inv <- predict(normobj.pulse, newdata = newdata.p$pulse_lwr, inverse = T)
newdata.p$temp_upr_inv <- predict(normobj.temp, newdata = newdata.p$temp_upr, inverse = T)
newdata.p$temp_lwr_inv <- predict(normobj.temp, newdata = newdata.p$temp_lwr, inverse = T)

dtps_sao2 <- dtps[,c("osler_id", "rev_day", "endpoint", "status", "sao2_fio2_ratio_m")]; dtps_sao2$biomarker = "sao2"
dtps_temp <- dtps[,c("osler_id", "rev_day", "endpoint", "status", "temp_c_m")];dtps_temp$biomarker = "temp"
dtps_pulse <- dtps[,c("osler_id", "rev_day", "endpoint", "status", "pulse_m")];dtps_pulse$biomarker = "pulse"
colnames(dtps_sao2) <- colnames(dtps_temp) <- colnames(dtps_pulse) <-
  c("osler_id", "rev_day", "endpoint", "status", "observed", "biomarker")
dtps_ref <- rbind(dtps_sao2, dtps_temp, dtps_pulse)
dtps_ref$biomarker = factor(dtps_ref$biomarker, levels = c("sao2", "temp", "pulse"))

newdata_sao2 <- newdata.p[,c("osler_id", "rev_day", "status","sao2_inv", "sao2_lwr_inv", "sao2_upr_inv")];
newdata_sao2$biomarker = "sao2"
newdata_temp <- newdata.p[,c("osler_id", "rev_day", "status","temp_inv", "temp_lwr_inv", "temp_upr_inv")];
newdata_temp$biomarker = "temp"
newdata_pulse <- newdata.p[,c("osler_id", "rev_day", "status","pulse_inv", "pulse_lwr_inv", "pulse_upr_inv")];
newdata_pulse$biomarker = "pulse"
colnames(newdata_sao2) <- colnames(newdata_temp) <- colnames(newdata_pulse) <- 
  c("osler_id", "rev_day", "status", "mean","lwr","upr", "biomarker")
newdata.p_ref <- rbind(newdata_sao2, newdata_temp, newdata_pulse)
newdata.p_ref$biomarker = factor(newdata.p_ref$biomarker, levels = c("sao2", "temp", "pulse"))

fontsize = 16
p_retro <- ggplot()+
  geom_line(data = dtps_ref, 
            aes(x = rev_day, y = observed, group = osler_id), color = "grey")+
  geom_point(data = dtps_ref, 
             aes(x = endpoint, y = observed, group = osler_id, color = factor(status)), size = 1)+
  geom_line(data = newdata.p_ref, 
            aes(x = rev_day, y = mean, group = status,color = factor(status)), size = 1)+
  geom_ribbon(data = newdata.p_ref, 
              aes(x = rev_day, ymin = lwr, ymax=upr, group = status,fill = factor(status)), alpha = 0.2)+
  facet_grid(biomarker~status,scales = "free_y")+ #, labeller = labeller(status = facet.lab)
  scale_color_manual(values = c("#00BFC4", "#F8766D","#7CAE00"))+
  scale_fill_manual(values = c("#00BFC4", "#F8766D","#7CAE00"))+
  theme_classic()+
  guides(fill = "none", color  = "none")+
  xlab("Day Until Event")+
  scale_y_continuous("", position = "right")+
  theme(axis.title = element_text(size = fontsize),
        axis.text = element_text(size = fontsize),
        strip.text = element_blank(),
        axis.text.y = element_blank())
p_retro

## forward in time -----
load("paperPrep/finalizedCode/model/alldata_jtgaussian_model_forw.RData")

fit.pre.forw <- fit.pre.forw
Bhat <- fixef(fit.pre.forw, use = "mean")
rn.Bhat <- rownames(Bhat)
order.Bhat <- c(grep("sao2_norm", rn.Bhat),
                grep("temp_norm", rn.Bhat),
                grep("pulse_norm", rn.Bhat))
Bhat <- matrix(Bhat[order.Bhat], ncol = 1)
rownames(Bhat) <- rn.Bhat[order.Bhat]

smy <- summary(fit.pre.forw)
smy.sol <- smy$solutions
nms <- rownames(smy.sol)
order.sol <- c(grep("sao2_norm", nms),
               grep("temp_norm", nms),
               grep("pulse_norm", nms))
sol <- (smy.sol[order.sol,])

grp <- 3
dt.baseline <- dt.pre %>% dplyr::select(resp_bl, temp_bl, pulse_bl, sao2fio2_bl, crp_bl, alc_bl,
                                        dd_bl, gfr_bl, charlson_cat, demo5, sex, bmi_cat, smoking) %>% unique()

newdata <- data.frame(osler_id = "Patient",
                      interval = rep(0:20, grp),
                      rev_day = rep(-20:0, grp),
                      status = rep(1:3, each  = 21),
                      resp_bl = mean(dt.baseline$resp_bl),
                      temp_bl = mean(dt.baseline$temp_bl),
                      pulse_bl = mean(dt.baseline$pulse_bl),
                      crp_bl = mean(dt.baseline$crp_bl),
                      alc_bl = mean(dt.baseline$alc_bl),
                      sao2fio2_bl = mean(dt.baseline$sao2fio2_bl),
                      dd_bl = mean(dt.baseline$dd_bl),
                      gfr_bl = 2,
                      charlson_cat = "0",
                      demo5 = "Young_White",
                      sex = "Female",
                      bmi_cat = "<30",
                      smoking ="Never Smoker",
                      sao2_fio2_ratio_m = NA,
                      temp_c_m = NA,
                      pulse_m = NA)
newdata <- newdata %>%  mutate(sex = factor(sex, levels = c("Female", "Male")),
                               smoking = factor(smoking, levels = c("Never Smoker", "Former Smoker", "Current Smoker")),
                               charlson_cat = factor(charlson_cat, levels = c("0", "1-2" , "3-4", ">5")),
                               demo5 = factor(demo5, levels =   c("Young_White", "Old_White" , "Young_Black", "Old_Black", "Latinx/Other")),
                               bmi_cat = factor(bmi_cat, levels = c("<30", ">=30")),
                               status = factor(status, levels = c("1", "2", "3")))


fm.Xi <- as.formula(paste0("~ns(interval,Boundary.knots = c(0,20), knots = c(1,3,6)) +
                                           resp_bl + temp_bl + pulse_bl + sao2fio2_bl +
                                           crp_bl + alc_bl + dd_bl + gfr_bl +
                                           charlson_cat + demo5 + sex + bmi_cat +
                                           smoking"))
Xi1 <- model.matrix(fm.Xi, data = newdata) #;colnames(Xi1) <- rownames(Bhat)[grep("sao2_fio2_ratio_m", rownames(Bhat))]
Xi2 <- model.matrix(fm.Xi, data = newdata)
Xi3 <- model.matrix(fm.Xi, data = newdata)
Xi <- as.matrix(bdiag(list(Xi1, Xi2, Xi3)))#Xi <- as.matrix(bdiag(list(Xi1, Xi3)))#

Yhati <-  Xi %*% Bhat
newdata$sao_fitted = Yhati[1:nrow(Xi1)]
newdata$temp_fitted = Yhati[(nrow(Xi1)+1):(2*nrow(Xi1))]
newdata$pulse_fitted = Yhati[(2*nrow(Xi1)+1):(3*nrow(Xi1))]

v.samp <- fit.pre.forw$Sol[,1:length(Bhat)] # 600 by 111
nms <- colnames(v.samp)
order.v <- c(grep("sao2_norm", nms),
             grep("temp_norm", nms),
             grep("pulse_norm", nms))
v.samp <- v.samp[order.v,order.v]
varFix <- var(v.samp)


predvar <- diag(Xi %*% varFix %*% t(Xi))
predse <- sqrt(predvar)
newdata$sao_se = predse[1:nrow(Xi1)]
newdata$temp_se = predse[(nrow(Xi1)+1):(2*nrow(Xi1))]
newdata$pulse_se = predse[(2*nrow(Xi1)+1):(3*nrow(Xi1))]


## get observations backwards in time
dt.joint <- readRDS("data/data_for_joint.rds")
tmp <- dt.joint %>% group_by(osler_id) %>% dplyr::summarize(maxstatus = max(status))
dt <- dt.joint %>% left_join(tmp, by = 'osler_id') %>% filter(complete.cases(status))
load("S:/COVID_LDA/code/zwang_mgbowring/data/vitals_data0809.RData")
dt.bls <- vitals_merge %>% dplyr::select(osler_id, status,first_vented, outcome) %>% unique()
dt <- dt %>% left_join(dt.bls)

dtp <- dt %>% filter(status %in% c(1,2,3))
tmp <- dtp %>% dplyr::group_by(osler_id) %>% dplyr::summarize(endpoint = max(interval)) %>% mutate(interval = endpoint)
dtp <- dtp %>% left_join(tmp)

dtp <- dtp %>% dplyr::filter(interval <= 15) %>% 
  mutate(max_interval_trim = ifelse(endpoint > 15, NA, endpoint))

## sample 50 persons
idlist <- unique(dtp$osler_id[!is.na(dtp$max_interval_trim)])
set.seed(2022)
ids <- dtp %>% 
  dplyr::select(osler_id, status) %>% 
  unique() %>% 
  group_by(status) %>% sample_n(100)
dtps <- dtp %>% filter(osler_id %in% ids$osler_id)

dtps <- dtps %>% filter(complete.cases(sao2_fio2_ratio_m, pulse_m, temp_c_m))
dtps$sao2_norm <- predict(normobj.sao2, dtps$sao2_fio2_ratio_m)
dtps$temp_norm <- predict(normobj.temp, dtps$temp_c_m)
dtps$pulse_norm <- predict(normobj.pulse, dtps$pulse_m)

facet.lab <- c("Discharge", "Die", "Ventilation")
names(facet.lab) <- c("1", "2", "3")

newdata.p <- newdata %>% filter(interval <= 15)

newdata.p$sao2_inv <- predict(normobj.sao2, newdata = newdata.p$sao_fitted, inverse = T)
newdata.p$temp_inv <- predict(normobj.temp, newdata = newdata.p$temp_fitted, inverse = T)
newdata.p$pulse_inv <- predict(normobj.pulse, newdata = newdata.p$pulse_fitted, inverse = T)

newdata.p <- newdata.p %>% mutate(sao2_lwr = sao_fitted - 1.96*sao_se,
                                  sao2_upr = sao_fitted + 1.96*sao_se,
                                  temp_lwr = temp_fitted - 1.96*temp_se,
                                  temp_upr = temp_fitted + 1.96*temp_se,
                                  pulse_lwr = pulse_fitted - 1.96*pulse_se,
                                  pulse_upr = pulse_fitted + 1.96*pulse_se)
newdata.p$sao2_upr_inv <- predict(normobj.sao2, newdata = newdata.p$sao2_upr, inverse = T)
newdata.p$sao2_lwr_inv <- predict(normobj.sao2, newdata = newdata.p$sao2_lwr, inverse = T)
newdata.p$pulse_upr_inv <- predict(normobj.pulse, newdata = newdata.p$pulse_upr, inverse = T)
newdata.p$pulse_lwr_inv <- predict(normobj.pulse, newdata = newdata.p$pulse_lwr, inverse = T)
newdata.p$temp_upr_inv <- predict(normobj.temp, newdata = newdata.p$temp_upr, inverse = T)
newdata.p$temp_lwr_inv <- predict(normobj.temp, newdata = newdata.p$temp_lwr, inverse = T)

dtps_sao2 <- dtps[,c("osler_id", "interval", "endpoint", "status", "sao2_fio2_ratio_m")]; dtps_sao2$biomarker = "sao2"
dtps_temp <- dtps[,c("osler_id", "interval", "endpoint", "status", "temp_c_m")];dtps_temp$biomarker = "temp"
dtps_pulse <- dtps[,c("osler_id", "interval", "endpoint", "status", "pulse_m")];dtps_pulse$biomarker = "pulse"
colnames(dtps_sao2) <- colnames(dtps_temp) <- colnames(dtps_pulse) <-
  c("osler_id", "interval", "endpoint", "status", "observed", "biomarker")
dtps_ref_pro <- rbind(dtps_sao2, dtps_temp, dtps_pulse)
dtps_ref_pro$biomarker = factor(dtps_ref_pro$biomarker, levels = c("sao2", "temp", "pulse"))

newdata_sao2 <- newdata.p[,c("osler_id", "interval", "status","sao2_inv", "sao2_lwr_inv", "sao2_upr_inv")];
newdata_sao2$biomarker = "sao2"
newdata_temp <- newdata.p[,c("osler_id", "interval", "status","temp_inv", "temp_lwr_inv", "temp_upr_inv")];
newdata_temp$biomarker = "temp"
newdata_pulse <- newdata.p[,c("osler_id", "interval", "status","pulse_inv", "pulse_lwr_inv", "pulse_upr_inv")];
newdata_pulse$biomarker = "pulse"
colnames(newdata_sao2) <- colnames(newdata_temp) <- colnames(newdata_pulse) <- 
  c("osler_id", "interval", "status", "mean","lwr","upr", "biomarker")
newdata.p_ref_pro <- rbind(newdata_sao2, newdata_temp, newdata_pulse)
newdata.p_ref_pro$biomarker = factor(newdata.p_ref_pro$biomarker, levels = c("sao2", "temp", "pulse"))

fontsize = 16
facet.lab <- c("SpO2/FiO2", "Temperature", "Pulse")
names(facet.lab) <- c("sao2", "temp", "pulse")

p_prosp <-  ggplot()+
  geom_line(data = dtps_ref_pro, aes(x = interval, y = observed, group = osler_id), color = "grey")+
  geom_jitter(data = dtps_ref_pro, aes(x = endpoint, y = observed, group = osler_id, color = factor(status)), size = 1)+
  geom_line(data = newdata.p_ref_pro, aes(x = interval, y = mean, group = status),color = 'gold4', size = 1)+
  geom_ribbon(data = newdata.p_ref_pro, aes(x = interval, ymin = lwr, 
                                    ymax=upr, group = status),fill = 'gold4', alpha = 0.1)+
  facet_grid(biomarker~.,scales = "free_y", switch = "y", labeller = labeller(biomarker = facet.lab))+ #
  scale_color_manual(values = c("#00BFC4", "#F8766D","#7CAE00"))+
  scale_fill_manual(values = c("#00BFC4", "#F8766D","#7CAE00"))+
  theme_classic()+
  guides(fill = "none", color  = "none")+
  xlab("Day of Hospitalization")+
  ylab("")+
  theme(axis.title = element_text(size = fontsize),
        axis.text = element_text(size = fontsize),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size = fontsize, colour = "black"))
p_prosp

#grid.arrange(p_prosp, p_retro, nrow = 1)

# p1 <- ggplot()+
#   geom_line(data = dtps, aes(x = interval, y = sao2_fio2_ratio_m, group = osler_id), color = "grey")+
#   geom_jitter(data = dtps, aes(x = endpoint, y = sao2_fio2_ratio_m, group = osler_id, color = factor(status)), size = 1)+
#   geom_line(data = newdata.p, aes(x = interval, y = sao2_inv, group = status),color = 'gold4', size = 1)+
#   geom_ribbon(data = newdata.p, aes(x = interval, ymin = sao2_lwr_inv, 
#                                     ymax=sao2_upr_inv, group = status),fill = 'gold4', alpha = 0.1)+
#   #facet_grid(~status, labeller = labeller(status = facet.lab))+
#   scale_color_manual(values = c("#00BFC4", "#F8766D","#7CAE00"))+
#   scale_fill_manual(values = c("#00BFC4", "#F8766D","#7CAE00"))+
#   theme_classic()+
#   guides(fill = F, color  = F)
# p2 <- ggplot()+
#   geom_line(data = dtps, aes(x = interval, y = temp_c_m, group = osler_id), color = "grey")+
#   geom_jitter(data = dtps, aes(x = endpoint, y = temp_c_m, group = osler_id, color = factor(status)), size = 1)+
#   geom_line(data = newdata.p, aes(x = interval, y = temp_inv, group = status),color = 'gold4', size = 1)+
#   geom_ribbon(data = newdata.p, aes(x = interval, ymin = temp_lwr_inv, 
#                                     ymax=temp_upr_inv, group = status),fill = 'gold4', alpha = 0.1)+
#   #facet_grid(~status, labeller = labeller(status = facet.lab))+
#   scale_color_manual(values = c("#00BFC4", "#F8766D","#7CAE00"))+
#   scale_fill_manual(values = c("#00BFC4", "#F8766D","#7CAE00"))+
#   theme_classic()+
#   guides(fill = F, color  = F)
# p3 <- ggplot()+
#   geom_line(data = dtps, aes(x = interval, y = pulse_m, group = osler_id), color = "grey")+
#   geom_jitter(data = dtps, aes(x = endpoint, y = pulse_m, group = osler_id, color = factor(status)), size = 1)+
#   geom_line(data = newdata.p, aes(x = interval, y = pulse_inv, group = status),color = 'gold4', size = 1)+
#   geom_ribbon(data = newdata.p, aes(x = interval, ymin = pulse_lwr_inv, 
#                                     ymax=pulse_upr_inv, group = status),fill = 'gold4', alpha = 0.1)+
#   #facet_grid(~status, labeller = labeller(status = facet.lab))+
#   scale_color_manual(values = c("#00BFC4", "#F8766D","#7CAE00"))+
#   scale_fill_manual(values = c("#00BFC4", "#F8766D","#7CAE00"))+
#   theme_classic()+
#   guides(fill = F, color  = F)
# pdf("paperPrep/finalizedCode/figs/jtgaussian_fig_prosp.pdf", width = 4, height = 8)
# grid.arrange(p1+ylab("SaO2/SpO2")+xlab("Day until event")+
#                theme(axis.title = element_text(size = 14),
#                      axis.text = element_text(size = 14),
#                      strip.text = element_text(size = 14)), 
#              p2+ylab("Temperature")+xlab("Day until event")+
#                theme(axis.title = element_text(size = 14),
#                      axis.text = element_text(size = 14),
#                      strip.text = element_text(size = 14)),
#              p3+ylab("Pulse")+xlab("Day until event")+
#                theme(axis.title = element_text(size = 14),
#                      axis.text = element_text(size = 14),
#                      strip.text = element_text(size = 14)), ncol = 1)
# dev.off()


toprow_prosp <- data.frame(type = c("                              Prospective"),
                           event = c(""))
toprow_retro <- data.frame(type = c( " ",  
                                     "Retrospective             ",
                                     "                 "),
                           event = c("            Discharge                     ", 
                                     "   Died                        ", 
                                     "Vent"))

toptheme <- ttheme_minimal(core = list(fg_params=list(cex=1.5)))
ptop_prosp<- grid.arrange(tableGrob(t(toprow_prosp), theme = toptheme, rows = NULL ,cols = NULL))
ptop_retro<- grid.arrange(tableGrob(t(toprow_retro), theme = toptheme, rows = NULL ,cols = NULL))
top_row = list(ptop_prosp, ptop_retro)
library(cowplot)
ptop <- plot_grid(plotlist = top_row, align = "h", axis  = "l", ncol = 2, nrow = 1,
          rel_widths = c(1,3))


## merge plots together ----
bottom_row <- list(p_prosp, p_retro)
pbottom <- plot_grid(plotlist = bottom_row, align = "h", axis  = "l", ncol = 2, nrow = 1,
                     rel_widths = c(1,2))
plist <- list(ptop, pbottom)
plot_grid(plotlist = plist, align="v", axis = "t", ncol = 1, nrow = 2,
          rel_heights = c(1,8))

pdf("paperPrep/finalizedCode/figs/jtgaussian_fig_merged.pdf", width = 12, height = 6)
plot_grid(plotlist = plist, align="v", axis = "t", ncol = 1, nrow = 2,
          rel_heights = c(1,8))
dev.off()

