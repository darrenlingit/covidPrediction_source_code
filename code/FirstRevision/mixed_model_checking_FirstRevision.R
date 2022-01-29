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

## residual vs fitted; backward in time-----
load("paperPrep/finalizedCode/model/alldata_jtgaussian_model.RData")
normobj.sao2 <- orderNorm(dt.pre$sao2_fio2_ratio_o)
normobj.temp <- orderNorm(dt.pre$temp_c_o)
normobj.pulse <- orderNorm(dt.pre$pulse_o)

# normobj <- list(normobj.sao2, normobj.temp, normobj.pulse)
# save(normobj, file = "paperPrep/finalizedCode/model/normobj.RData")

fit.pre <- fit.pre
library(broom.mixed)
Bhat <- fixef(fit.pre, use = "mean")
rn.Bhat <- rownames(Bhat)
order.Bhat <- c(grep("sao2_norm", rn.Bhat),
                grep("temp_norm", rn.Bhat),
                grep("pulse_norm", rn.Bhat))
Bhat <- matrix(Bhat[order.Bhat], ncol = 1)
rownames(Bhat) <- rn.Bhat[order.Bhat]

newdata = dt.pre%>% filter(interval!=0)
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

Yi <- c(newdata$sao2_norm, newdata$temp_norm, newdata$pulse_norm)
Residi <- Yi - Yhati
newdata$sao_resid = Residi[1:nrow(Xi1)]
newdata$temp_resid = Residi[(nrow(Xi1)+1):(2*nrow(Xi1))]
newdata$pulse_resid = Residi[(2*nrow(Xi1)+1):(3*nrow(Xi1))]
library(Hmisc)
quantiles <- quantile(Yhati, prob = seq(0,1, length = 11))
newdata$sao_fitted_decile <- cut2(newdata$sao_fitted, cuts = as.numeric(quantiles))
newdata$temp_fitted_decile <- cut2(newdata$temp_fitted, cuts = as.numeric(quantiles))
newdata$pulse_fitted_decile <- cut2(newdata$pulse_fitted, cuts = as.numeric(quantiles))
# newdata_saosumm <- newdata %>% 
#   group_by(sao_fitted_decile) %>% 
#   summarise(meanres_sao = mean(sao_resid),
#             sdres_sao = sd(sao_resid),
#             nres_sao = n())
# newdata_tempsumm <- newdata %>% 
#   group_by(temp_fitted_decile) %>% 
#   summarise(meanres_temp = mean(temp_resid),
#             sdres_temp = sd(temp_resid),
#             nres_temp = n())
# newdata_pulsesumm <- newdata %>% 
#   group_by(pulse_fitted_decile) %>% 
#   summarise(meanres_pulse = mean(pulse_resid),
#             sdres_pulse = sd(pulse_resid),
#             nres_pulse = n())
# 
# pb1 <- ggplot() + 
#   geom_pointrange(data = newdata_saosumm, aes(x=sao_fitted_decile, 
#                                               y = meanres_sao,
#                                               ymin = meanres_sao - 2*sdres_sao/sqrt(nres_sao),
#                                               ymax = meanres_sao + 2*sdres_sao/sqrt(nres_sao)))+
#   theme_classic()+
#   xlab("Fitted Values in Deciles") + 
#   ylab("Residuals") + 
#   ggtitle("SpO2/FiO2")+
#   geom_hline(yintercept = 0)
# pb2 <- ggplot() + 
#   geom_pointrange(data = newdata_tempsumm, aes(x=temp_fitted_decile, 
#                                               y = meanres_temp,
#                                               ymin = meanres_temp - 2*sdres_temp/sqrt(nres_temp),
#                                               ymax = meanres_temp + 2*sdres_temp/sqrt(nres_temp)))+
#   theme_classic()+
#   xlab("Fitted Values in Deciles") + 
#   ylab("Residuals") + 
#   ggtitle("Temperature")+
#   geom_hline(yintercept = 0)
# pb3 <- ggplot() + 
#   geom_pointrange(data = newdata_pulsesumm, aes(x=pulse_fitted_decile, 
#                                                y = meanres_pulse,
#                                                ymin = meanres_pulse - 2*sdres_pulse/sqrt(nres_pulse),
#                                                ymax = meanres_pulse + 2*sdres_pulse/sqrt(nres_pulse)))+
#   theme_classic()+
#   xlab("Fitted Values in Deciles") + 
#   ylab("Residuals") + 
#   ggtitle("Pulse")+
#   geom_hline(yintercept = 0)
# pb <- grid.arrange(pb1+ylim(c(-0.5,0.2)), 
#                    pb2+ylim(c(-0.5,0.2)), 
#                    pb3+ylim(c(-0.5,0.2)), nrow = 1, top = "Retrospective Linear Mixed Model")
xbreaks <- levels(newdata$sao_fitted_decile)
xlabs <- round(quantiles[-11],2)
names(xlabs) <- NULL
pb1 <- ggplot() + 
  geom_boxplot(data = newdata, aes(x=sao_fitted_decile, y = sao_resid))+
  theme_classic()+
  xlab("Fitted Values in Deciles") + 
  ylab("Residuals") + 
  ggtitle("SpO2/FiO2")+
  geom_hline(yintercept = 0)+
  scale_x_discrete(breaks = xbreaks, 
                   labels = xlabs)
pb2 <- ggplot() + 
  geom_boxplot(data = newdata, aes(x=temp_fitted_decile, y = temp_resid))+
  theme_classic()+
  xlab("Fitted Values in Deciles") + 
  ylab("Residuals") + 
  ggtitle("Temperature")+
  geom_hline(yintercept = 0)+
  scale_x_discrete(breaks = xbreaks, 
                   labels = xlabs)
pb3 <- ggplot() + 
  geom_boxplot(data = newdata, aes(x=pulse_fitted_decile, y = pulse_resid))+
  theme_classic()+
  xlab("Fitted Values in Deciles") + 
  ylab("Residuals") + 
  ggtitle("Pulse")+
  geom_hline(yintercept = 0)+
  scale_x_discrete(breaks = xbreaks, 
                   labels = xlabs)
pb <- grid.arrange(pb1, 
                   pb2, 
                   pb3, nrow = 1, top = "Retrospective Linear Mixed Model")





## residual vs fitted; forward in time -----
load("paperPrep/finalizedCode/model/alldata_jtgaussian_model_forw.RData")
fit.pre.forw <- fit.pre.forw
Bhat <- fixef(fit.pre.forw, use = "mean")
rn.Bhat <- rownames(Bhat)
order.Bhat <- c(grep("sao2_norm", rn.Bhat),
                grep("temp_norm", rn.Bhat),
                grep("pulse_norm", rn.Bhat))
Bhat <- matrix(Bhat[order.Bhat], ncol = 1)
rownames(Bhat) <- rn.Bhat[order.Bhat]

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
newdata$sao_fitted_f = Yhati[1:nrow(Xi1)]
newdata$temp_fitted_f = Yhati[(nrow(Xi1)+1):(2*nrow(Xi1))]
newdata$pulse_fitted_f = Yhati[(2*nrow(Xi1)+1):(3*nrow(Xi1))]

Yi <- c(newdata$sao2_norm, newdata$temp_norm, newdata$pulse_norm)
Residi <- Yi - Yhati
newdata$sao_resid_f = Residi[1:nrow(Xi1)]
newdata$temp_resid_f = Residi[(nrow(Xi1)+1):(2*nrow(Xi1))]
newdata$pulse_resid_f = Residi[(2*nrow(Xi1)+1):(3*nrow(Xi1))]

quantiles <- quantile(Yhati, prob = seq(0,1, length = 11))
newdata$sao_fitted_f_decile <- cut2(newdata$sao_fitted_f, cuts = as.numeric(quantiles))
newdata$temp_fitted_f_decile <- cut2(newdata$temp_fitted_f, cuts = as.numeric(quantiles))
newdata$pulse_fitted_f_decile <- cut2(newdata$pulse_fitted_f, cuts = as.numeric(quantiles))

# newdata_saosumm_f <- newdata %>% 
#   group_by(sao_fitted_f_decile) %>% 
#   summarise(meanres_sao = mean(sao_resid_f),
#             sdres_sao = sd(sao_resid_f),
#             nres_sao = n())
# newdata_tempsumm_f <- newdata %>% 
#   group_by(temp_fitted_f_decile) %>% 
#   summarise(meanres_temp = mean(temp_resid_f),
#             sdres_temp = sd(temp_resid_f),
#             nres_temp = n())
# newdata_pulsesumm_f <- newdata %>% 
#   group_by(pulse_fitted_f_decile) %>% 
#   summarise(meanres_pulse = mean(pulse_resid_f),
#             sdres_pulse = sd(pulse_resid_f),
#             nres_pulse = n())
# 
# 
# pf1 <- ggplot() + 
#   geom_pointrange(data = newdata_saosumm_f, aes(x=sao_fitted_f_decile, 
#                                               y = meanres_sao,
#                                               ymin = meanres_sao - 2*sdres_sao/sqrt(nres_sao),
#                                               ymax = meanres_sao + 2*sdres_sao/sqrt(nres_sao)))+
#   theme_classic()+
#   xlab("Fitted Values in Deciles") + 
#   ylab("Residuals") + 
#   ggtitle("SpO2/FiO2")+
#   geom_hline(yintercept = 0)
# pf2 <- ggplot() + 
#   geom_pointrange(data = newdata_tempsumm_f, aes(x=temp_fitted_f_decile, 
#                                                y = meanres_temp,
#                                                ymin = meanres_temp - 2*sdres_temp/sqrt(nres_temp),
#                                                ymax = meanres_temp + 2*sdres_temp/sqrt(nres_temp)))+
#   theme_classic()+
#   xlab("Fitted Values in Deciles") + 
#   ylab("Residuals") + 
#   ggtitle("Temperature")+
#   geom_hline(yintercept = 0)
# pf3 <- ggplot() + 
#   geom_pointrange(data = newdata_pulsesumm_f, aes(x=pulse_fitted_f_decile, 
#                                                 y = meanres_pulse,
#                                                 ymin = meanres_pulse - 2*sdres_pulse/sqrt(nres_pulse),
#                                                 ymax = meanres_pulse + 2*sdres_pulse/sqrt(nres_pulse)))+
#   theme_classic()+
#   xlab("Fitted Values in Deciles") + 
#   ylab("Residuals") + 
#   ggtitle("Pulse")+
#   geom_hline(yintercept = 0)
# pf <- grid.arrange(pf1+ylim(c(-0.9,0.7)), 
#                    pf2+ylim(c(-0.9,0.7)),
#                    pf3+ylim(c(-0.9,0.7)), nrow = 1, top = "Prospective Linear Mixed Model")
xbreaks <- levels(newdata$sao_fitted_f_decile)
xlabs <- round(quantiles[-11],2)
names(xlabs) <- NULL
pf1 <- ggplot() + 
  geom_boxplot(data = newdata, aes(x=sao_fitted_f_decile, y = sao_resid_f))+
  theme_classic()+
  xlab("Fitted Values in Deciles") + 
  ylab("Residuals") + 
  ggtitle("Pulse")+
  geom_hline(yintercept = 0)+
  scale_x_discrete(breaks = xbreaks, 
                   labels = xlabs)
pf2 <- ggplot() + 
  geom_boxplot(data = newdata, aes(x=temp_fitted_f_decile, y = temp_resid_f))+
  theme_classic()+
  xlab("Fitted Values in Deciles") + 
  ylab("Residuals") + 
  ggtitle("Temperature")+
  geom_hline(yintercept = 0)+
  scale_x_discrete(breaks = xbreaks, 
                   labels = xlabs)
pf3 <- ggplot() + 
  geom_boxplot(data = newdata, aes(x=pulse_fitted_f_decile, y = pulse_resid_f))+
  theme_classic()+
  xlab("Fitted Values in Deciles") + 
  ylab("Residuals") + 
  ggtitle("Pulse")+
  geom_hline(yintercept = 0)+  
  scale_x_discrete(breaks = xbreaks, 
                   labels = xlabs)

pf <- grid.arrange(pf1, 
                   pf2,
                   pf3, nrow = 1, top = "Prospective Linear Mixed Model")

pdf("paperPrep/finalizedCode/figs/resid_fitted_hist.pdf")
grid.arrange(pb, pf, ncol = 1)
dev.off()

pdf("paperPrep/finalizedCode/figs/resid_qq.pdf")
par(mfrow = c(2,3))
qqnorm(newdata$sao_resid, main = "SpO2/FiO2")
abline(0,1, col = "blue")

qqnorm(newdata$temp_resid, main = "Temperature")
abline(0,1, col = "blue")

qqnorm(newdata$pulse_resid, main = "Pulse")
abline(0,1, col = "blue")

mtext("Retrospective LMM Residuals Normal Q-Q Plot", outer = TRUE,
      side = 3, line = -1)

qqnorm(newdata$sao_resid_f, main = "SpO2/FiO2")
abline(0,1, col = "blue")

qqnorm(newdata$temp_resid_f, main = "Temperature")
abline(0,1, col = "blue")

qqnorm(newdata$pulse_resid_f, main = "Pulse")
abline(0,1, col = "blue")

mtext("Prospective LMM Residuals Normal Q-Q Plot", outer = TRUE,
      side = 3, line = -28)
dev.off()
