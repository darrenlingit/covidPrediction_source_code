rm(list=ls())
nbiom <- 3
df.re <- 3
nb = 2000; nt = 8000
source("00_loadlibs.R")
setwd(wd)

#load data
dt.pre <- data
fit.pre = MCMCglmm(cbind(sao2_norm, temp_norm, pulse_norm) ~
                     -1 + trait + trait:(status*ns(rev_day,knots = c(-7,-4,-2,-1), Boundary.knots = c(-20,0)) +
                                           resp_bl + temp_bl + pulse_bl + sao2fio2_bl +
                                           crp_bl + alc_bl + dd_bl + gfr_bl +
                                           charlson_cat + demo5 + sex + bmi_cat +
                                           smoking),
                   random = ~ us(trait + trait:ns(rev_day,knots = c(-3), Boundary.knots = c(-20,0))):osler_id,  #b0 + b1day for each outcome/trait
                   rcov = ~ us(trait):units, #unstructured, coursenote pg67
                   burnin = nb, nitt = nt, pr = T,pl=T,
                   family = rep("gaussian", nbiom),
                   data = dt.pre %>% filter(interval!=0))
save(fit.pre, file = "model/alldata_jtgaussian_model.RData")

set.seed(2021)
fit.pre.forw = MCMCglmm(cbind(sao2_norm, temp_norm, pulse_norm) ~
                     -1 + trait + trait:(ns(interval,Boundary.knots = c(0,20), knots = c(1,3,6)) +
                                           resp_bl + temp_bl + pulse_bl + sao2fio2_bl +
                                           crp_bl + alc_bl + dd_bl + gfr_bl +
                                           charlson_cat + demo5 + sex + bmi_cat +
                                           smoking),
                   random = ~ us(trait + trait:interval):osler_id,  #b0 + b1day for each outcome/trait
                   rcov = ~ us(trait):units, #unstructured, coursenote pg67
                   burnin = nb, nitt = nt, pr = T,
                   family = rep("gaussian", 3),
                   data = dt.pre %>% filter(interval != 0))
save(fit.pre.forw, file = "model/alldata_jtgaussian_model_forw.RData")

## residual vs fitted; backward in time-----
load("model/alldata_jtgaussian_model.RData")
normobj.sao2 <- orderNorm(dt.pre$sao2_fio2_ratio_o)
normobj.temp <- orderNorm(dt.pre$temp_c_o)
normobj.pulse <- orderNorm(dt.pre$pulse_o)

normobj <- list(normobj.sao2, normobj.temp, normobj.pulse)
save(normobj, file = "model/normobj.RData")

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
load("model/alldata_jtgaussian_model_forw.RData")
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

grid.arrange(pb, pf, ncol = 1)

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
