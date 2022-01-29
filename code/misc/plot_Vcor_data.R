## Residual auto-correlation plot for Y-XB (data based Vcor)
## compare model based to data based Vcor 

# residual auto-correlation plot Y-XB -----------
rm(list=ls())
nbiom <- 3
redf <- 3
dt.joint <- readRDS("data/data_for_joint.rds")
## create max status
tmp <- dt.joint %>% group_by(osler_id) %>% dplyr::summarize(maxstatus = max(status))
dt <- dt.joint %>% left_join(tmp, by = 'osler_id') %>% filter(complete.cases(status))
load("S:/COVID_LDA/code/zwang_mgbowring/data/vitals_data0809.RData")
dt.bls <- vitals_merge %>% dplyr::select(osler_id, status,first_vented, outcome, max_interval) %>% unique()
dt <- dt %>% left_join(dt.bls)

dt1 <- dt %>% filter(status %in% c(1,2,3))
dt1$long <- ifelse(dt1$rev_day < -20,'long','short')
id.check <- dt1 %>% dplyr::select(osler_id, long, status) %>% filter(long == "long") %>% unique()
dt.pre <- dt1 %>% dplyr::filter(!(osler_id %in% id.check$osler_id)) #length(unique(dt.pre$osler_id)) 


dt.pre <- dt.pre %>% 
  dplyr::rename(sao2_fio2_ratio_o = sao2_fio2_ratio_m,
                temp_c_o = temp_c_m,
                pulse_o = pulse_m)

library(bestNormalize)
dt.pre <- dt.pre %>% 
  filter(complete.cases(sao2_fio2_ratio_o, temp_c_o, pulse_o))  ## 1687
## build separate lme model for each biomarker
# attr(ns(dt.pre$rev_day,4), "knots"); attr(ns(dt.pre$rev_day,4), "Boundary.knots")
# attr(ns(dt.pre$rev_day,2), "knots"); attr(ns(dt.pre$rev_day,2), "Boundary.knots")
# get.fm.pre <-function(var) as.formula(paste0(var, "~ status * ns(rev_day,knots = c(-6,-3,-1),
#                                  Boundary.knots = c(-20,0)) + resp_bl + temp_bl +
#                                  pulse_bl + sao2fio2_bl + crp_bl + alc_bl + dd_bl + gfr_bl +
#                                  charlson_cat + demo5 + sex + bmi_cat + smoking"))
# lme1 <- lme(get.fm.pre("sao2_fio2_ratio_m"), random = ~1+ns(rev_day,knots = c(-3), Boundary.knots = c(-20,0))|osler_id,data = dt.pre, na.action = na.exclude)
# lme3 <- lme(get.fm.pre("temp_c_m"), random = ~1+ns(rev_day,knots = c(-3), Boundary.knots = c(-20,0))|osler_id,data = dt.pre, na.action = na.exclude)
# lme4 <- lme(get.fm.pre("pulse_m"), random = ~1+ns(rev_day,knots = c(-3), Boundary.knots = c(-20,0))|osler_id,data = dt.pre, na.action = na.exclude)
# fit.pre <- list(lme1,lme3, lme4)
# save(fit.pre, file = "paperPrep/fit.lme.pre.RData")
# 


load("paperPrep/model/fit.lme.pre.RData")
load("paperPrep/model/fit.lme.post.RData")
dt1 <- dt %>% filter(status %in% c(1,2,3))
dt1$long <- ifelse(dt1$rev_day < -20,'long','short')
id.check <- dt1 %>% dplyr::select(osler_id, long, status) %>% filter(long == "long") %>% unique()
dt.pre <- dt1 %>% dplyr::filter(!(osler_id %in% id.check$osler_id)) #length(unique(dt.pre$osler_id)) #768 patients

load("paperPrep/model/fit.lme.pre.RData")
dt.pre$res1 <- residuals(fit.pre[[1]], level = 0);
dt.pre$res2 <- residuals(fit.pre[[2]], level = 0);
dt.pre$res3 <- residuals(fit.pre[[3]], level = 0)
dtcorr1.pre <- dt.pre %>%  
  filter(complete.cases(sao2_fio2_ratio_m)) %>%
  group_by(osler_id) %>% 
  dplyr::select(osler_id, rev_day, res1) %>% 
  distinct(.keep_all = TRUE) %>% 
  spread(rev_day, res1) 
plot_corr1.pre <- cor(as.matrix(dtcorr1.pre[,-1]), use = 'pairwise.complete.obs')
#corrplot(plot_corr1.pre)
dtcorr2.pre <- dt.pre %>%  
  filter(complete.cases(temp_c_m)) %>%
  group_by(osler_id) %>% 
  dplyr::select(osler_id, rev_day, res2) %>% 
  distinct(.keep_all = TRUE) %>% 
  spread(rev_day, res2) 
plot_corr2.pre <- cor(as.matrix(dtcorr2.pre[,-1]), use = 'pairwise.complete.obs')

dtcorr3.pre <- dt.pre %>% 
  filter(complete.cases(pulse_m)) %>%
  group_by(osler_id) %>% 
  dplyr::select(osler_id, rev_day, res3) %>% 
  distinct(.keep_all = TRUE) %>% 
  spread(rev_day, res3)
plot_corr3.pre <- cor(as.matrix(dtcorr3.pre[,-1]), use = 'pairwise.complete.obs')

##jointly
coln <- c("osler_id", -15:0)
dtcorr1.pre <- dtcorr1.pre[, coln]
dtcorr2.pre <- dtcorr2.pre[, coln]
dtcorr3.pre <- dtcorr3.pre[, coln]
data_corr_all <- left_join(dtcorr1.pre, dtcorr2.pre, by = "osler_id") %>%
  left_join(dtcorr3.pre, by = "osler_id") 
plot_corr <- cor(as.matrix(data_corr_all[,-1]), use = "pairwise.complete.obs")
rownames(plot_corr) <- colnames(plot_corr) <- rep(c(-15:0), 3)
pdf("paperPrep/finalizedCode/figs/data_Vcor.pdf")
corrplot(plot_corr,  title = "", tl.col = "black", tl.cex = 0.6)
dev.off()
