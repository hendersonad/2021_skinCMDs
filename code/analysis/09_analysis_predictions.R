library(tidyverse)
library(here)
library(magrittr)
library(gt)
library(gtsummary)
library(survival)
library(readstata13)

if (Sys.info()["user"] == "lsh1510922") {
  if (Sys.info()["sysname"] == "Darwin") {
    #datapath <- "/Users/lsh1510922/Documents/Postdoc/2021_extract/"
    datapath <- "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
  if (Sys.info()["sysname"] == "Windows") {
    datapath <- "Z:/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
}

dir.create(file.path(here("out")), showWarnings = FALSE)


YY <- c("depression", "anxiety")
XX <- c("psoriasis", "eczema")

exposure <- XX[1]
outcome <- YY[2]
ABBRVexp <- substr(exposure, 1, 3)




# load data and model -----------------------------------------------------
df_model <-
  readRDS(paste0(datapath, "out/models_data/df_model", ABBRVexp, "_", outcome,".rds"))
cox_fit <- readRDS(paste0(datapath, "out/models_data/", ABBRVexp, "_", outcome, "_mod1_modeldata.rds"))
cox_fit3 <- readRDS(paste0(datapath, "out/models_data/", ABBRVexp, "_", outcome, "_mod3_modeldata.rds"))

set.seed(30)
sample_ids <- sample(unique(df_model$setid), size = 10000, replace = F)
df_sample <- df_model %>% 
  filter(setid %in% sample_ids)

smpl_fit_strat <- coxph(Surv(t, out) ~ exposed + strata(setid), data = df_sample)
smpl_fit <- coxph(Surv(t, out) ~ exposed, data = df_sample)

if (ABBRVexp == "ecz") {
  smpl_fit3 <-
    coxph(
      Surv(t, out) ~ exposed + carstairs + cal_period + comorbid + cci + bmi2 + sleep + alc + smokstatus + gc90days + strata(setid),
      data = df_sample
    ) 
} else if (ABBRVexp == "pso") {
  smpl_fit3 <-
    coxph(
      Surv(t, out) ~ exposed + carstairs + cal_period + comorbid + cci + bmi2 + alc + smokstatus + strata(setid),
      data = df_sample
    ) 
}
table(df_sample$exposed)
summary(cox_fit, conf.int = F)
summary(smpl_fit, conf.int = F)
summary(smpl_fit_strat, conf.int = F)



n_ids <- 100
subset_sets <- sample_ids[1:n_ids]
data_new_1 <- data.frame(exposed = c(rep("Psoriasis",n_ids)), setid = rep(subset_sets,1))
data_new_0 <- data.frame(exposed = c(rep("Unexposed",n_ids)), setid = rep(subset_sets,1))

fit_strat_1 <- survfit(cox_fit, newdata = data_new_1)
fit_strat_0 <- survfit(cox_fit, newdata = data_new_0)
dim(fit_strat_1)
dim(fit_strat_0)
summary(fit_strat_0)
fit_strat_0$strata
strata_ids0 <- data.frame(nrows = fit_strat_0$strata, id = subset_sets, exposed = 0) %>% 
  uncount(nrows)
strata_ids1 <- data.frame(nrows = fit_strat_1$strata, id = subset_sets, exposed = 1) %>% 
  uncount(nrows)
df_predicted0 <- data.frame(t = fit_strat_0$time, 
                           surv = fit_strat_0$surv, 
                           surv_lo = fit_strat_0$lower, 
                           surv_hi = fit_strat_0$upper)
df_predicted1 <- data.frame(t = fit_strat_1$time, 
                           surv = fit_strat_1$surv, 
                           surv_lo = fit_strat_1$lower, 
                           surv_hi = fit_strat_1$upper)


df_predicted_combo <- df_predicted0 %>% bind_cols(strata_ids0) %>% 
  bind_rows(
    bind_cols(df_predicted1, strata_ids1)
  )
df_predicted_combo$setid <- factor(df_predicted_combo$id)
df_predicted_combo$exposed <- factor(df_predicted_combo$exposed)
#df_sample %>% filter(setid == 22346748) %>% select(setid, patid, exposed, out, t)

#plot(fit_strat)

ggplot(df_predicted_combo, aes(x = t, y = surv, ymin = surv_lo, ymax = surv_hi, group = setid)) + 
  geom_ribbon(aes(fill = exposed), alpha = 0.02) +
  geom_line(aes(col = exposed), lwd = 0.5, alpha = 0.5) +
  facet_wrap(~exposed)

dev.copy(pdf, here("out/supplementary/survivor_curves_test.pdf"), width = 6, height = 6)
dev.off()

# specify the survival scenarios ------------------------------------------
## choose confounder covariates
## Year 2010, age 30, white, IMD 3, no comorbidity

## then add in mediators - alcohol/smok/BMI/sleep/gc90





