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
outcome <- YY[1]
ABBRVexp <- substr(exposure, 1, 3)

# load df_model -----------------------------------------------------------
df_model <-
  readRDS(paste0(datapath, "out/models_data/df_model", ABBRVexp, "_", outcome,".rds"))

cox_fit <- readRDS(paste0(datapath, "out/models_data/", ABBRVexp, "_", outcome, "_mod1_modeldata.rds"))
cox_fit3 <- readRDS(paste0(datapath, "out/models_data/", ABBRVexp, "_", outcome, "_mod3_modeldata.rds"))

set.seed(30)
sample_ids <- sample(unique(df_model$setid), size = 10000, replace = F)
df_sample <- df_model %>% 
  filter(setid %in% sample_ids)

smpl_fit <- coxph(Surv(t, out) ~ exposed + strata(setid), data = df_sample)
smpl_fit <- coxph(Surv(t, out) ~ exposed, data = df_sample)
km_fit <- survfit(Surv(t, out) ~ exposed + carstairs, data = df_sample)
km_fit_full <- survfit(Surv(t, out) ~ exposed, data = df_sample)

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


plot(km_fit,fun="cloglog",xlab="time (log scale)",ylab="log(-log S(t))",
     col=c("blue","red"))
legend("topleft",c("Placebo","Active"),col=c("blue","red"),lty=1)


zph <- cox.zph(cox_fit)
szph <- cox.zph(strat_cox_fit)
par(mfrow = c(1, 1))
plot(zph, var = 1)


fit <- survfit(smpl_fit, newdata = data.frame(exposed = c("Unexposed", "Psoriasis")))

plot(fit, fun = "cloglog", ,col=c("black","grey"),xlab="time",ylab="Estimated survivor function")

