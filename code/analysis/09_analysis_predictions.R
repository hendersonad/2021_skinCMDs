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




# load data and model -----------------------------------------------------
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

fit <- survfit(smpl_fit, newdata = data.frame(exposed = c("Unexposed", "Psoriasis")))
fitsum <- summary(fit)
fitsum$year <- fitsum$time/365.25

plot(fitsum$year, fitsum$surv[,1], type="l",ylim=c(0,1),col=c("black"),xlab="Year",ylab="Estimated survivor function")
lines(fitsum$year, fitsum$surv[,2], type="l",ylim=c(0,1),col=c("grey"),xlab="Year",ylab="Estimated survivor function")
legend(0,0.2,c("Unexposed",str_to_title(exposure)),col=c("black","grey"),lty=1,cex=1)

# estimate baseline hazrd -------------------------------------------------
vector_of_times <- df_sample$t[df_sample$out == 1]
sort(vector_of_times)

predict_df <- data.frame(t = sort(vector_of_times))

predict_df <- predict_df %>% 
  group_by(t) %>% 
  summarise(dj = n())

test <- predict_df %>% 
  mutate()

h0t <- function(
  dj <- sum()
)

# specify the survival scenarios ------------------------------------------
## choose confounder covariates
## Year 2010, age 30, white, IMD 3, no comorbidity

## then add in mediators - alcohol/smok/BMI/sleep/gc90





