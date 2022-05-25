library(tidyverse)
library(here)
library(labelled)
library(stringr)
library(janitor)
library(timetk)
library(gtsummary)
library(magrittr)
library(survival)
library(survminer)
library(htmltools)
library(lubridate)
library(gt)

if (Sys.info()["user"] == "lsh1510922") {
  if (Sys.info()["sysname"] == "Darwin") {
    datapath <- "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
  if (Sys.info()["sysname"] == "Windows") {
    datapath <- "Z:/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
}
dir.create(file.path(here("out", "analysis")), showWarnings = FALSE)
dir.create(file.path(datapath, "out/models_data"), showWarnings = FALSE)

XX <- c("psoriasis", "eczema")
YY <- c("anxiety", "depression")


st_time <- Sys.time()
for (exposure in XX) {
  #exposure <- XX[2]
  ABBRVexp <- str_sub(exposure, 1 , 3)
  if (exposure == "eczema") {
    df_anx_split <-
      readRDS(paste0(datapath, "out/df_model", ABBRVexp, "_anxiety_noghosts_3yrs.rds"))
    df_dep_split <-
      readRDS(paste0(datapath, "out/df_model", ABBRVexp, "_depression_noghosts_3yrs.rds"))
  } else if (exposure == "psoriasis") {
    df_anx_split <-
      readRDS(paste0(datapath, "out/df_model", ABBRVexp, "_anxiety_noghosts_3yrs.rds"))
    df_dep_split <-
      readRDS(paste0(datapath, "out/df_model", ABBRVexp, "_depression_noghosts_3yrs.rds"))
  }
  
  .dib(exposure)
  
  for (outcome in YY) {
    #outcome = YY[2]
    .dib(outcome)
    
    if (outcome == "anxiety") {
      df_validsets_noghosts <- df_anx_split
    } else if (outcome == "depression") {
      df_validsets_noghosts <- df_dep_split
    }
    
    # Run a simple Cox regression ----------------------------------------
    .dib("Running model 1 (crude)")
    mod1 <-
      coxph(Surv(t, out) ~ exposed + strata(setid), data = df_validsets_noghosts)
    
    # Run a confounder Cox regression ----------------------------------------
    .dib("Running model 2 (confounder)")
    mod2 <-
      coxph(Surv(t, out) ~ exposed + carstairs + cal_period + comorbid + cci + strata(setid),
            data = df_validsets_noghosts)
    
    # Run a mediator Cox regression -------------------------------------------
    .dib("Running model 3 (mediator)")
    if (ABBRVexp == "ecz") {
      mod3 <-
        coxph(
          Surv(t, out) ~ exposed + carstairs + cal_period + comorbid + cci + bmi_cat + sleep + alc + smokstatus + gc90days + strata(setid),
          data = df_validsets_noghosts
        ) 
    } else if (ABBRVexp == "pso") {
      mod3 <-
        coxph(
          Surv(t, out) ~ exposed + carstairs + cal_period + comorbid + cci + bmi_cat + alc + smokstatus + strata(setid),
          data = df_validsets_noghosts
        ) 
    }
    
    saveRDS(
      mod1,
      file = paste0(
        datapath,
        "out/models_data/",
        ABBRVexp,
        "_",
        outcome,
        "_mod1_modeldata_noghosts-3yrs.rds"
      )
    )
    saveRDS(
      mod2,
      file = paste0(
        datapath,
        "out/models_data/",
        ABBRVexp,
        "_",
        outcome,
        "_mod2_modeldata_noghosts-3yrs.rds"
      )
    )
    saveRDS(
      mod3,
      file = paste0(
        datapath,
        "out/models_data/",
        ABBRVexp,
        "_",
        outcome,
        "_mod3_modeldata_noghosts-3yrs.rds"
      )
    )
  }
  
  if (sum(grepl(pattern = "df_model_noghosts", x = ls())) == 1) {
    rm(df_model)
  }
  if (sum(grepl(pattern = "df_dep_split", x = ls())) == 1) {
    rm(df_dep_split)
  }
  if (sum(grepl(pattern = "df_anx_split", x = ls())) == 1) {
    rm(df_anx_split)
  }
}