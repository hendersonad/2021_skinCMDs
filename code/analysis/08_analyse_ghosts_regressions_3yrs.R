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
    datapath <- "/Users/lsh1510922/Documents/Postdoc/2021_extract/"
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
      readRDS(paste0(datapath, "out/ecz-anxiety_split.rds"))
    df_dep_split <-
      readRDS(paste0(datapath, "out/ecz-depression_split.rds"))
    var_consbeforeindex <-
      haven::read_dta(paste0(datapath, "out/variables-ecz-consultations-yrbeforeindex-3yrs.dta"))
  } else if (exposure == "psoriasis") {
    df_anx_split <-
      readRDS(paste0(datapath, "out/pso-anxiety_split.rds"))
    df_dep_split <-
      readRDS(paste0(datapath, "out/pso-depression_split.rds"))
    var_consbeforeindex <-
      haven::read_dta(paste0(datapath, "out/variables-pso-consultations-yrbeforeindex-3yrs.dta"))
  }
  
  .dib(exposure)
  
  df_anx_split$t <-
    as.numeric(df_anx_split$tstop - df_anx_split$tstart)
  df_dep_split$t <-
    as.numeric(df_dep_split$tstop - df_dep_split$tstart)
  
  df_anx_split$gc90days <-
    factor(df_anx_split$gc90days, levels = c("not active", "active"))
  df_dep_split$gc90days <-
    factor(df_dep_split$gc90days, levels = c("not active", "active"))
  
  df_anx_split$bmi2 <-
    df_anx_split$bmi - mean(df_anx_split$bmi, na.rm = T)
  df_dep_split$bmi2 <-
    df_dep_split$bmi - mean(df_dep_split$bmi, na.rm = T)
  
  for (outcome in YY) {
    #outcome = YY[2]
    .dib(outcome)
    
    if (outcome == "anxiety") {
      df_model <- df_anx_split
    } else if (outcome == "depression") {
      df_model <- df_dep_split
    }
    
    # Bit of variable formatting for output in regression tables --------------
    df_model <- df_model %>% 
      mutate(exposed = case_when(
        exposed == 0 ~ "Unexposed",
        exposed == 1 ~ stringr::str_to_title(paste0(exposure))
      ))
    
    df_model$exposed <- factor(df_model$exposed, levels = c("Unexposed", str_to_title(exposure)))
    var_label(df_model$exposed) <- "Exposure"
    var_label(df_model$carstairs) <- "Carstairs index of deprivation"
    levels(df_model$carstairs)[1] <- "1 (least deprived)"
    levels(df_model$carstairs)[5] <- "5 (most deprived)"
    var_label(df_model$cal_period) <- "Calendar Period"
    var_label(df_model$cci) <- "Charlson's comorbidity index"
    levels(df_model$cci) <- c("Low", "Moderate", "Severe")
    var_label(df_model$agegroup) <- "Age group"
    df_model$agegroup <- relevel(df_model$agegroup, ref = "50-59")
    
    # Mediators
    var_label(df_model$bmi2) <- "BMI (centred)"
    var_label(df_model$alc) <- "Alcohol misuse"
    levels(df_model$alc) <- c("No", "Yes")
    var_label(df_model$smokstatus) <- "Smoking status"
    levels(df_model$smokstatus) <- stringr::str_to_title(levels(df_model$smokstatus))
    var_label(df_model$gc90days) <- "Recent high dose glucocorticoid steroid use (<30days)"
    levels(df_model$gc90days) <- c("No", "Yes")
    df_model$sleep <- factor(df_model$sleep, levels = 0:1, labels = c("No", "Yes"))
    var_label(df_model$sleep) <- "Sleep problems"
    var_label(df_model$severity) <- paste(str_to_title(exposure), "severity", sep = " ")
    levels(df_model$severity) <- str_to_title(levels(df_model$severity))
    if (ABBRVexp == "ecz") {
      df_model$comorbid <- df_model$asthma
      var_label(df_model$comorbid) <- "Asthma"
      levels(df_model$comorbid) <- c("No", "Yes")
    } else{
      df_model$comorbid <- df_model$arthritis
      var_label(df_model$comorbid) <- "Arthritis"
      levels(df_model$comorbid) <- c("No", "Yes")
    }
    
    # add in filtering for consultation <= 1year before index -----------------
    df_model_noghosts <- df_model %>% 
      left_join(var_consbeforeindex, by = "patid") %>% 
      mutate(pre_cons = replace_na(consyrbeforeindex, 0))
    table_noghosts <- df_model_noghosts %>% 
      count(exposed, pre_cons) %>% 
      group_by(exposed) %>% 
      mutate(prop = prop.table(n)*100 %>% signif(digits = 2)) %>% 
      pivot_wider(exposed, names_from = pre_cons, values_from = c(n, prop))

    saveRDS(table_noghosts, file = paste0(datapath, "out/table", ABBRVexp, "_", outcome, "_noghosts.rds"))
    
    df_model_noghosts <- df_model_noghosts %>% 
      filter(pre_cons == 1)

    # only keep valid sets ----------------------------------------------------
    validsets_noghosts <- df_model_noghosts %>% 
      group_by(setid) %>% 
      summarise(mean = mean(exposed == str_to_title(exposure))) %>% 
      filter(mean > 0 & mean < 1)
    validsets <- validsets_noghosts$setid
    
    df_validsets_noghosts <- df_model_noghosts %>% 
      filter(setid %in% validsets)
    
    saveRDS(df_validsets_noghosts, file = paste0(datapath, "out/models_data/df_model", ABBRVexp, "_", outcome, "_noghosts-3yrs.rds"))
    
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
          Surv(t, out) ~ exposed + carstairs + cal_period + comorbid + cci + bmi2 + sleep + alc + smokstatus + gc90days + strata(setid),
          data = df_validsets_noghosts
        ) 
    } else if (ABBRVexp == "pso") {
      mod3 <-
        coxph(
          Surv(t, out) ~ exposed + carstairs + cal_period + comorbid + cci + bmi2 + alc + smokstatus + strata(setid),
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
  
  if (sum(grepl(pattern = "df_model", x = ls())) == 1) {
    rm(df_model)
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