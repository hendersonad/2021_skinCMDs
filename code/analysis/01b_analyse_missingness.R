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
library(naniar)
library(visdat)

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
exposure = XX[2]
outcome = YY[2]

for(exposure in XX){
  ABBRVexp <- substr(exposure, 1, 3)
  .dib(exposure)
  for(outcome in YY){
    df_model <- readRDS(paste0(datapath, "out/df_model", ABBRVexp, "_", outcome, ".rds"))
    df_model_1row <- df_model %>% 
      group_by(setid, patid) %>% 
      slice(1) %>% 
      ungroup()
      
    ## visualise missing data
    if (ABBRVexp == "ecz") {
      df_model_select <- df_model_1row %>% 
        select(exposed, carstairs, cal_period, comorbid, cci, bmi, obese_cat, alc, smokstatus, smokstatus_nomiss, sleep, gc90days) %>% 
        rename(smoking_imputed = smokstatus_nomiss, obesity_categorised = obese_cat)
    } else if (ABBRVexp == "pso") {
      df_model_select <- df_model_1row %>% 
        select(exposed, carstairs, cal_period, comorbid, cci, bmi, obese_cat, alc, smokstatus, smokstatus_nomiss) %>% 
        rename(smoking_imputed = smokstatus_nomiss, obesity_categorised = obese_cat)
    }
    
    df_missingsum <- df_model_select %>% 
      ungroup() %>% 
      group_by(exposed) %>% 
      summarise_all(~sum(is.na(.)))
    
    df_missingN <- df_model_select %>% 
      group_by(exposed) %>% 
      summarise(n = n())
    
    df_missingplot <- df_missingsum %>% 
      pivot_longer(!exposed) %>% 
      arrange(value) %>% 
      filter(name == "obesity_categorised" | value>0)  %>% 
      left_join(df_missingN, by = "exposed") %>% 
      mutate(pc = (value/n)*100)
    
    fct_levels <- df_missingplot %>% 
      ungroup() %>% 
      group_by(name) %>% 
      summarise(val = sum(value)) %>% 
      arrange(val) %>% 
      pull(name)
    fct_labels <- stringr::str_replace(fct_levels, "_", " ")
    fct_labels <- stringr::str_to_title(fct_labels)
    
    df_missingplot$plotname <- factor(df_missingplot$name, levels = fct_levels, labels = fct_labels)
    
    pdf(paste0(here::here("out/supplementary/"), "missing_", ABBRVexp, "_", substr(outcome,1,3), ".pdf"), width = 8, height = 6)
      p1 <- ggplot(df_missingplot, aes(x = plotname, y = pc, ymax = pc, ymin = 0, group = exposed, col = exposed)) +
        geom_linerange() +
        geom_point() +
        coord_flip() +
        facet_wrap(~exposed) +
        labs(x = "Variable", y = "% missing") +
        theme_ali() +
        theme(strip.background = element_blank(),
              legend.position = "none")
      print(p1)
    dev.off()
  }
  
  ## describe characteristics by missing status 
  if (ABBRVexp == "ecz") {
    df_model_select <- df_model_1row %>% 
      select(setid, patid, dob, indexdate, gender, exposed, carstairs, cal_period, comorbid, cci, bmi, obese_cat, alc, smokstatus, smokstatus_nomiss, sleep, gc90days) %>% 
      rename(smoking_imputed = smokstatus_nomiss, obesity_categorised = obese_cat) %>% 
      mutate(age = as.numeric(indexdate - dob)/365.25) %>% 
      select(-dob, -indexdate)
  } else if (ABBRVexp == "pso") {
    df_model_select <- df_model_1row %>% 
      select(setid, patid, dob, indexdate, gender, exposed, carstairs, cal_period, comorbid, cci, bmi, obese_cat, alc, smokstatus, smokstatus_nomiss) %>% 
      rename(smoking_imputed = smokstatus_nomiss, obesity_categorised = obese_cat) %>% 
      mutate(age = as.numeric(indexdate - dob)/365.25) %>% 
      select(-dob, -indexdate)
  }
  
  missing_data <- NULL
  for(var in c("carstairs", "smokstatus", "bmi")){
    patids_missing <- df_model_select %>% 
      select(setid, patid, all_of(var)) %>% 
      filter(is.na(get(var))) %>% 
      distinct(setid, patid)
    missing_data <- missing_data %>% 
      bind_rows(
        bind_cols(
          patids_missing,
          var_missing = var
          )
      )
  }
  missing_data <- missing_data %>% 
    distinct(setid, patid) %>% 
    mutate(missing = 1)
  
  missing_data_impute <- NULL
  for(var in c("obesity_categorised", "smoking_imputed")){
    patids_missing <- df_model_select %>% 
      select(setid, patid, all_of(var)) %>% 
      filter(is.na(get(var))) %>% 
      distinct(setid, patid)
    missing_data_impute <- missing_data_impute %>% 
      bind_rows(
        bind_cols(
          patids_missing,
          var_missing = var
          )
      )
  }
  missing_data_impute <- missing_data_impute %>% 
    distinct(setid, patid) %>% 
    mutate(missing = 1)
  
  missing_either <- missing_data %>% 
    left_join(missing_data_impute, by = c("setid", "patid")) %>% 
    rename(missing_original = missing.x, missing_impute = missing.y) %>% 
    mutate_at(c("missing_original", "missing_impute"), ~ifelse(is.na(.), 0, .)) %>% 
    mutate(missing = missing_original + missing_impute) %>% 
    select(-missing_original, -missing_impute)
  
  df_model_miss_1obs <- df_model_select %>% 
    left_join(missing_either, by = c("setid", "patid")) %>% 
    mutate_at("missing", ~ifelse(is.na(.), 0, .))
  
  # df_model_miss_1obs <- df_model_miss %>% 
  #   group_by(setid,patid) %>% 
  #   slice(1)
  
  var_label(df_model_miss_1obs$smoking_imputed) <- "Smoking status (imputed)"
  df_model_miss_1obs$gender <- factor(as.character(df_model_miss_1obs$gender), levels = c("Male", "Female"))
  var_label(df_model_miss_1obs$gender) <- "Sex"
  var_label(df_model_miss_1obs$age) <- "Age at index"
  
  tab_missing <- df_model_miss_1obs %>% 
    ungroup() %>% 
    mutate(missing = factor(missing, levels = 0:2, labels = c("Complete", "Missing data", "Missing data after imputation"))) %>% 
    select(-patid, -setid) %>%
    select(age, gender, everything()) %>% 
    tbl_strata(strata = exposed,
               .tbl_fun = 
                 ~ .x %>% 
                 tbl_summary(by = missing,
                             statistic = list(all_continuous() ~ "{p50} ({p25}-{p75})",
                                              all_categorical() ~ "{n} ({p}%)"),
                             digits = all_continuous() ~ 1
                             ) %>%
                 bold_labels() %>%
                 modify_footnote(
                   all_stat_cols() ~ "Median (IQR) or Frequency (%)"
                 )
    )
  tab_missing
  tab_missing %>%
    as_gt() %>%
    gt::gtsave(filename = paste0("tab4_",ABBRVexp,"_missingness.html"), path = here::here("out/tables"))
  
}
