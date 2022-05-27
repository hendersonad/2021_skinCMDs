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

if (Sys.info()["user"] == "lsh1510922") {
  if (Sys.info()["sysname"] == "Darwin") {
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
  #for(outcome in YY){
  df_model_anx_imputed <- readRDS(paste0(datapath, "out/df_model", ABBRVexp, "_anxiety.rds"))
  df_model_dep_imputed <- readRDS(paste0(datapath, "out/df_model", ABBRVexp, "_depression.rds"))
  
  ## combine these datasets
  ## filter to one row
  fn_onerow <- function(df, exp = ABBRVexp){
    temp <- df %>% 
      group_by(setid, patid) %>% 
      slice(1) %>% 
      ungroup() 
    
    if(exp == "ecz"){
      temp %>% 
        select(setid, patid, dob, indexdate, gender, exposed, eth_edited, carstairs, cal_period, comorbid, cci, bmi_cat, bmi, alc, smokstatus, smokstatus_original, sleep, gc90days) %>% 
        rename(smoking_imputed = smokstatus, obesity_categorised = bmi_cat, ethnicity = eth_edited) %>% 
        mutate(age = as.numeric(indexdate - dob)/365.25) %>% 
        select(-dob, -indexdate)
    } else if(exp == "pso"){
      temp %>% 
        select(setid, patid, dob, indexdate, gender, exposed, eth_edited, carstairs, cal_period, comorbid, cci, bmi_cat, bmi, alc, smokstatus, smokstatus_original) %>% 
        rename(smoking_imputed = smokstatus, obesity_categorised = bmi_cat, ethnicity = eth_edited) %>% 
        mutate(age = as.numeric(indexdate - dob)/365.25) %>% 
        select(-dob, -indexdate)
    }
  }
  df_model_anx_1row <- fn_onerow(df_model_anx_imputed)
  df_model_dep_1row <- fn_onerow(df_model_dep_imputed)
  
  rm(df_model_anx_imputed, df_model_dep_imputed)
  
  df_model_merge <- df_model_anx_1row %>% 
    full_join(df_model_dep_1row)
  
  df_missingsum <- df_model_merge %>% 
    ungroup() %>% 
    group_by(exposed) %>% 
    summarise_all(~sum(is.na(.)))
  
  df_missingN <- df_model_merge %>% 
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
  
  
  df_missingplot
  
  pdf(paste0(here::here("out/supplementary"), "/missing_", ABBRVexp, ".pdf"), width = 8, height = 6)
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
  
  ## describe characteristics by missing status 
  missing_data <- NULL
  for(var in c("carstairs", "smokstatus_original", "bmi","obesity_categorised", "smoking_imputed", "ethnicity")){
    patids_missing <- df_model_merge %>% 
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
    distinct(setid, patid, .keep_all = TRUE) %>% 
    mutate(missing = 1+as.numeric(var_missing == "ethnicity")) ## missing = 1 or 2 if only missing ethnicity
  
  df_model_miss_1obs <- df_model_merge %>% 
    left_join(missing_data, by = c("setid", "patid")) %>% 
    mutate_at("missing", ~ifelse(is.na(.), 0, .))
  
  var_label(df_model_miss_1obs$smoking_imputed) <- "Smoking status (imputed)"
  df_model_miss_1obs$gender <- factor(as.character(df_model_miss_1obs$gender), levels = c("Male", "Female"))
  var_label(df_model_miss_1obs$gender) <- "Sex"
  var_label(df_model_miss_1obs$age) <- "Age at index"
  var_label(df_model_miss_1obs$ethnicity) <- "Ethnicity"
  var_label(df_model_miss_1obs$carstairs) <- "Carstairs index of deprivation"
  var_label(df_model_miss_1obs$smokstatus_original) <- "Original smoking data"
  var_label(df_model_miss_1obs$bmi) <- "Body Mass Index (BMI)"
  var_label(df_model_miss_1obs$obesity_categorised) <- "Obesity status"
  var_label(df_model_miss_1obs$cal_period) <- "Calendar Period"
  var_label(df_model_miss_1obs$cci) <- "Charlson's comorbidity index"
  if (ABBRVexp == "ecz") {
    var_label(df_model_miss_1obs$comorbid) <- "Asthma"
  } else{
    var_label(df_model_miss_1obs$comorbid) <- "Arthritis"
  }
  
  
  
  tab_missing <- df_model_miss_1obs %>% 
    ungroup() %>% 
    mutate(missing = factor(missing, levels = 0:2, labels = c("Complete", "Missing other data", "Missing ethnicity data"))) %>% 
    select(-patid, -setid, -var_missing) %>%
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
