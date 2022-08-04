library(tidyverse)
library(here)
library(flextable)
library(gt)
library(gtsummary)
library(janitor)
library(timetk)
library(skimr)
library(glue)


if (Sys.info()["user"] == "lsh1510922") {
  if (Sys.info()["sysname"] == "Darwin") {
    datapath <- "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
  if (Sys.info()["sysname"] == "Windows") {
    datapath <- "Z:/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
}
XX <- c("psoriasis", "eczema")
#exposure <- "eczema"

for (exposure in XX) {
  ABBRVexp <- str_sub(exposure, 1 , 3)
  
  .dib(exposure)
  if (exposure == "eczema") {
    df_anx_split <-
      readRDS(paste0(datapath, "out/df_modelecz_anxiety.rds"))
    df_dep_split <-
      readRDS(paste0(datapath, "out/df_modelecz_depression.rds"))
  } else if (exposure == "psoriasis") {
    df_anx_split <-
      readRDS(paste0(datapath, "out/df_modelpso_anxiety.rds"))
    df_dep_split <-
      readRDS(paste0(datapath, "out/df_modelpso_depression.rds"))
  }
  
  # psoriasis
  # levels(df_anx$severity) <- c("severe", "mild")
  # levels(df_dep$severity) <- c("severe", "mild")
  table(df_anx_split$severity)
  table(df_dep_split$severity)
  
  ####
  # There seem to be some people with severe records while unexposed controls (before becoming cases). Which seems mental? What is a happening?????
  # Currently just suppressing the severe in unexposed but need to investigate
  # There also seems to be an absolute mess with the factor levels (from build1.R)
  
  # Build flat baseline char. -----------------------------------------------------------------
  build_baseline <- function(df_in = slice(df_anx_split, 1:10000)) {
    df_in$severity[df_in$exposed == "Unexposed"] <- "None" ## MANUAL SUPPRESSION: not ideal but 94 unexposed~severe in psoriasis (at baseline) so not too big an issue
    if (exposure == "eczema") {
      df_in <- df_in %>%
        mutate(comorbid = asthma)
      #severity = factor(severity, levels = c(0, "mild", "severe"), labels=c("Mild", "Moderate", "Severe")))
    } else if (exposure == "psoriasis") {
      df_in <- df_in %>%
        mutate(comorbid = arthritis)
      #severity = factor(severity, levels = c(0, "mild"), labels=c("Mild", "Severe")))
    }
    df_out <- df_in %>%
      arrange(setid, patid, tstart) %>%
      group_by(setid, patid) %>%
      slice(1)
    
    tab1 <- df_out %>%
      ungroup() %>%
      mutate(fup = (enddate - indexdate) / 365.25) %>%
      dplyr::select(setid,patid,exposed,gender,age,agegroup,country,ruc,fup,cal_period,eth_edited,carstairs,ruc,bmi,bmi_cat,alc,smokstatus,sleep,sleep_all,comorbid,cci,severity,out) ## need to add DEATH here once it is working properly but seems to have been corrupted by the stsplit (death=1 being copied over multiple lines which is non-sensical)
    table(tab1$exposed, tab1$severity, useNA = "always")
    
    ## make table with nice gtsummary
    table1 <- tab1 %>%
      ungroup() %>%
      dplyr::select(-patid,-setid,-out) %>%
      tbl_summary(
        by = exposed,
        statistic = list(
          all_continuous() ~ "{p50} ({p25}-{p75})",
          all_categorical() ~ "{n} ({p}%)"
        ),
        digits = all_continuous() ~ 1,
        label = list(
          gender = "Sex",
          fup = "Follow-up time (years)",
          cal_period = "Calendar period",
          age = "Age",
          agegroup = "Age (categorised)",
          eth_edited = "Ethnicity",
          bmi = "BMI",
          bmi_cat = "Obesity (categorised)",
          alc = "Harmful alcohol use",
          comorbid = ifelse(
            exposure == "eczema",
            "Asthma diagnosis",
            "Arthritis diagnosis"
          ),
          sleep = "Sleep problems",
          sleep_all = "Sleep problems (incl. benzos)",
          cci = "Charlson's comorbidity index",
          smokstatus = "Smoking status",
          carstairs = "Carstairs deprivation quintile",
          ruc = "Rural/Urban",
          country = "Country",
          severity = paste0("Severe ", exposure)
        )
      ) %>%
      add_overall() %>%
      bold_labels() %>%
      modify_footnote(all_stat_cols() ~ "Median (IQR) or Frequency (%)")
    
    table1
  }
  .dib("Building Table 1")
  table1_anx <- build_baseline(df_anx_split)
  table1_dep <- build_baseline(df_dep_split)
  
  table1 <- tbl_merge(
    tbls = list(table1_anx, table1_dep),
    tab_spanner = c("**Anxiety cohort**", "**Depression cohort**")
  )
  table1
  table1 %>%
    as_gt() %>%
    gt::gtsave(
      filename = paste0("tab1_", ABBRVexp, ".html"),
      path = here::here("out/tables")
    )
  
  # Build table 2 - any exposure during follow up ---------------------------
  build_followup <- function(df_in) {
    df_out <- df_in %>%
      ungroup() %>% 
      mutate(smok_recategorised = factor(smokstatus, levels = c("Non-Smoker", "Ex-Smoker", "Current Or Ex-Smoker", "Current Smoker"))) %>% 
      mutate(severity_recategorised = factor(severity, levels = c("None", "Mild",  "Moderate", "Severe"))) %>% 
      mutate(smokNum = as.numeric(smok_recategorised)) %>% 
      mutate(severityNum = as.numeric(severity_recategorised)) %>% 
      group_by(setid, patid) %>%
      mutate(everSev = factor(max(severityNum, na.rm = T), levels = 1:4, labels = c("None", "Mild",  "Moderate", "Severe")),
        everOut = ifelse(any(out == 1), 1, 0),
        everAlc = ifelse(any(alc == "Yes"), 1, 0),
        everSmok = factor(max(smokNum, na.rm = T), levels = 1:4, labels = c("Non-Smoker", "Ex-Smoker", "Current Or Ex-Smoker", "Current Smoker")),
        everSleep = ifelse(any(sleep == "Yes"), 1, 0),
        everSteroid = ifelse(any(gc90days == "Yes"), 1, 0),
        everComorbid = ifelse(any(comorbid == "Yes"), 1, 0)
      ) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(fup = (enddate - indexdate) / 365.25) %>%
      dplyr::select(
        setid,
        patid,
        exposed,
        everAlc,
        everSmok,
        everSleep,
        everComorbid,
        everSev,
        everSteroid,
        everOut
      )
    
    tab2 <- df_out 
    
    ## investigating the unexposed people with severe psoriasis Rx codes
    table(tab2$exposed, tab2$everSev)
    
    ## make table with nice gtsummary
    if (exposure == "psoriasis") {
      tab2 <- dplyr::select(tab2,-everSteroid)
    }
    table2 <- tab2 %>%
      dplyr::select(-patid,-setid,-everOut) %>%
      tbl_summary(
        by = exposed,
        statistic = list(
          all_continuous() ~ "{p50} ({p25}-{p75})",
          all_categorical() ~ "{n} ({p}%)"
        ),
        digits = all_continuous() ~ 1,
        label = list(
          gender = "Sex",
          fup = "Follow-up time (years)",
          cci = "Charlson's comorbidity index",
          age = "Age",
          agegroup = "Age (categorised)",
          eth_edited = "Ethnicity",
          bmi = "BMI",
          bmi_cat = "Obesity (categorised)",
          #everOut = outcome_cohort,
          everAlc = "Harmful alcohol use",
          everComorbid = ifelse(
            exposure == "eczema",
            "Asthma diagnosis",
            "Arthritis diagnosis"
          ),
          everSleep = "Sleep problems",
          everSmok = "Smoking status",
          carstairs = "Carstairs deprivation quintile",
          ruc = "Rural/Urban",
          everSteroid = "High dose oral steroid prescription",
          everSev = paste0("Severe ", exposure)
        )
      ) %>%
      add_overall() %>%
      bold_labels() %>%
      modify_footnote(all_stat_cols() ~ "Median (IQR) or Frequency (%)")
    
    table2
  }
  .dib("Building Table 2")
  
  table2_anx <- build_followup(df_anx_split)
  table2_dep <- build_followup(df_dep_split)
  
  table2 <- tbl_merge(
    tbls = list(table2_anx, table2_dep),
    tab_spanner = c("**Anxiety cohort**", "**Depression cohort**")
  )
  table2
  table2 %>%
    as_gt() %>%
    gt::gtsave(
      filename = paste0("tab2_", ABBRVexp, "_fup.html"),
      path = here::here("out/tables")
    )
  
  
  # Build table 3 - person years ---------------------------------------------------
  ## will need this little function later to summ follow up by factor vars
  summ_pyars <- function(V = "alc", df = tab1) {
    df %>%
      group_by(exposed, get(V)) %>%
      summarise(pyar = sum(t)) %>%
      rename(var = `get(V)`) %>%
      mutate(name = V,
             group_pyar = sum(pyar)) %>%
      ungroup()
  }
  
  ## and this is the main function
  build_pyears <- function(df_in) {
    df_in$severity[df_in$exposed=="Unexposed"] <- "None"
    tab1 <- df_in %>%
      as_tibble() %>%
      ungroup() %>%
      arrange(setid, patid)
    
    x <- tab1 %>%
      select_if(is.factor) %>%
      dplyr::select(-arthritis,-asthma,-age_cat, -exposed) %>%
      names() %>% 
      as.list()
    
    pyars_table <-
      map(x, summ_pyars, df = tab1) %>% ## to use the little function defined above
      do.call(rbind, .)
    
    if(exposure == "psoriasis"){
      pyars_table <- pyars_table %>% 
        filter(name != "gc90days")
    }
    
    gt_pyars_table <- pyars_table %>%
      group_by(exposed) %>%
      mutate(pc_pyar = (pyar / group_pyar) * 100) %>%
      pivot_wider(
        id_cols = c(name, var),
        names_from = c(exposed),
        values_from = c(pyar, pc_pyar)
      ) %>%
      clean_names() %>%
      mutate_at("var", ~ stringr::str_to_title(.)) %>%
      mutate(
        new_lab = case_when(
          name == "exposed" ~ "Exposed",
          name == "gender" ~ "Sex",
          name == "agegroup" ~ "Age (categorised)",
          name == "eth_edited" ~ "Ethnicity",
          name == "bmi_cat" ~ "Obesity (categorised)",
          name == "alc" ~ "Harmful alcohol use",
          name == "comorbid" ~ ifelse(
            exposure == "eczema",
            "Asthma diagnosis",
            "Arthritis diagnosis"
          ),
          name == "cci" ~ "Charlson's comorbidity index",
          name == "sleep" ~ "Sleep problems",
          name == "smokstatus" ~ "Smoking status",
          name == "carstairs" ~ "Carstairs deprivation quintile",
          name == "cal_period" ~ "Calendar period",
          name == "gc90day" ~ "Oral steroid prescription",
          name == "ruc" ~ "Rural/Urban",
          name == "country" ~ "Country",
          name == "severity" ~ "Severity",
        )
      ) %>%
      dplyr::select(-name) %>%
      arrange(new_lab)
    
    gt_pyars_table %>%
      rename(pyar_exposure = paste0("pyar_", exposure),
             pc_pyar_exposure = paste0("pc_pyar_", exposure)) %>% 
      gt(rowname_col = "var",
         groupname_col = "new_lab") %>%
      tab_stubhead(label = "Characteristic") %>%
      fmt_number(
        columns = where(is.numeric),
        decimals = 1,
        use_seps = T
      ) %>%
      cols_merge(columns = contains("unexposed"),
                 pattern = "{1} ({2}%)") %>%
      cols_merge(columns = contains("exposure"),
                 pattern = "{1} ({2}%)") %>%
      tab_style(style = cell_text(weight = "bold"),
                locations = cells_row_groups(groups = everything())) %>%
      cols_label(
        pyar_unexposed = md("**Matched controls**"),
        pyar_exposure = md(paste0("**With ", exposure, "**"))
      )
  }
  .dib("Building Table 3")
  
  table3_anx <- build_pyears(df_in = df_anx_split)
  table3_dep <- build_pyears(df_in = df_dep_split)

  table3_anx %>%
    gt::gtsave(
      filename = paste0("tab3_", ABBRVexp, "_pyars_anxiety.html"),
      path = here::here("out/tables")
    )
  
  table3_dep %>%
    gt::gtsave(
      filename = paste0("tab3_", ABBRVexp, "_pyars_depression.html"),
      path = here::here("out/tables")
    )
  
}