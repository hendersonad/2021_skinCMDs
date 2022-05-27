library(here)
library(haven)
library(tidyverse)
library(survival)
library(survminer)
library(ggplot2)

if(Sys.info()["user"]=="lsh1510922"){
  if(Sys.info()["sysname"]=="Darwin"){
    datapath <- "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
  if(Sys.info()["sysname"]=="Windows"){
    datapath <- "Z:/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
}

XX <- c("psoriasis", "eczema")

  exposure <- XX[1]
  
carstairs_avail <- function(X = XX[1]) {
  ABBRVexp <- str_sub(X,1 ,3)
  .dib(X) 

  # import data -------------------------------------------------------------
  patient <- haven::read_dta(paste0(datapath, "in/Patient_extract_", ABBRVexp, "_extract3_1.dta"))
  cohort <- haven::read_dta(paste0(datapath, "out/getmatchedcohort-", X, "-main-mhealth.dta"))
  
  ## linked data
  var_catstairs_pt <- readr::read_delim(file = paste0(datapath, "linked/patient_carstairs_20_051.txt"), delim = "\t")
  var_catstairs_prac <- readr::read_delim(file = paste0(datapath, "linked/practice_carstairs_20_051.txt"), delim = "\t")
  var_ruc <- readr::read_delim(file = paste0(datapath, "linked/practice_urban_rural_20_051.txt"), delim = "\t")
  
  # get one variable for RUC
  var_ruc <- var_ruc %>%
    group_by(pracid) %>%
    mutate(ruc = max(ni2015_urban_rural, e2011_urban_rural, w2011_urban_rural, s2016_urban_rural, na.rm = T)) %>%
    ungroup()
  
  pt_level <- cohort %>% 
    left_join(var_catstairs_pt, by = "patid") %>% 
    mutate(carstairs_pt = as.numeric(!is.na(carstairs2011_5)))  
  pt_prac_level <- pt_level %>% 
    dplyr::select(-carstairs2011_5) %>%
    mutate(pracid = patid %% 100) %>% 
    left_join(var_catstairs_prac, by = "pracid") %>% 
    mutate(carstairs_prac = as.numeric(!is.na(carstairs2011_5)))  %>% 
    dplyr::select(-carstairs2011_5) 
  
  twoXtwo <- function(df, exp, out){
    df1 <- df %>% 
      ungroup() %>% 
      dplyr::select(exp = {{ exp }}, out = {{ out }})
    tab <- table(df1$exp, df1$out, useNA = "always")
    tab_p <- prop.table(tab,1)
    tibble(
      exposure = exp,
      val = rownames(tab),
      No = tab[, 1],
      No_pc = tab_p[, 1] * 100,
      Yes = tab[, 2],
      Yes_pc = tab_p[, 2] * 100,
      Miss = tab[, 3]
    )
  }
  
  tab1 <- twoXtwo(pt_prac_level, "exposed", "carstairs_pt") %>% mutate(level = "patient")
  tab2 <- twoXtwo(pt_prac_level, "exposed", "carstairs_prac") %>% mutate(level = "practice")
  carstairs_avail <- tab1 %>% 
    bind_rows(tab2) %>%
    mutate(exposure = X) 
  carstairs_avail
}

pso_carstairs <- carstairs_avail("psoriasis")
ecz_carstairs <- carstairs_avail("eczema")

carstairs_data <- pso_carstairs %>% 
  bind_rows(ecz_carstairs) %>% 
  drop_na()

gt_carstairs <- carstairs_data %>%
  arrange(level, exposure) %>% 
  mutate(val = ifelse(val == 0, paste0("Unexposed (", exposure, " cohort)"), str_to_title(exposure))) %>% 
  dplyr::select(val, level, everything(), -level, -Miss, -exposure) %>%
  gt::gt() %>%
  tab_row_group(label = "Patient level",
                rows = 1:4) %>%
  tab_row_group(label = "Practice level",
                rows = 5:8) %>%
  gt::tab_header(title = "Carstairs data availability in matched cohorts") %>%
  gt::fmt_number(columns = c(3, 5), decimals = 1) %>%
  gt::fmt_number(columns = c(2, 4), decimals = 0) %>%
  gt::data_color(
    columns = c(Yes_pc),
    colors = scales::col_numeric(
      palette = paletteer::paletteer_c(palette = "viridis::inferno",
                                       n = 100) %>% as.character(),
      domain = c(
        min(carstairs_data$Yes_pc, na.rm = T),
        max(carstairs_data$Yes_pc, na.rm = T)
      )
    )
  ) %>%
  cols_label(val = "",
             No = "Missing", 
             Yes = "Recorded",
             No_pc = "%",
             Yes_pc = "%") 
gt_carstairs %>% 
  gt::gtsave(
    filename =  paste0("carstairs_data_availability.html"),
    path = here::here("out/supplementary")
  )
