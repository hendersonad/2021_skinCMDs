library(tidyverse)
library(data.table)
library(here)
library(magrittr)
library(survival)
library(gt)

if (Sys.info()["user"] == "lsh1510922") {
  if (Sys.info()["sysname"] == "Darwin") {
    datapath <- "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
  if (Sys.info()["sysname"] == "Windows") {
    datapath <- "Z:/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
}
# Estimates rate differences using HR outputs -----------------------------
XX <- c("psoriasis", "eczema")
exposure <- XX[1]
get_rate_difference <- function(exposure){
  
  ABBRVexp <- str_sub(exposure, 1 , 3)
  # load models -------------------------------------------------------------
  mod2_anx <-
    readRDS(paste0(
      datapath,
      "out/models_data/",
      ABBRVexp,
      "_anxiety_mod2_modeldata.rds"
    ))
  
  # load models -------------------------------------------------------------
  mod2_dep <-
    readRDS(paste0(
      datapath,
      "out/models_data/",
      ABBRVexp,
      "_depression_mod2_modeldata.rds"
    ))
  
  # load data ---------------------------------------------------------------
  df_model_anx <-
    readRDS(paste0(
      datapath,
      "out/df_model",
      ABBRVexp,
      "_anxiety.rds"
    ))
  df_model_dep <-
    readRDS(paste0(
      datapath,
      "out/df_model",
      ABBRVexp,
      "_depression.rds"
    ))
  
  rd_calculation <- function(outcome){
    df_model <- get(paste0("df_model_", substr(outcome,1,3)))
    model_output <- get(paste0("mod2_", substr(outcome,1,3)))
    # get rid of patids with missing data
    df_model_noNA <- df_model %>% 
      dplyr::select(setid, patid, t, out, exposed, carstairs, cal_period, comorbid, cci) %>% 
      drop_na()
    
    ## get p-yars at risk
    ## get n events
    ## calculate incidence rate per 100,000 years
    exposed_data <- df_model_noNA %>% 
      filter(exposed==str_to_title(exposure))  
    exposed_summary <- exposed_data %>% 
      summarise(pyears = sum(t/365.25), events = sum(out)) %>% 
      mutate(pyears_000 = pyears/1e5,
             rate_exposed = events/pyears_000)
    exposed_rate <- exposed_summary$rate_exposed
    
    ## take confounder adjusted HR and invert it
    model_estimates <- model_output %>% 
      broom::tidy(conf.int = T) %>% 
      filter(str_detect("exposed", string = term)) %>% 
      mutate(HR = exp(estimate))
    invertedHR <- 1/model_estimates$HR
    
    ## estimate unexposed rate as exposed_rate * 1/HR
    unexposed_rate <- exposed_rate * invertedHR
    
    ## calc rate difference 
    rate_diff_simple <- exposed_rate - unexposed_rate
    
    ## bootstrap confidence interval 
    patids <- unique(exposed_data$patid) 
    n <- patids %>% length()
    b <- 5000
    
    exp_btsp <- vector(length = b)
    unexp_btsp <- vector(length = b)
    rd_btsp <- vector(length = b)
    set.seed(153342)
    for(btsp in 1:b){
      if(btsp%%1000 == 0){print(c("Run ", btsp))}
      hr_draw <- rnorm(1, mean = model_estimates$estimate,
                       sd = model_estimates$std.error) 
      index <- sample(patids, size = n, replace = TRUE)
      exposed_btsp <- exposed_data %>% 
        filter(patid %in% index)
      exposed_btsp_summ <- exposed_btsp %>% 
        summarise(pyears = sum(t/365.25), events = sum(out)) %>% 
        mutate(pyears_000 = pyears/1e5,
               rate_exposed = events/pyears_000)
      exp_btsp[btsp] <- exposed_btsp_summ$rate_exposed
      unexp_btsp[btsp] <- exp_btsp[btsp] * 1/exp(hr_draw)
      rd_btsp[btsp] <- exp_btsp[btsp] - unexp_btsp[btsp]
    }
    
    sigD <- 1
    rate_exposed_CI <- paste(round(quantile(exp_btsp, prob = c(0.025, 0.975)),sigD), collapse = " - ")
    rate_unexposed_CI <- paste(round(quantile(unexp_btsp, prob = c(0.025, 0.975)),sigD), collapse = " - ")
    rate_diff_CI <- paste(round(quantile(rd_btsp, prob = c(0.025, 0.975)),sigD), collapse = " - ")
    hr_CI <- paste(round(exp(c(model_estimates$conf.low, model_estimates$conf.high)), 2), collapse = " - ")
    
    data.frame(
      exp = exposure, 
      out = outcome,
      pyears = exposed_summary$pyears_000,
      events = exposed_summary$events, 
      rate_exposed = exposed_summary$rate_exposed,
      rate_exposed_CI = rate_exposed_CI,
      hr = model_estimates$HR,
      hr_CI = hr_CI,
      rate_unexposed = unexposed_rate,
      rate_unexposed_CI = rate_unexposed_CI,
      rate_diff = rate_diff_simple,
      rate_diff_CI = rate_diff_CI
    )
  }
  rd_anx <- rd_calculation("anxiety")
  rd_dep <- rd_calculation("depression")
  
  bind_rows(rd_anx, rd_dep)
}

pso_ratediff <- get_rate_difference("psoriasis")
ecz_ratediff <- get_rate_difference("eczema")

rate_differences <- pso_ratediff %>% 
  bind_rows(ecz_ratediff)

gt_rates <- rate_differences %>%
  dplyr::select(-exp) %>% 
  mutate(out = str_to_title(out)) %>% 
  gt::gt() %>%
  tab_row_group(label = "Eczema",
                rows = 3:4) %>%
  tab_row_group(label = "Psoriasis",
                rows = 1:2) %>%
  gt::fmt_number(columns = c(2,4,6,8,10), decimals = 1) %>%
  gt::fmt_number(columns = 3, decimals = 0) %>%
  cols_label(out = "Outcome",
             pyears = "Person-years (100,000)",
             events = "Events",
             rate_exposed = "Rate (exposed group)",
             hr = "Hazard ratio",
             rate_unexposed = "Rate (unexposed group)",
             rate_diff = "Rate difference"
  ) 

rate_differences %>% write_csv(here::here("out/supplementary/rate_difference.csv"))
gt_rates %>% gt::gtsave(filename = "rate_difference.html", path = here::here("out/supplementary/"))

gt_rates_publication <- rate_differences %>%
  dplyr::select(-exp, -pyears, -events, -hr, -hr_CI) %>% 
  mutate(out = str_to_title(out)) %>% 
  gt::gt() %>%
  tab_row_group(label = "Eczema",
                rows = 3:4) %>%
  tab_row_group(label = "Psoriasis",
                rows = 1:2) %>%
  gt::fmt_number(columns = c(2,4,6), decimals = 1) %>%
  cols_merge(columns = c(2,3), pattern = "{1} ({2})") %>% 
  cols_merge(columns = c(4,5), pattern = "{1} ({2})") %>% 
  cols_merge(columns = c(6,7), pattern = "{1} ({2})") %>% 
  cols_label(out = "Outcome",
             rate_exposed = "Rate (exposed group)",
             rate_unexposed = "Rate (unexposed group)",
             rate_diff = "Rate difference"
  ) 
gt_rates_publication %>% gt::gtsave(filename = "rate_difference_pub.html", path = here::here("out/supplementary/"))

