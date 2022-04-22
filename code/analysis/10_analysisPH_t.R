library(tidyverse)
library(here)
library(survival)
library(lmtest)
library(broom)
library(biostat3)
library(multcomp)
library(rms)

if (Sys.info()["user"] == "lsh1510922") {
  if (Sys.info()["sysname"] == "Darwin") {
    #datapath <- "/Users/lsh1510922/Documents/Postdoc/2021_extract/"
    datapath <- "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
  if (Sys.info()["sysname"] == "Windows") {
    datapath <- "Z:/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
}

source(here("code/programs/schonfeld_plot.R")) ## adapted plot code for Scales Schoenfeld residuals
dir.create(file.path(here("out")), showWarnings = FALSE)
dir.create(file.path(here("out", "PHchecks")), showWarnings = FALSE)

YY <- c("depression", "anxiety")
XX <- c("psoriasis", "eczema")

#exposure <- XX[1]
#outcome <- YY[1]

exposure <- XX[1]
outcome <- YY[2]

for(exposure in XX) {
  ABBRVexp <- substr(exposure, 1, 3)
  
  for(outcome in YY) {
    # load data ---------------------------------------------------------------
    df_model <-
      readRDS(paste0(datapath, "out/models_data/df_model", ABBRVexp, "_", outcome,".rds"))
    # 
    # set.seed(30)
    # sample_ids <- sample(unique(df_model$setid), size = 10000, replace = F)
    # df_sample <- df_model %>% 
    #   filter(setid %in% sample_ids)
    # 
    # 
    # fit model with time interaction -----------------------------------------
    fit1_int <-
      coxph(
        Surv(t, out) ~ exposed * t + strata(setid),
        data = df_model
      ) 
    fit2_int <-
      coxph(
        Surv(t, out) ~ exposed * t + carstairs + cal_period + comorbid + cci + strata(setid),
        data = df_model
      ) 
    fit3_int <-
      coxph(
        Surv(t, out) ~ exposed * t + carstairs + cal_period + comorbid + cci + bmi2 + alc + smokstatus + strata(setid),
        data = df_model
      ) 
    
    time_series <- seq(0, 365*5, 365.25/12)
    interactions <- c(
      paste0("exposed", str_to_title(exposure), " + ", time_series, "*exposed", str_to_title(exposure), ":t = 0")
    )
    
    lincom_out_min <- lincom(fit1_int, 
                         interactions, eform = T, level = 0.95)
    lincom_out_con <- lincom(fit2_int, 
                         interactions, eform = T, level = 0.95)
    lincom_out_med <- lincom(fit3_int, 
                         interactions, eform = T, level = 0.95)
    
    
    lincom_out_min <- as_tibble(lincom_out_min, rownames = "lincom") %>% mutate(model = "Minimal")
    lincom_out_con <- as_tibble(lincom_out_con, rownames = "lincom") %>% mutate(model = "Confounder")
    lincom_out_med <- as_tibble(lincom_out_med, rownames = "lincom") %>% mutate(model = "Mediator")
    tibble_int <- lincom_out_min %>%
      bind_rows(lincom_out_con) %>% 
      bind_rows(lincom_out_med) %>% 
      janitor::clean_names() %>%
      dplyr::select(model, 
                    lincom,
                    estimate,
                    conf.low = x2_5_percent,
                    conf.high = x97_5_percent) %>%
      mutate_at(c("estimate", "conf.low", "conf.high"), ~ unlist(.)) %>% 
      mutate(y = outcome, x = exposure, z = "t")
    
    tibble_plot <- tibble_int %>% 
      mutate(t = rep(time_series, 3)) %>% 
      mutate(yr = t/365.25)
    
    plot <- ggplot(tibble_plot, aes(x = yr, y = estimate, ymin = conf.low, ymax = conf.high, group = model, fill = model, col = model)) + 
      geom_line() + 
      geom_ribbon(lty = 0, alpha = 0.2) +
      geom_hline(yintercept = 1, lty = 4, col = 1) +
      ylim(c(0,2)) +
      labs(x = "Year", 
           y = expression(hat(beta) ~ "exposure:t"), 
           title = paste0(str_to_title(exposure), " ~ ", str_to_title(outcome))) +
      guides(fill = guide_legend(title = "Model"), col = guide_legend(title = "Model")) +
      theme_ali()
    
    assign(paste0("plot_", ABBRVexp, "_", substr(outcome, 1, 3)), plot)
  }
}

pdf(paste0(here("out//supplementary//"), "int_t_beta_ecz_anx.pdf"), 8, 8)
  plot_ecz_anx
dev.off()

pdf(paste0(here("out//supplementary//"), "int_t_beta_ecz_dep.pdf"), 8, 8)
  plot_ecz_dep
dev.off()

pdf(paste0(here("out//supplementary//"), "int_t_beta_pso_anx.pdf"), 8, 8)
  plot_pso_anx
dev.off()

pdf(paste0(here("out//supplementary//"), "int_t_beta_pso_dep.pdf"), 8, 8)
  plot_pso_dep
dev.off()





# newdata_plot <- data.frame(
#   exposed = c("Unexposed", "Psoriasis"),
#   t = rep(365.25,2),
#   carstairs = rep("3",2),
#   cal_period = rep("2008-14",2),
#   comorbid = rep("No",2),
#   cci = rep("Low",2),
#   bmi2 = rep(0,2), 
#   alc = rep("No",2),
#   smokstatus = rep("Non-Smoker",2)
# )
# test <- survfit(smpl_fit3_int, newdata = newdata_plot)
# test_fit <- summary(test)
# 
# 
# 

