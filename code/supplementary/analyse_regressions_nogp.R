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
    datapath <- "/Users/lsh1510922/Documents/Postdoc/2021_extract/"
  }
  if (Sys.info()["sysname"] == "Windows") {
    datapath <- "Z:/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
}
dir.create(file.path(here("out", "nogp")), showWarnings = FALSE)
dir.create(file.path(here("out", "nogp", "analysis")), showWarnings = FALSE)
dir.create(file.path(datapath, "out/nogp/"), showWarnings = FALSE)
dir.create(file.path(datapath, "out/nogp/models_data"), showWarnings = FALSE)

XX <- c("psoriasis")
YY <- c("anxiety", "depression")

export_plots <- F

st_time <- Sys.time()
for (exposure in XX) {
  #exposure <- XX[1]
  ABBRVexp <- str_sub(exposure, 1 , 3)
    df_anx_split <-
      readRDS(paste0(datapath, "out/nogp/pso-anxiety_split.rds"))
    df_dep_split <-
      readRDS(paste0(datapath, "out/nogp/pso-depression_split.rds"))
  
  .dib(exposure)
  
  df_anx_split$t <-
    as.numeric(df_anx_split$tstop - df_anx_split$tstart)
  df_dep_split$t <-
    as.numeric(df_dep_split$tstop - df_dep_split$tstart)
  
  # annual period prevalence plot -------------------------------------------
  if (export_plots) {
    temp <- df_anx_split %>%
      select(setid, patid, exposed, dob, indexdate, tstart, tstop, out) %>%
      mutate(eligible_from = dob + tstart,
             eligible_to = dob + tstop)
    
    #Convert tstart and tstop back to dates
    df_split_cal <- survSplit(
      Surv(
        time = as.numeric(eligible_from),
        time2 = as.numeric(eligible_to),
        event = out
      ) ~ .,
      data = temp,
      cut = as.numeric(as.Date("1996-01-02") + c(365.25 * 1:23)),
      episode = "year"
    )
    
    
    df_split_cal <- df_split_cal %>%
      group_by(setid, patid, year) %>%
      slice(1)
    df_split_cal$year2 <- factor(df_split_cal$year,
                                 levels = 1:24,
                                 labels = as.character(1997:2020))
    df_split_sum <- df_split_cal %>%
      group_by(exposed, year2) %>%
      summarise(at_risk = n(), outcome = sum(out)) %>%
      mutate(prop = outcome / at_risk,
             exposed = factor(
               exposed,
               levels = 0:1,
               labels = c("Matched controls", paste0("With ", exposure))
             )) %>%
      ungroup() %>%
      rowwise() %>%
      mutate(tst = list(broom::tidy(
        prop.test(outcome, at_risk, conf.level = 0.95)
      ))) %>%
      tidyr::unnest(tst)
    
    p1 <-
      ggplot(
        df_split_sum,
        aes(
          x = year2,
          y = estimate * 100,
          ymin = conf.low * 100,
          ymax = conf.high * 100,
          group = exposed,
          colour = exposed,
          fill = exposed
        )
      ) +
      geom_line(lty = 2) +
      geom_ribbon(lty = 0, alpha = 0.2) +
      labs(title = paste0(exposure, " ~ anxiety"),
           y = "Period prevalence",
           x = "Year") +
      theme_ali() +
      theme(axis.text.x = element_text(angle = 30))
    
    dodge <- position_dodge(width = 0.9)
    p2 <-
      ggplot(
        df_split_sum,
        aes(
          x = year2,
          y = estimate * 100,
          ymin = conf.low * 100,
          ymax = conf.high * 100,
          group = exposed,
          colour = exposed,
          fill = exposed
        )
      ) +
      geom_point(pch = 16, position = dodge) +
      geom_errorbar(lty = 1,
                    alpha = 0.8,
                    position = dodge) +
      labs(title = paste0(exposure, " ~ anxiety"),
           y = "Period prevalence",
           x = "Year") +
      theme_ali() +
      theme(axis.text.x = element_text(angle = 30))
    
    
    temp <- df_dep_split %>%
      select(setid, patid, exposed, dob, indexdate, tstart, tstop, out) %>%
      mutate(eligible_from = dob + tstart,
             eligible_to = dob + tstop)
    
    #Convert tstart and tstop back to dates
    df_split_cal <- survSplit(
      Surv(
        time = as.numeric(eligible_from),
        time2 = as.numeric(eligible_to),
        event = out
      ) ~ .,
      data = temp,
      cut = as.numeric(as.Date("1996-01-02") + c(365.25 * 1:23)),
      episode = "year"
    )
    
    
    df_split_cal <- df_split_cal %>%
      group_by(setid, patid, year) %>%
      slice(1)
    df_split_cal$year2 <- factor(df_split_cal$year,
                                 levels = 1:24,
                                 labels = as.character(1997:2020))
    df_split_sum <- df_split_cal %>%
      group_by(exposed, year2) %>%
      summarise(at_risk = n(), outcome = sum(out)) %>%
      mutate(prop = outcome / at_risk,
             exposed = factor(
               exposed,
               levels = 0:1,
               labels = c("Matched controls", paste0("With ", exposure))
             )) %>%
      ungroup() %>%
      rowwise() %>%
      mutate(tst = list(broom::tidy(
        prop.test(outcome, at_risk, conf.level = 0.95)
      ))) %>%
      tidyr::unnest(tst)
    
    p3 <-
      ggplot(
        df_split_sum,
        aes(
          x = year2,
          y = estimate * 100,
          ymin = conf.low * 100,
          ymax = conf.high * 100,
          group = exposed,
          colour = exposed,
          fill = exposed
        )
      ) +
      geom_line(lty = 2) +
      geom_ribbon(lty = 0, alpha = 0.2) +
      labs(
        title = paste0(exposure, " ~ depression"),
        y = "Period prevalence",
        x = "Year"
      ) +
      theme_ali() +
      theme(axis.text.x = element_text(angle = 30))
    
    
    #dodge <- position_dodge(width=0.9)
    p4 <-
      ggplot(
        df_split_sum,
        aes(
          x = year2,
          y = estimate * 100,
          ymin = conf.low * 100,
          ymax = conf.high * 100,
          group = exposed,
          colour = exposed,
          fill = exposed
        )
      ) +
      geom_point(pch = 16, position = dodge) +
      geom_errorbar(lty = 1,
                    alpha = 0.8,
                    position = dodge) +
      labs(
        title = paste0(exposure, " ~ depression"),
        y = "Period prevalence",
        x = "Year"
      ) +
      theme_ali() +
      theme(axis.text.x = element_text(angle = 30))
    
    pAll <- cowplot::plot_grid(p2, p4, ncol = 1)
    
    pdf(here::here(
      "out/nogp/analysis",
      paste0("ann-prevalence-", ABBRVexp, ".pdf")
    ), 10, 8)
    print(pAll)
    dev.off()
    
    rm(df_split_sum, df_split_cal, temp, p1, p2, p3, p4)
  }
  
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
    #rename(!!paste0(exposure) := exposed) %>% 
    df_model$exposed <- factor(df_model$exposed, levels = c("Unexposed", str_to_title(exposure)))
    var_label(df_model$exposed) <- "Exposure"
    var_label(df_model$agegroup) <- "Age group"
    df_model$agegroup <- relevel(df_model$agegroup, ref = "50-59")
    
    saveRDS(df_model, file = paste0(datapath, "out/nogp/models_data/df_model", ABBRVexp, "_", outcome, ".rds"))
    
    # Run a simple Cox regression ----------------------------------------
    .dib("Running model 1 (crude)")
    mod1 <-
      coxph(Surv(t, out) ~ exposed + strata(setid), data = df_model) #%>%
    
    # Run a confounder Cox regression ----------------------------------------
    # .dib("Running model 2 (confounder)")
    # mod2 <-
    #   coxph(Surv(t, out) ~ exposed + carstairs + cal_period + comorbid + cci + strata(setid),
    #         data = df_model)
    
    saveRDS(
      mod1,
      file = paste0(
        datapath,
        "out/nogp/models_data/",
        ABBRVexp,
        "_",
        outcome,
        "_mod1_modeldata.rds"
      )
    )
    
  }
}
end_tim <- Sys.time()
end_tim - st_time



# compare to gp-matched results -------------------------------------------
exposure <- "psoriasis"
ABBRVexp <- str_sub(exposure, 1 , 3)

mod1_gp <- readRDS(paste0(datapath,"out/models_data/",ABBRVexp,"_",outcome,"_mod1_modeldata.rds"))
mod1_nogp <- readRDS(paste0(datapath,"out/nogp/models_data/",ABBRVexp,"_",outcome,"_mod1_modeldata.rds"))

summary(mod1_gp, conf.level = 0.99)
summary(mod1_nogp, conf.level = 0.99)

