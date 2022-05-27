library(tidyverse)
library(here)
library(magrittr)
library(survival)
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
exposure = XX[1]
outcome = YY[1]

final_calresults <- NULL

for(exposure in XX){
  ABBRVexp <- substr(exposure, 1, 3)
  for (outcome in YY){
    # load df_model -----------------------------------------------------------
    df_model <-
      readRDS(paste0(datapath, "out/df_model", ABBRVexp, "_", outcome,".rds"))
    
    # make a nice series of the days of each month in the study period (1997-2020)
    leapyearsindata <- leap_year(seq.Date(as.Date("1997-01-01"), as.Date("2019-01-01"), "1 year"))
    little_fn <- function(leap = FALSE){
      oneyear <- c(as.Date("1997-01-01"),as.Date("1997-02-01"),as.Date("1997-03-01"),as.Date("1997-04-01"),as.Date("1997-05-01"),as.Date("1997-06-01"),as.Date("1997-07-01"),as.Date("1997-08-01"),as.Date("1997-09-01"),as.Date("1997-10-01"),as.Date("1997-11-01"),as.Date("1997-12-01"))
      oneyear_leap <- c(as.Date("2000-01-01"),as.Date("2000-02-01"),as.Date("2000-03-01"),as.Date("2000-04-01"),as.Date("2000-05-01"),as.Date("2000-06-01"),as.Date("2000-07-01"),as.Date("2000-08-01"),as.Date("2000-09-01"),as.Date("2000-10-01"),as.Date("2000-11-01"),as.Date("2000-12-01"))
      oneyeardays <- days_in_month(oneyear) 
      oneyearleapdays <- days_in_month(oneyear_leap) 
      if(leap){
        oneyearleapdays
      }else{
        oneyeardays
      }
    }
    split_series <- sapply(leapyearsindata, little_fn) %>% as.vector() %>% cumsum()
    split_series <- c(0, split_series)
    names(split_series) <- as.Date("1997-01-01")+split_series
    # as.Date("1997-01-01") + max(split_series) # check that this is 1jan2020
    
    temp <- df_model %>%
      filter(exposed == str_to_title(exposure)) %>% 
      select(setid, patid, exposed, dob, indexdate, tstart, tstop, out) %>%
      mutate(eligible_from = dob + tstart,
             eligible_to = dob + tstop)
    temp$tstartdate <- temp$dob + temp$tstart
    temp$tstopdate <- temp$dob + temp$tstop
    
    temp$tstartrel97 <- as.numeric(temp$tstartdate - as.Date("1997-01-01"), units = "days")
    temp$tstoprel97 <- as.numeric(temp$tstopdate - as.Date("1997-01-01") , units = "days")
    
    #Convert tstart and tstop back to dates
    df_split_cal <- survSplit(
      Surv(
        time = as.numeric(tstartrel97),
        time2 = as.numeric(tstoprel97),
        event = out
      ) ~ .,
      data = temp,
      cut = split_series,
      episode = "month"
    )
    months <- as.Date("1997-01-01")+split_series %>% as.vector()
    months_merge <- data.frame(
      x = 1:length(split_series),
      m = month(months)
    )
    df_split_cal12 <- df_split_cal %>% 
      left_join(months_merge, by = c("month"="x"))
    rm(df_split_cal, temp)
    
    head(df_split_cal12)
    
    sum_cal12 <- df_split_cal12 %>% 
      group_by(m) %>% 
      summarise(at_risk = n(), outcome = sum(out)) %>% 
      mutate(prop = outcome / at_risk) %>% 
      ungroup() %>%
      rowwise() %>%
      mutate(tst = list(broom::tidy(
        prop.test(outcome, at_risk, conf.level = 0.95)
      ))) %>%
      tidyr::unnest(tst)
    
    sum_cal12$month2 <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
    sum_cal12$exp <- exposure
    sum_cal12$out <- outcome
    
    final_calresults <- bind_rows(
      final_calresults,
      sum_cal12
    )
  }
}

final_calresults
final_calresults$month3 <- factor(final_calresults$month2, levels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
dodge <- position_dodge(width = 0.9)
p2 <-
  ggplot(
    final_calresults,
    aes(
      x = month3,
      y = estimate * 100,
      ymin = conf.low * 100,
      ymax = conf.high * 100,
      group = out, 
      col = out
    )
  ) +
  geom_point(pch = 16, position = dodge) +
  geom_errorbar(lty = 1,
                alpha = 0.8,
                position = dodge) +
  ylim(c(0,0.3)) +
  labs(y = "Incidence",
       x = "Month",
       col = "Outcome") +
  facet_wrap(~exp, ncol = 2) +
  theme_ali() +
  theme(axis.text.x = element_text(angle = 30),
        strip.background = element_blank())
p2
pdf(here::here("out/supplementary",paste0("month-prevalence-v2.pdf")), 10, 6)
  print(p2)
dev.off()

# fit a model with month in it and see what happens -----------------------
tibble_out <- NULL
for(exposure in XX){
  ABBRVexp <- substr(exposure, 1, 3)
  for (outcome in YY){
    # load df_model -----------------------------------------------------------
    df_model <-
      readRDS(paste0(datapath, "out/models_data/df_model", ABBRVexp, "_", outcome,".rds"))
    ## load model 2 (confounder adjusted)
    original_confounder <- readRDS(paste0(datapath, "out/models_data/", ABBRVexp, "_", outcome,"_mod2_modeldata.rds"))
    
    df_model$month <- factor(lubridate::month(df_model$dob + df_model$tstop), labels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
    
    
    .dib("Running model 2 (confounder + month)")
    mod2 <-
      coxph(Surv(t, out) ~ exposed*month + carstairs + cal_period + comorbid + cci + strata(setid),
            data = df_model)
    broom::tidy(mod2, exp = T)
    broom::tidy(original_confounder, exp = T)
    coeffs <- mod2$coefficients %>% names()
    int_var <- "month"
    int_levels <- df_model[, int_var] %>% levels()
    n_int_levels <- length(int_levels) - 1 # -1 because of reference category
    coeffs_interaction <- coeffs[str_detect(coeffs, int_var)]
    
    interactions <- c(coeffs[1])
    for (ii in 1:n_int_levels) {
      X <-
        paste0(c(coeffs[1], coeffs_interaction[ii + n_int_levels]),
               collapse = "+")
      interactions <- c(interactions, X)
    }
    interactions
    lincom_out <- biostat3::lincom(mod2,
                         interactions,
                         eform = TRUE,
                         level = 0.99)
    rownames(lincom_out)[1] <-
      paste0(rownames(lincom_out)[1], "+", int_var, int_levels[1])
    
    lincom_out <- as_tibble(lincom_out, rownames = "lincom")
    tibble_int <- lincom_out %>%
      janitor::clean_names() %>%
      dplyr::select(lincom,
                    estimate,
                    conf.low = x0_5_percent,
                    conf.high = x99_5_percent) %>%
      mutate_at(c("estimate", "conf.low", "conf.high"), ~ unlist(.)) %>% 
      mutate(y = outcome, x = exposure, z = "month") 
    tibble_int$month <- factor(levels(df_model$month), levels = levels(df_model$month))
    tibble_out <- tibble_out %>% 
      bind_rows(tibble_int)
  }
}

tibble_out$y <- str_to_title(tibble_out$y)
tibble_out$x <- str_to_title(tibble_out$x)
ggplot(tibble_out, aes(x = month, y = estimate, ymin = conf.low, ymax = conf.high, col = month)) +
  geom_errorbar() +
  geom_point() +
  facet_grid(y~x) +
  labs(y = "HR", x = "Month of year") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic"),
        legend.position = "none")

dev.copy(pdf, here::here("out/supplementary/forest_month_interaction.pdf"), width = 8, height = 6); dev.off()
