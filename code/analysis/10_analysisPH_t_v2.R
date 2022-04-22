library(tidyverse)
library(here)
library(survival)
library(lmtest)
library(broom)
library(biostat3)
library(multcomp)
library(rms)
library(splines)

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

exposure <- XX[2]
outcome <- YY[2]

for(exposure in XX) {
  ABBRVexp <- substr(exposure, 1, 3)
  for(outcome in YY) {
    # load data ---------------------------------------------------------------
    df_model <- readRDS(paste0(datapath, "out/models_data/df_model", ABBRVexp, "_", outcome,".rds"))
    
    ## fit model with time interaction 
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
    
    plot_int <- ggplot(tibble_plot, aes(x = yr, y = estimate, ymin = conf.low, ymax = conf.high, group = model, fill = model, col = model)) + 
      geom_line() + 
      geom_ribbon(lty = 0, alpha = 0.2) +
      geom_hline(yintercept = 1, lty = 4, col = 1) +
      ylim(c(0,2)) +
      labs(x = "Year", 
           y = expression(hat(beta) ~ "exposure:t"), 
           title = paste0(str_to_title(exposure), " ~ ", str_to_title(outcome))) +
      guides(fill = guide_legend(title = "Model"), col = guide_legend(title = "Model")) +
      theme_ali()
    
    pdf(paste0(here::here("out/PHchecks/"), "linear_time_", ABBRVexp, substr(outcome,1,3),".pdf"), 6, 6)
    plot_int
    dev.off()
    
    ## Spline
    df_model$exp <- as.numeric(df_model$exposed)-1
    
    ## Therneau and time dependent variables 
    cox_fit <- readRDS(paste0(datapath, "out/models_data/", ABBRVexp, "_", outcome, "_mod1_modeldata.rds"))
    
    zpI1 <- cox.zph(cox_fit, transform = "identity")
    plot(zpI1, resid = F)
    
    # zp1 <- cox.zph(cox_fit, transform = function(time) log(time + 20))
    # plot(zp1, resid = F)
    # terry <- coxph(Surv(t, out) ~ exp + tt(exp) + strata(setid),
    #                data = df_model,
    #                tt = function(x, t, ...) x * log(t+20))
    # terry
    # plot(zp[1], col = 12, lwd = 2, resid = F)
    # abline(h = coef(cox_fit)[1], lwd = 2, lty = 1, col = 9)
    # abline(coef(terry)[1:2], lwd = 2, lty = 1, col = 2)
    # 
    ## try with pspline
    
    pspline <- coxph(Surv(t, out) ~ exp + tt(exp) + strata(setid), 
                     data = df_model, 
                     tt = function(x, t, ...) x * pspline(t/365.25))
    
    pspline$coefficients[1] %>% exp()
    pspline$coefficients[-1]
    basis_fn <- length(pspline$coefficients[!is.na(pspline$coefficients)])
    output <- data.frame(fup = seq(min(df_model$exp + df_model$t) / 365.25,
                                   max(df_model$exp + df_model$t) / 365.25,
                                   0.01))
    
    pspline_time <- pspline(output$fup)
    output$HR <- pspline$coefficients[1] + (pspline_time %*% pspline$coefficients[-1])
    
    ## try with splines::bs
    kk <- c(1, 2, 3)
    spline <- coxph(Surv(t, out) ~ exp + tt(exp) + strata(setid), 
                    data = df_model, 
                    tt = function(x, t, ...) x * splines::bs(t/365.25, degree = 3, knots = kk))
    
    spline %>% summary()
    
    spline$coefficients
    basis_fn <- length(spline$coefficients[!is.na(spline$coefficients)])-1
    output <- data.frame(fup = seq(min(df_model$exp) + min(df_model$t) / 365.25,
                                   max(df_model$exp + df_model$t) / 365.25,
                                   0.01))
    spline_time <- splines::bs(output$fup, degree = 3, knots = kk)
    
    
    output$HRspline <- spline$coefficients[1] + (spline_time %*% spline$coefficients[-1])
    
    pdf(paste0(here::here("out/PHchecks/"), "spline_time_", ABBRVexp, substr(outcome,1,3),".pdf"), 8, 8)
    par(mfrow = c(2,1))
    plot_schonfeld(zpI[1], col = "darkgreen", df = 5, se = F,
                   lwd = 1.5, resid = F, xlab = "Time (in days)", hr = T, ylab = "")
    mtext(expression(e^{hat(beta)(t)} ~ "for" ~ exposed), side = 2, padj = -2, cex = 0.7)
    abline(h = exp(coef(cox_fit)[1]), lwd = 2, lty = 2, col = 2)
    abline(h = 1, lwd = 2, lty = 2, col = 9)
    
    plot(range(output$fup), range(c(exp(output$HR),exp(output$HRspline))), type = "n", log = "y",
         xlab="Time (in years)", ylab="")
    mtext(expression(e^{hat(beta)(t)} ~ "for" ~ exposed), side = 2, padj = -2, cex = 0.7)
    lines(output$fup, exp(output$HR))
    lines(output$fup, exp(output$HRspline), col = 4)
    abline(h = exp(coef(cox_fit)[1]), lwd = 2, lty = 2, col = 2)
    abline(h = 1, lwd = 2, lty = 2, col = 9)
    legend(0, 1.2, legend = c("Penalised spline", "Basis spline (1,2,3 years)"), col = c(1,4), lty = 1)
    dev.off()
  }
}