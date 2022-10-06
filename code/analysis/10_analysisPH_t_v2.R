source(here::here("code/packages.R"))
source(here::here("code/file_paths.R"))

source(here("code/programs/schonfeld_plot.R")) ## adapted plot code for Scales Schoenfeld residuals
dir.create(file.path(here("out")), showWarnings = FALSE)
dir.create(file.path(here("out", "PHchecks")), showWarnings = FALSE)

YY <- c("depression", "anxiety")
XX <- c("psoriasis", "eczema")

exposure <- XX[1]
outcome <- YY[1]


pdf(paste0(here::here("out/analysis"), "/fig2_spline_time_estimates_withBootstrap.pdf"), 8, 8)
par(mfrow = c(2,2))
ii=1
matOut <- matrix(NA, nrow = 4, ncol = 4*3)
matj <- 1
for(exposure in XX) {
  ABBRVexp <- substr(exposure, 1, 3)
  for(outcome in YY) {
    # load data ---------------------------------------------------------------
    df_model <- readRDS(paste0(datapath, "out/df_model", ABBRVexp, "_", outcome,".rds"))
    .dib(exposure)
    
    ## fit model with time interaction 
    # fit1_int <-
    #   coxph(
    #     Surv(t, out) ~ exposed * t + strata(setid),
    #     data = df_model
    #   ) 
    # fit2_int <-
    #   coxph(
    #     Surv(t, out) ~ exposed * t + carstairs + cal_period + comorbid + cci + strata(setid),
    #     data = df_model
    #   ) 
    # fit3_int <-
    #   coxph(
    #     Surv(t, out) ~ exposed * t + carstairs + cal_period + comorbid + cci + bmi2 + alc + smokstatus + strata(setid),
    #     data = df_model
    #   ) 
    # 
    # time_series <- seq(0, 365*6, 365.25/12)
    # interactions <- c(
    #   paste0("exposed", str_to_title(exposure), " + ", time_series, "*exposed", str_to_title(exposure), ":t = 0")
    # )
    # 
    # lincom_out_min <- lincom(fit1_int, 
    #                          interactions, eform = T, level = 0.95)
    # lincom_out_con <- lincom(fit2_int, 
    #                          interactions, eform = T, level = 0.95)
    # lincom_out_med <- lincom(fit3_int, 
    #                          interactions, eform = T, level = 0.95)
    # 
    # 
    # lincom_out_min <- as_tibble(lincom_out_min, rownames = "lincom") %>% mutate(model = "Minimal")
    # lincom_out_con <- as_tibble(lincom_out_con, rownames = "lincom") %>% mutate(model = "Confounder")
    # lincom_out_med <- as_tibble(lincom_out_med, rownames = "lincom") %>% mutate(model = "Mediator")
    # tibble_int <- lincom_out_min %>%
    #   bind_rows(lincom_out_con) %>% 
    #   bind_rows(lincom_out_med) %>% 
    #   janitor::clean_names() %>%
    #   dplyr::select(model, 
    #                 lincom,
    #                 estimate,
    #                 conf.low = x2_5_percent,
    #                 conf.high = x97_5_percent) %>%
    #   mutate_at(c("estimate", "conf.low", "conf.high"), ~ unlist(.)) %>% 
    #   mutate(y = outcome, x = exposure, z = "t")
    # 
    # tibble_plot <- tibble_int %>% 
    #   mutate(t = rep(time_series, 3)) %>% 
    #   mutate(yr = t/365.25)
    # 
    # plot_int <- ggplot(tibble_plot, aes(x = yr, y = estimate, ymin = conf.low, ymax = conf.high, group = model, fill = model, col = model)) + 
    #   geom_line() + 
    #   geom_ribbon(lty = 0, alpha = 0.2) +
    #   geom_hline(yintercept = 1, lty = 4, col = 1) +
    #   ylim(c(0,2)) +
    #   labs(x = "Year", 
    #        y = expression(hat(beta) ~ "exposure:t"), 
    #        title = paste0(str_to_title(exposure), " ~ ", str_to_title(outcome))) +
    #   guides(fill = guide_legend(title = "Model"), col = guide_legend(title = "Model")) +
    #   theme_ali()
    # assign(paste0("interaction_plot", ABBRVexp, substr(outcome,1,3)), plot_int)

    ## Spline
    cox_fit <- readRDS(paste0(datapath, "out/models_data/", ABBRVexp, "_", outcome, "_mod3_modeldata.rds"))
    
    df_model$exp <- as.numeric(df_model$exposed)-1 ## need a numeric exposure variable for tt() to work
    
    ## Therneau and time dependent variables 
    mediator_est <- broom::tidy(cox_fit, conf.int = T, conf.level = 0.95, exp = T) %>% slice(1)
    
    ## Bootstrap to get a confidence interval
    output <- data.frame(fup = seq(min(df_model$exp + df_model$t) / 365.25,
                                   max(df_model$exp + df_model$t) / 365.25,
                                   0.01))
    ## get time status at 0, 1, 3, 5 years
    get_match <- function(xx){
      which.min(abs(output$fup-xx))
    }
    matching_row_indices <- sapply(c(0,1,3,5), get_match)
    
    b <- 200
    data_ids <- unique(df_model$setid) 
    n <- length(data_ids)
    
    output_Pspline <- matrix(NA, nrow = b, ncol = length(output$fup))
    
    tstart <- Sys.time()
    set.seed(105351)
    for(btsp in 1:b){
      if(btsp %% 10 == 0){print(paste("run", btsp, "..."))}
      
      if(btsp == 1){ # for the first run just use the actual data 
        sample_data <- df_model
      }else{ # then sample and jumble up the data on each iteration of the bootstrap 
        sample <- sample(data_ids, size = n, replace = TRUE) 
        sample_ids <- tibble(setid = sample) %>% arrange(setid)
        sample_data <- sample_ids %>% 
          left_join(df_model, by = "setid")
      }
      
      # Penalised spline interaction model for time varying HR
      kk <- c(1, 2, 3)
      if (ABBRVexp == "ecz") {
        pspline_model <-
          coxph(
            Surv(t, out) ~ exp + tt(exp) + carstairs + cal_period + comorbid + cci + bmi_cat + sleep + alc + smokstatus + gc90days + strata(setid),
            data = sample_data,
            tt = function(x, t, ...) x * pspline(t/365.25)
          )
      } else if (ABBRVexp == "pso") {
        pspline_model <-
          coxph(
            Surv(t, out) ~ exp + tt(exp) + carstairs + cal_period + comorbid + cci + bmi_cat + alc + smokstatus + strata(setid),
            data = sample_data,
            tt = function(x, t, ...) x * pspline(t/365.25)
          ) 
      }
      ## get time-varying estimate of HR
      pspline_time <- pspline(output$fup)
      pspline_time_coeffs <- pspline_model$coefficients[str_detect(string = names(pspline_model$coefficients), pattern = "ps\\(")]
      out_Pspline <- pspline_model$coefficients[1] + (pspline_time %*% pspline_time_coeffs)
    
      output_Pspline[btsp, ] <- out_Pspline
    }
    tstop <- Sys.time()
    print(tstop-tstart)
    
    medP <- apply(output_Pspline,2,function(x){median(x, na.rm=T)}) %>% exp()
    ciP1 <- apply(output_Pspline,2,function(x){quantile(x,0.025, na.rm=T)}) %>% exp()
    ciP2 <- apply(output_Pspline,2,function(x){quantile(x,0.975, na.rm=T)}) %>% exp()
    
    main_Pspline <- output_Pspline[1,] %>% exp()
    
    ## get predicted HR at 0, 1, 3 and 5 years
    medTab <- medP[matching_row_indices]
    ci1Tab <- ciP1[matching_row_indices]
    ci2Tab <- ciP2[matching_row_indices]
    
    matOut[matj,] <- c(medTab, ci1Tab, ci2Tab)
    matj <- matj+1
    
    ## PLOT
    ymax <- max(ciP2)
    plot(range(output$fup), range(c(main_Pspline)), type = "n",
         xlab="Time (in years)", ylab="", ylim = c(0.8, 1.4))
    lines(range(output$fup), rep(mediator_est$estimate,2), lwd = 2, lty = 1, col = 1)
    polygon(c(range(output$fup), rev(range(output$fup))),
            c(rep(mediator_est$conf.low, 2),rep(mediator_est$conf.high, 2)),
            col = ggplot2::alpha(1, 0.2), lty = 0)
    mtext(paste0("HR(t) for ", exposure), side = 2, padj = -4, cex = 0.9)
    # pspline
    lines(output$fup, medP, col = 2, lty = 2)
    polygon(c(output$fup, rev(output$fup)),
            c(ciP1,rev(ciP2)),
            col = ggplot2::alpha(2, 0.2), lty = 0)
    
    abline(h = 1, lwd = 2, lty = 2, col = 9)
    legend(0, 0.99, legend = c("Proportional hazards", "Penalised spline"), col = c(1, 2), lty = 1, bty = "n")
    mtext(paste0(LETTERS[ii], ": ", str_to_title(exposure), " ~ ", str_to_title(outcome)), side=3, line=2, col=1, cex=1, font=2, adj = 0)
    ii <- ii+1
    
    saveRDS(output_Pspline, file = paste0(here::here("out/data/"), "/outputPspline_", exposure, "_", outcome, ".rds"))
  }
}
dev.off()

outputPspline_psoriasis_depression <- readRDS(here::here("out/data/outputPspline_psoriasis_depression.rds"))
outputPspline_psoriasis_anxiety <- readRDS(here::here("out/data/outputPspline_psoriasis_anxiety.rds"))
outputPspline_eczema_depression <- readRDS(here::here("out/data/outputPspline_eczema_depression.rds"))
outputPspline_eczema_anxiety <- readRDS(here::here("out/data/outputPspline_eczema_anxiety.rds"))

saveRDS(object = list(
  outputPspline_psoriasis_depression,
  outputPspline_psoriasis_anxiety,
  outputPspline_eczema_depression,
  outputPspline_eczema_anxiety),
  file = here::here("out/data/list_penalisedspline_bootstrapsamples.rds")
  )


## need to put exposure, outcome and time periods on this and save, then output as nice table 
matOut
dfOut <- as.data.frame(matOut)

colnames(dfOut) <- c(paste0("t", c(0,1,3,5), "med"),paste0("t", c(0,1,3,5), "ci1"),paste0("t", c(0,1,3,5), "ci2"))
dfOut$x <- c(rep(XX[1], 2), rep(XX[2], 2)) %>% str_to_title()
dfOut$y <- c(rep(YY, 2)) %>% str_to_title()

plotOut <- dfOut %>% 
  pivot_longer(cols = -c(x, y), names_pattern = "t(.)(.*)", names_to = c("year", ".value"))

pd <- position_dodge(width = 0.5)
ggplot(plotOut, aes(x = year, y = med, ymin = ci1, ymax = ci2, colour = x)) +
  geom_errorbar(position = pd, width = 0.25) +
  geom_point(position = pd, size = 3, shape = 1) +
  geom_hline(yintercept = 1, lty = 2, col = 1) +
  facet_wrap(~y) + 
  scale_y_log10(breaks=seq(0.5,2,0.1),limits=c(0.8,1.35)) +
  #ylim(c(0.8, NA)) +
  labs(y = "Hazard ratio", x = "Year since diagnosis", colour = "Exposure") +
  theme_ali() +
  theme(strip.background = element_blank())
ggsave(filename = here::here("out/supplementary/timevaryingHRestimates.pdf"), width = 6, height = 6)

tabOut <- dfOut
gtOut <- tabOut %>% 
  dplyr::select(x, y, everything()) %>% 
  gt() %>% 
  gt::fmt_number(columns = where(is.numeric), decimals = 2) %>% 
  cols_merge(columns = c(3, 7, 11), pattern = "{1} ({2}-{3})") %>% 
  cols_merge(columns = c(4, 8, 12), pattern = "{1} ({2}-{3})") %>% 
  cols_merge(columns = c(5, 9, 13), pattern = "{1} ({2}-{3})") %>% 
  cols_merge(columns = c(6, 10, 14), pattern = "{1} ({2}-{3})") %>% 
  cols_label(
    x = "Exposure", 
    y = "Outcome",
    t0med = "At diagnosis",
    t1med = "+ 1 year",
    t3med = "+ 3 years",
    t5med = "+ 5 years"
  )
gtOut
write.csv(gtOut$`_data`, here::here("out/supplementary/df_timevaryingHR_estimates.csv"))
gtOut %>% gt::gtsave(filename = "tab15_timevaryingHR_estimates.html", path = here::here("out/tables/"))

# 
# pdf(paste0(here::here("out/analysis"), "/fig3_linear_time_estimates.pdf"), 8, 8)
# p1 <- cowplot::plot_grid(
#   interaction_ploteczanx,
#   interaction_ploteczdep,
#   interaction_plotpsoanx,
#   interaction_plotpsodep
# )
# print(p1)
# dev.off()
