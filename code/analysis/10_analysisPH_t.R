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
    
    assign(paste0("plot_", ABBRVexp, "_", substr(outcome, 1, 3)), plot)
    
    ## interaction total fup
    df_model_fup <- df_model %>% 
      group_by(patid) %>% 
      mutate(rownum = 1:n(), 
             sumfup = cumsum(t)) %>% 
      ungroup()

    fit1_int_fup <-
      coxph(
        Surv(t, out) ~ exposed * sumfup + strata(setid),
        data = df_model_fup
      ) 
    fit2_int_fup <-
      coxph(
        Surv(t, out) ~ exposed * sumfup + carstairs + cal_period + comorbid + cci + strata(setid),
        data = df_model_fup
      ) 
    fit3_int_fup <-
      coxph(
        Surv(t, out) ~ exposed * sumfup + carstairs + cal_period + comorbid + cci + bmi2 + alc + smokstatus + strata(setid),
        data = df_model_fup
      ) 
    
    time_series <- seq(0, 365*20, 365.25/6)
    interactions <- c(
      paste0("exposed", str_to_title(exposure), " + ", time_series, "*exposed", str_to_title(exposure), ":sumfup = 0")
    )
    
    lincom_out_min <- lincom(fit1_int_fup,
                             interactions, eform = T, level = 0.95)
    lincom_out_con <- lincom(fit2_int_fup,
                             interactions, eform = T, level = 0.95)
    lincom_out_med <- lincom(fit3_int_fup,
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
      mutate(sumfup = rep(time_series, 3)) %>% 
      mutate(yr = sumfup/365.25)
    
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
    
    ## Spline 
    df_model_spline <- df_model %>% 
      group_by(patid) %>% 
      mutate(rownum = 1:n(), 
             sumfup = cumsum(t)) %>% 
      ungroup()
    create_splinetime <- splines::bs(df_model_spline$sumfup, degree = 3, knots = c(365, 365*3, 365*5)) %>%
      as.data.frame() %>% 
      janitor::clean_names()
  
    test <- df_model_spline %>% 
      bind_cols(create_splinetime) 
    test %>% 
      slice(1:10) %>% 
      dplyr::select(patid, t, sumfup, starts_with("x")) %>% 
      View
    
    test2 <- test %>% 
      slice(1:1000) %>% 
      arrange(sumfup)
    plot(test2$sumfup,  test2$x1, type = 'l', col = 1, ylim = c(0, 1))
    lines(test2$sumfup, test2$x2, type = 'l', col = 2)
    lines(test2$sumfup, test2$x3, type = 'l', col = 3)
    lines(test2$sumfup, test2$x4, type = 'l', col = 4)
    lines(test2$sumfup, test2$x5, type = 'l', col = 5)
    lines(test2$sumfup, test2$x6, type = 'l', col = 6)
    
    test3 <- df_model_fup  %>% 
      slice(1:100000)
    fit1_int_fup <-
      coxph(
        Surv(t, out) ~ exposed * splines::bs(t, degree = 3, knots = c(365, 1095)) + strata(setid),
        data = test3
      ) 
    fit1_int_fup
    
    
    library(Greg)
    test2$x1
    fit.coxph <- coxph(Surv(t, out) ~ exposed*x1 + exposed*x2 + exposed*x3 + exposed*x4 + exposed*x5+ exposed*x6 , data = test2)
    predict(fit.coxph, exposed)
    
    
    
    
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



n <- 1000
set.seed(731)
age <- 50 + 12*rnorm(n)
label(age) <- "Age"
sex <- factor(sample(c('Male','Female'), n, 
                     rep=TRUE, prob=c(.6, .4)))
cens <- 15*runif(n)
h <- .02*exp(.04*(age-50)+.8*(sex=='Female'))
dt <- -log(runif(n))/h
label(dt) <- 'Follow-up Time'
e <- ifelse(dt <= cens,1,0)
dt <- pmin(dt, cens)
units(dt) <- "Year"
dd <- datadist(age, sex)
options(datadist='dd')
S <- Surv(dt,e)

f <- cph(S ~ rcs(age,4) + sex, x=TRUE, y=TRUE)
cox.zph(f, "rank")             # tests of PH
anova(f)
plot(Predict(f, age, sex)) # plot age effect, 2 curves for 2 sexes


#install.packages(:)
library(rstpm2)
test2
mod_tvc <- stpm2(Surv(t,out==1)~exposed,data=test2,df=3)

plot(mod_tvc, newdata = data.frame(exposed = 1), type = "hr",
     var = "arm", ci = TRUE, rug = FALSE,
     main = "Time dependent hazard ratio", xlim=c(1,24), ylim=c(0,2.5),
     ylab = "Hazard ratio", xlab = "Time"
)


library(mfp)
library(visreg)

data(GBSG)
GBSG %>% head()
testmfp <- mfp(Surv(rfst, cens) ~ fp(age)+prm+esm+tumsize+menostat+tumgrad+strata(htreat), 
    family = cox, data = GBSG, select = 0.05)

testmfp_summ <- testmfp %>% summary()
testmfp_summ$call

refitmfp <- coxph(formula = Surv(rfst, cens) ~ strata(htreat) + prm + tumsize + 
                    tumgrad + I((age/100)^-1) + I((age/100)^-1 * log((age/100))), 
                  data = GBSG)
visreg(refitmfp, "age", ylab = "log(HR)")

plot(survfit(testmfp$fit, newdata = data.frame(age = c(30, 50, 70), 
                                               prm = rep(40,3),
                                               esm = rep(60,3),
                                               tumsize = rep(20,3),
                                               menostat = rep(2,3),
                                               tumgrad = rep(2,3)
                                               )
             ),
     col = c(2,4)
     )
## repeat for my data
testmfp <- mfp(Surv(t, out) ~ fp(I(exposed*t))+strata(setid), 
    family = cox, data = test2, select = 0.05)

testmfp_summ <- testmfp %>% summary()
testmfp_summ$call 
do.call(testmfp_summ$call)

refitmfp <- coxph(formula = Surv(t, out) ~ 1, data = test2)

visreg(refitmfp, "t", ylab = "log(HR)")
