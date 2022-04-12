library(here)
library(survival)
library(lmtest)
library(broom)
library(biostat3)
library(multcomp)
library(rms)
library(splines)
library(mgcv)
library(pROC)
library(gt)
library(tidyverse)
library(tidycat)

if (Sys.info()["user"] == "lsh1510922") {
  if (Sys.info()["sysname"] == "Darwin") {
    #datapath <- "/Users/lsh1510922/Documents/Postdoc/2021_extract/"
    datapath <- "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
  if (Sys.info()["sysname"] == "Windows") {
    datapath <- "Z:/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
}

dir.create(file.path(here("out")), showWarnings = FALSE)
dir.create(file.path(here("out", "predictions")), showWarnings = FALSE)

YY <- c("depression", "anxiety")
XX <- c("psoriasis", "eczema")

exposure <- XX[2]
outcome <- YY[2]
#


# 0 - little function to load data and summarise as static df -----------------
load_data_fn <- function(X, Y){
  ABBRVexp <- substr(X, 1, 3)
  
  # load data ---------------------------------------------------------------
  df_model <- readRDS(paste0(datapath, "out/models_data/df_model", ABBRVexp, "_", Y,".rds"))
  
  # restrict to skin disease pop --------------------------------------------
  df_exp <- df_model %>%
    filter(exposed == str_to_title(X))
  
  ## can't have time-updated covariates so collapse
  df_exp_select <- df_exp %>%
    dplyr::select(setid, patid, exposed, indexdate, enddate, dob, gender, comorbid, alc, smokstatus,
                  severity, sleep, sleep_all, gc90days, death, eth_edited, bmi, bmi2, country, ruc,
                  carstairs, cci, age, cal_period, out)
  
  ## recode variables with `max` value during follow up (bmi, comorbidity, alc, sleep)
  df_exp_tuc <- df_exp_select %>%
    select(setid, patid, bmi, bmi2, comorbid, cci, severity, alc, sleep, sleep_all, gc90days, out) %>%
    mutate_if(is.factor, ~ordered(.)) %>%
    group_by(setid, patid) %>%
    summarise(across(everything(), max))
  
  ## special for smoking because of weird categories
  df_exp_smok <- df_exp_select %>%
    group_by(setid, patid) %>%
    summarise(smoker = ifelse(any(smokstatus %in% c("Current Smoker", "Ex-Smoker", "Current Or Ex-Smoker")), 1, 0))
  
  ## variables we just want the value at indexdate
  df_exp_index <- df_exp_select %>%
    group_by(setid, patid) %>%
    select(indexdate, enddate, exposed, gender, dob, age, eth_edited, country, ruc, carstairs, cal_period) %>%
    slice(1)
  
  df_exp_static <- df_exp_index %>%
    left_join(df_exp_tuc, by = c("setid", "patid")) %>%
    left_join(df_exp_smok, by = c("setid", "patid"))
  
  ## need to add duration of disease (in 1-year increments)
  df_exp_fup <- df_exp %>%
    select(setid, patid, tstart, tstop, t) %>%
    group_by(setid, patid) %>%
    mutate(t = tstop[n()] - tstart[1]) %>% 
    mutate(years = t/365.25) %>% 
    slice(1) %>% 
    ungroup()
  
  df_exp_static <- df_exp_static %>%
    left_join(df_exp_fup, by = c("setid", "patid"))
  
  df_exp_static$gender <- factor(df_exp_static$gender, levels = c(NA, "Male", "Female", "Indeterminate", NA))
  df_exp_static$age <- (df_exp_static$tstart)/365.25
  mean_age <- mean(df_exp_static$age, na.rm = T)
  df_exp_static$age <- df_exp_static$age - mean_age
  
  df_exp_static$smoker <- factor(df_exp_static$smoker)
  
  ### remove ordering of factor variables (this was used to select the max(var) per patid but will mess up the regression presentation)
  ordered_vars <- df_exp_static %>% 
    ungroup() %>% 
    select_if(is.ordered) %>% 
    names()
  if(length(ordered_vars) > 0){
    for(varname in ordered_vars) {
      #print(varname)
      x_ordered <- df_exp_static[[varname]]
      x_fact <- factor(x_ordered, ordered = FALSE)
      df_exp_static[[varname]] <- x_fact
    }
  }
  
  df_exp_static
}

# 1- Logistic regression method -----------------------------------------------
for(exposure in XX) {
  ABBRVexp <- substr(exposure, 1, 3)
  for(outcome in YY) {
    df_exp_static <- load_data_fn(X = exposure, Y = outcome)
    
    # sample 80% train --------------------------------------------------------
    length_data <- dim(df_exp_static)[1]

    set.seed(12)
    patid_sample <- sample(df_exp_static$patid, size = round(length_data*0.8))
    
    df_exp_train <- df_exp_static %>%
      ungroup() %>%
      filter(patid %in% patid_sample)

    df_exp_test <- df_exp_static %>%
      ungroup() %>%
      filter(!patid %in% patid_sample)

    # univariable logistic regression with covariates -------------------------
    covars <- c("age", "gender", "carstairs", 
                "cci", "bmi2", "smoker")
    univ_roc <- function(covariate){
      uni_1 <- glm(out ~ get(covariate), data = df_exp_train, family = "binomial")
      
      pred_vals <- predict(uni_1, type = "link", data = df_exp_train)
      
      df_predictions <- df_exp_train %>% 
        filter(!is.na(get(covariate))) %>% 
        mutate(lp = pred_vals) %>% 
        select(patid, all_of(covariate), out, lp) %>% 
        mutate(risk = 1/(1 + exp(-lp)))
      
      # report AUC  -------------------------------------------------------------
      roc_calc <- pROC::roc(df_predictions$out, df_predictions$risk)
      roc_calc$auc
      plot(roc_calc, main = covariate, xlim = c(1,0), ylim = c(0,1))
      text(0.8, 0.8, round(roc_calc$auc,2), font = 2, pos = 4)
      
      df_calibration <- df_predictions %>% 
        ungroup() %>% 
        mutate(risk_dec = ntile(risk, 10)) %>% 
        group_by(risk_dec) %>% 
        summarise(n = n(), observed = mean(out), predicted = mean(risk))
      smoothfit <- loess(df_calibration$observed ~ df_calibration$predicted, degree = 2)
      
      scatter.smooth(df_calibration$predicted, df_calibration$observed, 
                     col = 7, type = "p", xlim = c(0,0.2), ylim = c(0,0.2), 
                     xlab = "Predicted probability", ylab = "Observed probability",
                     lpars = list(col = 4, lty = 2))
      abline(coef = c(0,1))
    }
    
    pdf(paste0(here("out/predictions/"), "01_univ_auc", ABBRVexp, "_", substr(outcome, 1, 3), ".pdf"), 8, 8)
    par(mfrow = c(3,4))
      sapply(covars, FUN = univ_roc)
    dev.off()
  }
}

# 1b - build multivariable logistic regression models --------------------------
pdf(paste0(here("out/predictions/"), "02_multimodel_logisticpredict.pdf"), 10, 10)
par(mfrow = c(4,4), mgp=c(3,1,0))
ii <- 0
for(exposure in XX) {
  ABBRVexp <- substr(exposure, 1, 3)
  for(outcome in YY) {
    ii <- ii+1
    
    df_exp_static <- load_data_fn(X = exposure, Y = outcome)
    
    # sample 80% train 
    length_data <- dim(df_exp_static)[1]
    
    set.seed(12)
    patid_sample <- sample(df_exp_static$patid, size = round(length_data*0.8))
    
    df_exp_train <- df_exp_static %>%
      ungroup() %>%
      filter(patid %in% patid_sample)
    
    df_exp_test <- df_exp_static %>%
      ungroup() %>%
      filter(!patid %in% patid_sample)
    
    if(ABBRVexp == "pso"){
      multi_1 <- glm(out ~ age + gender + carstairs + cci + bmi2 + smoker + alc, data = df_exp_train, family = "binomial")
      model_covars <- c("Intercept","age", "gender", "carstairs", "cci", "bmi2", "smoker", "alc")
      pretty_model_covars <- cbind.data.frame(
        model_covars, 
        pretty = c("Intercept", "Age (centred)", "Gender", "Carstairs index of deprivation", "CCI", "BMI (centred)", "Smoker", "Harmful alcohol use")
      )
    }
    if(ABBRVexp == "ecz"){
      multi_1 <- glm(out ~ age + gender + carstairs + cci + bmi2 + smoker + alc + sleep + gc90days, data = df_exp_train, family = "binomial")
      model_covars <- c("Intercept", "age", "gender", "carstairs", "cci", "bmi2", "smoker", "alc","sleep", "gc90days")
      pretty_model_covars <- cbind.data.frame(
        model_covars, 
        pretty = c("Intercept", "Age (centred)", "Gender", "Carstairs index of deprivation", "CCI", "BMI (centred)", "Smoker", "Harmful alcohol use", "Sleep problems", "Oral GC use (90 day risk window)")
      )
    }
    
    predict_cis <- confint.default(multi_1) %>% 
      as.data.frame(row.names = F) %>% 
      janitor::clean_names()
    predict_gt <- broom::tidy(multi_1, conf.int = F) %>% 
      bind_cols(predict_cis) %>% 
      select(variable = term, logOR = estimate, conf.low = x2_5_percent, conf.high = x97_5_percent, p.value) %>%
      drop_na() %>% 
      mutate(conf_int = paste0(signif(conf.low, 2), " , ", signif(conf.high, 2)),
             p = ifelse(p.value < 0.0001, "*", paste0(signif(p.value, 1)))) %>% 
      select(-conf.low, -conf.high, -p.value) %>% 
      separate(variable, into = c("delete", "level"), paste(model_covars, collapse = "|"), remove = FALSE) %>% 
      mutate(level=str_remove(level, "\\)")) %>% 
      mutate(temp = str_extract(variable, paste(model_covars, collapse = "|"))) %>% 
      left_join(pretty_model_covars, by = c("temp" = "model_covars")) %>% 
      mutate(OR = exp(logOR)) %>% 
      select(var = pretty, level, OR, logOR, conf_int, p) %>% 
      gt() %>% 
      cols_align(columns = 3:6, align = "right") %>% 
      fmt_number(n_sigfig = 3, columns = where(is.numeric)) %>% 
      cols_label(
        var = "Variable",
        level = "Level",
        logOR = "log(OR)",
        conf_int = "95% CI",
        p = md("*p*")
      ) %>% 
      tab_footnote("*  p < 0.0001", locations = cells_column_labels("p"))
    predict_gt
    gt::gtsave(
      predict_gt,
      filename =  paste0("tab1_", ABBRVexp, "_", substr(outcome, 1, 3), "_predictmodel_logistic.rtf"),
      path = here::here("out/predictions/")
    )
    gt::gtsave(
      predict_gt,
      filename =  paste0("tab1_", ABBRVexp, "_", substr(outcome, 1, 3), "_predictmodel_logistic.html"),
      path = here::here("out//predictions//")
    )  
    
    pred_vals_train <- predict(multi_1, type = "link", newdata = df_exp_train)
    pred_vals_test <- predict(multi_1, type = "link", newdata = df_exp_test)
    
    df_predictions_train <- df_exp_train %>% 
     #filter_at(all_of(model_covars[-1]), all_vars(!is.na(.))) %>% 
      mutate(lp = pred_vals_train) %>% 
      select(patid, all_of(model_covars[-1]), out, lp) %>% 
      mutate(risk = 1/(1 + exp(-lp)))
    
    df_predictions_test <- df_exp_test %>% 
      #filter_at(all_of(model_covars[-1]), all_vars(!is.na(.))) %>% 
      mutate(lp = pred_vals_test) %>% 
      select(patid, all_of(model_covars[-1]), out, lp) %>% 
      mutate(risk = 1/(1 + exp(-lp)))
    
    # PLOT PLOTS PLOTS 
    # report AUC  
    plot_roc <- function(test_train) {
      roc_df <- get(paste0("df_predictions_", test_train))
      roc_calc <- pROC::roc(roc_df$out, roc_df$risk)
      roc_calc$auc
      
      par(new = T)
      ii <- ifelse(test_train=="train", 0, 0.2)
      col_plot <- ifelse(test_train == "train", 4, 2)
      lines(roc_calc$specificities, roc_calc$sensitivities, col = col_plot)
      text(0.9, 0.9-ii, round(roc_calc$auc,2), col = col_plot, font = 2, pos = 4)
    }
    plot(1:0, 0:1, xlim = c(1,0), ylim = c(0,1), col = 0,
         ylab = "Sensitivity", xlab = "Specificity", 
         main = paste0(exposure, " ~ ", outcome))
    abline(coef = c(1,-1))
    plot_roc("train")
    plot_roc("test")
    legend("bottomright", legend = c("Train", "Test"), col = c(4,2), lty = 1, bty = "n")
    mtext(paste0(ii,"A"), side=3, adj=0, font=2)
    
    ## plot_boxplot of risk scores 
    risk_nonoutcome_train <- df_predictions_train[df_predictions_train$out == 0, "risk"] %>% pull()
    risk_withoutcome_train <- df_predictions_train[df_predictions_train$out == 1, "risk"] %>% pull()
    
    risk_nonoutcome_test <- df_predictions_test[df_predictions_test$out == 0, "risk"] %>% pull()
    risk_withoutcome_test <- df_predictions_test[df_predictions_test$out == 1, "risk"] %>% pull()
    
    # Make a list of these 2 vectors
    risk_list <- list(
      risk_nonoutcome_train,
      risk_withoutcome_train,
      risk_nonoutcome_test,
      risk_withoutcome_test
      )
    
    # Change the names of the elements of the list :
    names(risk_list) <- c(paste("Train data \n Control \n n=", length(risk_nonoutcome_train), sep = ""),  
                          paste("Train data \n Case \n n=", length(risk_withoutcome_train), sep = ""),  
                          paste("Test data \n Control \n n=", length(risk_nonoutcome_test), sep = ""),  
                          paste("Test data \n Case \n n=", length(risk_withoutcome_test), sep = "")
                          )
    
    # Change the mgp argument: avoid text overlaps axis
    # Final Boxplot
    mu1a <- signif(mean(risk_nonoutcome_train, na.rm = T), digits = 3)
    mu1b <- signif(mean(risk_nonoutcome_test, na.rm = T), digits = 3)
    text1_train <- bquote(mu ~ "=" ~ .(mu1a))
    text1_test <- bquote(mu ~ "=" ~ .(mu1b))
    mu2a <- signif(mean(risk_withoutcome_train, na.rm = T), digits = 2)
    mu2b <- signif(mean(risk_withoutcome_test, na.rm = T), digits = 2)
    text2_train <- bquote(mu ~ "=" ~ .(mu2a))
    text2_test <- bquote(mu ~ "=" ~ .(mu2a))
    
    col1 <- 1
    par(mgp = c(3,2,0), tck = NA, tcl = -0.25)
    boxplot(risk_list , 
            col= ggplot2::alpha(c(2,4,2,4), 0.2),
            ylab="Survival risk", outline = FALSE, ylim = c(0,0.5),
            pars=list(mgp=c(4,2,.5)))
    text(0.75, 0.4, text1_train, pos = 4, cex = 0.7, col =2)
    text(1.75, 0.4, text2_train, pos = 4, cex = 0.7, col = 4)
    text(2.75, 0.4, text1_test, pos = 4, cex = 0.7, col = 2,)
    text(3.75, 0.4, text2_test, pos = 4, cex = 0.7, col = 4)
    mtext(paste0(ii,"B"), side=3, adj=0, font=2)
    par(mgp = c(3,1,0), tck = NA, tcl = -0.5)
    
    #### GAM plot 
    plot_calibration <- function(test_train) {
      df_calibration <- get(paste0("df_predictions_", test_train))

      gam1 <- gam(out ~ s(risk, k=4) , data = df_calibration, family = "binomial")
      
      sample_plot <- sample(1:dim(df_calibration)[1], size = 0.1*dim(df_calibration)[1])
      
      col_plot <- ifelse(test_train =="train", 4, 2)
      plot_adjust <- ifelse(test_train =="train", 0.01, -0.01)
      axismax <- max(df_predictions_train$risk, na.rm = T)
      df_calibration$out_plot <- ifelse(df_calibration$out==1, axismax, 0)
      points(df_calibration$risk[sample_plot], df_calibration$out_plot[sample_plot]+plot_adjust, col = ggplot2::alpha(col_plot,0.025), cex = 0.2)
      
      tt <- seq(range(df_calibration$risk, na.rm = T)[1],range(df_calibration$risk, na.rm = T)[2],0.001)
      preds <- predict(gam1, newdata = list(risk=tt), type = "link", se.fit = TRUE)
      
      critval <- 1.96;
      upperCI <- preds$fit + (critval * preds$se.fit); 
      lowerCI <- preds$fit - (critval * preds$se.fit)
      
      fit <- preds$fit
      fitPlotF <- gam1$family$linkinv(fit); 
      CI1plotF <- gam1$family$linkinv(upperCI);  
      CI2plotF <- gam1$family$linkinv(lowerCI)
      
      ## Plot GAM fits
      polygon(c(tt,rev(tt)),c(CI1plotF,rev(CI2plotF)),col=ggplot2::alpha(col_plot,0.2),lty=0)
      lines(tt, fitPlotF ,col=col_plot,lwd=1)
    }
    axismax <- max(df_predictions_train$risk, na.rm = T)
    plot(c(0,axismax), c(0,axismax), 
         ylim = c(0-0.02, axismax+0.02), xlim = c(0, axismax),
         xlab = "Predicted probability", ylab = "Observed outcome", 
         col = 0)
    abline(coef = c(0,1), col = ggplot2::alpha(1,0.2))
    plot_calibration("train")
    plot_calibration("test")
    legend("left", legend = c("Train", "Test"), col = c(4,2), lty = 2, bty = "n")
    axis(side = 4, at = c(0-0.02, axismax+0.02), labels = c("Control","Case"), tick = FALSE, padj = -1)
    mtext(paste0(ii,"C"), side=3, adj=0, font=2)
    
    plot_validation <- function(test_train) {
      cal_df <- get(paste0("df_predictions_", test_train))
      
      df_calibration <- cal_df %>% 
        ungroup() %>% 
        mutate(risk_dec = ntile(risk, 10)) %>% 
        group_by(risk_dec) %>% 
        summarise(n = n(), observed = mean(out), predicted = mean(risk))
      
      col_plot <- ifelse(test_train == "train", 4, 2)
      xy <- xy.coords(df_calibration$predicted, df_calibration$observed, "Predicted probability", "Observed probability")
      x <- xy$x
      y <- xy$y
      pred <- loess.smooth(x, y, span = 2/3, degree = 2)
      
      if(test_train == "test"){ 
        par(new = T)
      }
      points(x, y, col = col_plot, cex = 1.2)
      lines(pred$x, pred$y, lty = 2, col = col_plot)
    }
    axismax <- max(df_predictions_train$risk, na.rm = T)
    plot(0:axismax, 0:axismax, 
         ylim = c(0, axismax), xlim = c(0, axismax),
         xlab = "Predicted probability", ylab = "Observed probability", 
         col = 0)
    abline(coef = c(0,1))
    plot_validation("train")
    plot_validation("test")
    legend("bottomright", legend = c("Train", "Test"), col = c(4,2), lty = 2, bty = "n")
    mtext(paste0(ii,"D"), side=3, adj=0, font=2)
  }
}
dev.off()

# 
# # 2 - Redo with Cox and survival models to get curves ----------------------
# par(mfrow = c(2,2))
# for(exposure in XX) {
#   ABBRVexp <- substr(exposure, 1, 3)
#   for(outcome in YY) {
#     
#     # load data ---------------------------------------------------------------
#     df_model <- readRDS(paste0(datapath, "out/models_data/df_model", ABBRVexp, "_", outcome,".rds"))
#     
#     # restrict to skin disease pop --------------------------------------------
#     df_exp_surv <- df_model %>% 
#       filter(exposed == str_to_title(exposure))
#     df_exp_surv$gender <- factor(df_exp_surv$gender, levels = c(NA, "Male", "Female", "Indeterminate", NA))
#     df_exp_surv$age <- (df_exp_surv$tstart)/365.25
#     
#     # sample 80% train --------------------------------------------------------
#     length_data <- length(unique(df_exp_surv$patid))
#     
#     set.seed(12)
#     patid_sample <- sample(df_exp_surv$patid, size = round(length_data*0.8))
#     
#     df_exp_train <- df_exp_surv %>% 
#       ungroup() %>% 
#       filter(patid %in% patid_sample)
#     
#     df_exp_test <- df_exp_surv %>% 
#       ungroup() %>% 
#       filter(!patid %in% patid_sample)
#     
#     
#     # univariable logistic regression with covariates -------------------------
#     # covars <- c("age", "gender", "carstairs", 
#     #             "eth_edited", "cci", "ruc",
#     #             "bmi2", "smokstatus", "alc")
#     # univ_roc <- function(covariate){
#     #   uni_1 <- coxph(Surv(t, out) ~ get(covariate), data = df_exp_train)
#     #   
#     #   base_hazards <- basehaz(uni_1)
#     #   for(tt in seq(365, 365*2, 365*5)){
#     #     base_hazard_t <- base_hazards %>% filter(time == tt) %>%  pull(hazard) 
#     #     base_surv_t <- 1 - base_hazard_t
#     #     
#     #     predict_surv <- df_exp_train %>% 
#     #       filter(!is.na(covariate)) %>% 
#     #       select(patid,out, dob,indexdate, t, tstart, tstop, all_of(covariate)) %>% 
#     #       group_by(patid) %>% 
#     #       mutate(fupstart = tstart[1],
#     #              fup = tstop - fupstart) %>% 
#     #       ungroup()
#     #     
#     #     lin_predictor <- predict(uni_1, type = "lp", newdata = predict_surv)
#     #     
#     #     predict_surv <- predict_surv %>% 
#     #       mutate(basesurv = base_surv_t,
#     #              lp = lin_predictor,
#     #              risk = 1-(basesurv^(exp(lin_predictor))))
#     #     ## who had the outcome at time t 
#     #     predict_out <- df_exp_train %>% 
#     #       filter(!is.na(covariate)) %>% 
#     #       mutate(out_t = ifelse(t<=tt & out == 1, 1, 0)) %>% 
#     #       group_by(patid) %>% 
#     #       summarise(out_t = max(out_t))
#     #     
#     #     predict_roc <- predict_surv %>% 
#     #       filter(fup<=tt)
#     #     
#     #     # report AUC  -------------------------------------------------------------
#     #     roc_calc <- pROC::roc(predict_roc$out, predict_roc$risk)
#     #     roc_calc$auc
#     #     
#     #     if(tt==365){
#     #       ii <- 1
#     #       plot(roc_calc, main = str_to_title(covariate))
#     #       text(1, 1, round(roc_calc$auc,2), font = 2, pos = 4)
#     #       legend("bottomright", legend = c("1 year", "3 year", "5 year"), col = 1:3, lty = 1, bty = "n")
#     #     }else{
#     #       lines(roc_calc$specificities, roc_calc$sensitivities, col = ii+1)
#     #       text(1, 1-(0.2*ii), round(roc_calc$auc,2), font = 2, pos = 4, col = ii+1)
#     #       ii <- ii+1
#     #     }
#     #   }
#     # }
#     # 
#     # pdf(paste0(here("out/predictions/"), "03_Cox_univ_auc", ABBRVexp, "_", substr(outcome, 1, 3), ".pdf"), 8, 8)
#     # par(mfcol = c(3,3))
#     # sapply(covars, FUN = univ_roc)
#     # dev.off()
#     
#     # build multivariable logistic regression models --------------------------
#     if(ABBRVexp == "pso"){
#       multi_1 <- coxph(Surv(t, out) ~ age + gender + carstairs + cci + bmi2 + smokstatus + alc, data = df_exp_train)
#       model_covars <- c("age", "gender", "carstairs", "cci", "bmi2", "smokstatus", "alc")
#       pretty_model_covars <- cbind.data.frame(
#         model_covars, 
#         pretty = c("Age", "Gender", "Carstairs index of deprivation", "CCI", "BMI (centred)", "Smoker", "Harmful alcohol use")
#       )
#     }
#     if(ABBRVexp == "ecz"){
#       multi_1 <- coxph(Surv(t, out) ~ age + gender + carstairs + cci + bmi2 + smokstatus + alc + sleep + gc90days, data = df_exp_train)
#       model_covars <- c("age", "gender", "carstairs", "cci", "bmi2", "smokstatus", "alc","sleep", "gc90days")
#       pretty_model_covars <- cbind.data.frame(
#         model_covars, 
#         pretty = c("Age", "Gender", "Carstairs index of deprivation", "CCI", "BMI (centred)", "Smoker", "Harmful alcohol use", "Sleep problems", "Oral GC use (90 day risk window)")
#       )
#     }
#     predict_gt <- broom::tidy(multi_1, conf.int = T, conf.level = 0.95) %>% 
#       select(variable = term, logHR = estimate, conf.low, conf.high, p.value) %>% 
#       drop_na() %>% 
#       mutate(conf_int = paste0(signif(conf.low, 2), " , ", signif(conf.high, 2)),
#              p = ifelse(p.value < 0.0001, "*", paste0(signif(p.value, 1)))) %>% 
#       select(-conf.low, -conf.high, -p.value) %>% 
#       separate(variable, into = c("delete", "level"), paste(model_covars, collapse = "|"), remove = FALSE) %>% 
#       mutate(temp = str_extract(variable, paste(model_covars, collapse = "|"))) %>% 
#       left_join(pretty_model_covars, by = c("temp" = "model_covars")) %>% 
#       mutate(HR = exp(logHR)) %>% 
#       select(var = pretty, level, HR, logHR, conf_int, p) %>% 
#       gt() %>% 
#       cols_align(columns = 3:6, align = "right") %>% 
#       fmt_number(n_sigfig = 3, columns = where(is.numeric)) %>% 
#       cols_label(
#         var = "Variable",
#         level = "Level",
#         logHR = "log(HR)",
#         conf_int = "95% CI",
#         p = md("*p*")
#       ) %>% 
#       tab_footnote("*  p < 0.0001", locations = cells_column_labels("p"))
#     predict_gt
#     gt::gtsave(
#       predict_gt,
#       filename =  paste0("tab1_", ABBRVexp, "_", substr(outcome, 1, 3), "_predictmodel.rtf"),
#       path = here::here("out/predictions/")
#     )
#     gt::gtsave(
#       predict_gt,
#       filename =  paste0("tab1_", ABBRVexp, "_", substr(outcome, 1, 3), "_predictmodel.html"),
#       path = here::here("out/predictions/")
#     )  
#     
#     create_prediction_df <- function(df){
#       df %>% 
#         filter_at(all_of(model_covars), all_vars(!is.na(.))) %>% 
#         select(patid,out, dob,indexdate, t, tstart, tstop, all_of(model_covars)) %>% 
#         group_by(patid) %>% 
#         mutate(fupstart = tstart[1],
#                fup = (tstart - fupstart)+1) %>% 
#         ungroup()
#     }
#     
#     predict_surv_train <- create_prediction_df(df_exp_train)
#     predict_surv_test <- create_prediction_df(df_exp_test)
#     
#     predict_df <- data.frame(
#       age = 0,
#       bmi2 = 0,
#       gender = levels(df_exp_train$gender)[1],
#       carstairs = levels(df_exp_train$carstairs)[1],
#       cci = levels(df_exp_train$cci)[1],
#       smokstatus = levels(df_exp_train$smokstatus)[1],
#       alc = levels(df_exp_train$alc)[1],
#       sleep = levels(df_exp_train$sleep)[1],
#       gc90days = levels(df_exp_train$gc90days)[1]
#     )
#     
#     base_surv_train <- survfit(multi_1, newdata = predict_df)
#     temp <- cbind.data.frame(t = base_surv_train$time, surv = base_surv_train$surv)
#     
#     lin_predictor_train <- predict(multi_1, type = "lp", newdata = predict_surv_train)
#     lin_predictor_test <- predict(multi_1, type = "lp", newdata = predict_surv_test)
#     
#     risk_predict <- function(test_train){
#       predict_surv <- get(paste0("predict_surv_", test_train)) %>%
#         left_join(temp, by = "t") %>% 
#         mutate(lp = get(paste0("lin_predictor_", test_train)),
#                risk = 1-(surv^(exp(lp))),
#                surv = (surv^(exp(lp)))) %>% 
#         group_by(patid) %>% 
#         mutate(cumrisk = cumprod(risk),
#                cumsurv = cumprod(surv)) %>% 
#         ungroup()
#     
#       # predict_surv <- predict_surv %>% 
#       #   group_by(patid) %>% 
#       #   slice(n()) %>% 
#       #   mutate(cumsurv_1 = 1-cumsurv)
#         
#       predict_surv
#     }
#     
#     
#     risk_predictions_train <- risk_predict("train")
#     risk_predictions_test <- risk_predict("test")
#     # risk_predictions_train %>% filter(patid == 1081) %>% View()
#     # risk_predictions_train %>% filter(patid == 1081) %>% select(out, t, age, gender, alc, carstairs, sleep, surv, lp, risk, surv, cumrisk, cumsurv) %>% View()
#     # risk_predictions_train %>% filter(patid == 4706) %>% select(out, t, age, gender, alc, carstairs, sleep, surv, lp, risk, surv, cumrisk, cumsurv) %>% View()
#     # risk_predictions_train %>% filter(patid == 89052) %>% select(out, t, age, gender, alc, carstairs, sleep, surv, lp, risk, surv, cumrisk, cumsurv) %>% View()
#     
#     plot_prediction <- function(test_train){
#       df_plot <- get(paste0("risk_predictions_", test_train))
#       df_plot <- df_plot %>% 
#         ungroup() %>% 
#         mutate(pred_dec = ntile(risk, 10))
#       df_plot$pred_dec %>% table()
#       
#       ggplot(df_plot, aes(y = risk, x = as.factor(out))) +
#         geom_boxplot()
#       dev.copy(pdf, paste0(here("out/predictions/"), "06_predictions", ABBRVexp,"_", substr(outcome,1,3),".pdf"), 8, 8)
#         dev.off()
#       sample_plot <- sample(1:dim(df_plot)[1], size = 0.1*dim(df_plot)[1])
#       
#       gam1 <- gam(out ~ s(risk, k=4) , data = df_plot, family = "binomial")
#       if(test_train == "train"){
#         plot(df_plot$risk[sample_plot], df_plot$out[sample_plot]-0.02, col = ggplot2::alpha(4,0.025), cex = 0.2, ylab = str_to_title(outcome), xlab = "Estimated survival risk", ylim = c(-0.02, 1.02))
#       }else{
#         points(df_plot$risk[sample_plot], df_plot$out[sample_plot]+0.02, col = ggplot2::alpha(2,0.025), cex = 0.2)
#       }
#       tt <- seq(0,1,0.01)
#       preds <- predict(gam1, newdata = list(risk=tt), type = "link", se.fit = TRUE)
#       
#       critval <- 1.96;
#       upperCI <- preds$fit + (critval * preds$se.fit); 
#       lowerCI <- preds$fit - (critval * preds$se.fit)
#       
#       fit <- preds$fit
#       fitPlotF <- gam1$family$linkinv(fit); 
#       CI1plotF <- gam1$family$linkinv(upperCI);  
#       CI2plotF <- gam1$family$linkinv(lowerCI)
#       
#       ## Plot GAM fits
#       if(test_train == "train"){
#         polygon(c(tt,rev(tt)),c(CI1plotF,rev(CI2plotF)),col=ggplot2::alpha(4,0.2),lty=0)
#         lines(tt, fitPlotF ,col=4,lwd=1)
#       }else{
#         polygon(c(tt,rev(tt)),c(CI1plotF,rev(CI2plotF)),col=ggplot2::alpha(2,0.2),lty=0)
#         lines(tt, fitPlotF ,col=2,lwd=1)
#       }
#     }
#     
#     plot_prediction("train")
#       mtext(side = 3, paste0(exposure, " ~ ", outcome))
#     plot_prediction("test")  
#   }
# }
# dev.copy(pdf, paste0(here("out/predictions/"), "05_predictions.pdf"), 8, 8)
#     dev.off()
# 
#     
# 
#     
#     
# # 3 - Cox model but remove time-updating vars and make constant ---------------------------------
# ## identify subgroups whose risk increases by >50% 
# pdf(paste0(here("out/predictions/"), "05_predictions_v2.pdf"), 10, 10)
# par(mfrow = c(4,3), mgp=c(3,2,0))
# for(exposure in XX) {
#   ABBRVexp <- substr(exposure, 1, 3)
#   for(outcome in YY) {
#     
#     df_exp_surv <- load_data_fn(X = exposure, Y = outcome)
#     
#     # sample 80% train --------------------------------------------------------
#     length_data <- length(unique(df_exp_surv$patid))
#     
#     set.seed(12)
#     patid_sample <- sample(df_exp_surv$patid, size = round(length_data*0.8))
#     
#     df_exp_train <- df_exp_surv %>% 
#       ungroup() %>% 
#       filter(patid %in% patid_sample)
#     
#     df_exp_test <- df_exp_surv %>% 
#       ungroup() %>% 
#       filter(!patid %in% patid_sample)
#     
#     # build multivariable logistic regression models --------------------------
#     if(ABBRVexp == "pso"){
#       multi_1 <- coxph(Surv(t, out) ~ age + gender + carstairs + cci + bmi2 + smoker + alc, data = df_exp_train)
#       model_covars <- c("age", "gender", "carstairs", "cci", "bmi2", "smoker", "alc")
#       pretty_model_covars <- cbind.data.frame(
#         model_covars, 
#         pretty = c("Age (centred)", "Gender", "Carstairs index of deprivation", "CCI", "BMI (centred)", "Smoker", "Harmful alcohol use")
#       )
#     }
#     if(ABBRVexp == "ecz"){
#       multi_1 <- coxph(Surv(t, out) ~ age + gender + carstairs + cci + bmi2 + smoker + alc + sleep + gc90days, data = df_exp_train)
#       model_covars <- c("age", "gender", "carstairs", "cci", "bmi2", "smoker", "alc","sleep", "gc90days")
#       pretty_model_covars <- cbind.data.frame(
#         model_covars, 
#         pretty = c("Age (centred)", "Gender", "Carstairs index of deprivation", "CCI", "BMI (centred)", "Smoker", "Harmful alcohol use", "Sleep problems", "Oral GC use (90 day risk window)")
#       )
#     }
#     
#     predict_gt <- broom::tidy(multi_1, conf.int = T, conf.level = 0.95) %>% 
#       select(variable = term, logHR = estimate, conf.low, conf.high, p.value) %>% 
#       drop_na() %>% 
#       mutate(conf_int = paste0(signif(conf.low, 2), " , ", signif(conf.high, 2)),
#              p = ifelse(p.value < 0.0001, "*", paste0(signif(p.value, 1)))) %>% 
#       select(-conf.low, -conf.high, -p.value) %>% 
#       separate(variable, into = c("delete", "level"), paste(model_covars, collapse = "|"), remove = FALSE) %>% 
#       mutate(temp = str_extract(variable, paste(model_covars, collapse = "|"))) %>% 
#       left_join(pretty_model_covars, by = c("temp" = "model_covars")) %>% 
#       mutate(HR = exp(logHR)) %>% 
#       select(var = pretty, level, HR, logHR, conf_int, p) %>% 
#       gt() %>% 
#       cols_align(columns = 3:6, align = "right") %>% 
#       fmt_number(n_sigfig = 3, columns = where(is.numeric)) %>% 
#       cols_label(
#         var = "Variable",
#         level = "Level",
#         logHR = "log(HR)",
#         conf_int = "95% CI",
#         p = md("*p*")
#       ) %>% 
#       tab_footnote("*  p < 0.0001", locations = cells_column_labels("p"))
#     predict_gt
#     gt::gtsave(
#       predict_gt,
#       filename =  paste0("tab1_", ABBRVexp, "_", substr(outcome, 1, 3), "_predictmodel_v2.rtf"),
#       path = here::here("out/predictions/")
#     )
#     gt::gtsave(
#       predict_gt,
#       filename =  paste0("tab1_", ABBRVexp, "_", substr(outcome, 1, 3), "_predictmodel_v2.html"),
#       path = here::here("out/predictions/")
#     )  
#     
#     predict_surv_train <- df_exp_train
#     predict_surv_test <- df_exp_test
#     
#     predict_df <- data.frame(
#       age = 0,
#       bmi2 = 0,
#       gender = levels(df_exp_train$gender)[1],
#       carstairs = levels(df_exp_train$carstairs)[1],
#       cci = levels(df_exp_train$cci)[1],
#       smoker = levels(df_exp_train$smoker)[1],
#       alc = levels(df_exp_train$alc)[1],
#       sleep = levels(df_exp_train$sleep)[1],
#       gc90days = levels(df_exp_train$gc90days)[1]
#     )
#     
#     base_surv_train <- survfit(multi_1, newdata = predict_df)
#     surv_probs_baseline <- cbind.data.frame(t = base_surv_train$time, surv = base_surv_train$surv)
#     
#     lin_predictor_train <- predict(multi_1, type = "lp", newdata = predict_surv_train)
#     lin_predictor_test <- predict(multi_1, type = "lp", newdata = predict_surv_test)
#     
#     risk_predict <- function(test_train){
#       predict_surv <- get(paste0("predict_surv_", test_train)) %>%
#         left_join(surv_probs_baseline, by = "t") %>% 
#         mutate(lp = get(paste0("lin_predictor_", test_train)),
#                risk = 1-(surv^(exp(lp))),
#                surv = (surv^(exp(lp)))) 
#     
#       # predict_surv <- predict_surv %>% 
#       #   group_by(patid) %>% 
#       #   slice(n()) %>% 
#       #   mutate(cumsurv_1 = 1-cumsurv)
#         
#       predict_surv
#     }
#     
#     risk_predictions_train <- risk_predict("train")
#     risk_predictions_test <- risk_predict("test")
#     # risk_predictions_train %>% filter(patid == 1081) %>% View()
#     # risk_predictions_train %>% filter(patid == 1081) %>% select(out, t, age, gender, alc, carstairs, sleep, surv, lp, risk, surv, cumrisk, cumsurv) %>% View()
#     # risk_predictions_train %>% filter(patid == 4706) %>% select(out, t, age, gender, alc, carstairs, sleep, surv, lp, risk, surv, cumrisk, cumsurv) %>% View()
#     # risk_predictions_train %>% filter(patid == 89052) %>% select(out, t, age, gender, alc, carstairs, sleep, surv, lp, risk, surv) %>% View()
#     
#     plot_prediction <- function(test_train){
#       df_plot <- get(paste0("risk_predictions_", test_train))
#       
#       df_plot <- df_plot %>% 
#         ungroup() %>% 
#         mutate(pred_dec = ntile(risk, 10))
#       #df_plot$pred_dec %>% table()
#       
#       #ggplot(df_plot, aes(y = risk, x = as.factor(out))) +
#       #  geom_boxplot()
#       #dev.copy(pdf, paste0(here("out/predictions/"), "06_predictions", ABBRVexp,"_", substr(outcome,1,3),".pdf"), 8, 8)
#       #  dev.off()
#       
#       gam1 <- gam(out ~ s(risk, k=4) , data = df_plot, family = "binomial")
#       
#       sample_plot <- sample(1:dim(df_plot)[1], size = 0.1*dim(df_plot)[1])
#       if(test_train == "train"){
#         plot(df_plot$risk[sample_plot], df_plot$out[sample_plot]-0.02, col = ggplot2::alpha(4,0.025), cex = 0.2, ylab = str_to_title(outcome), xlab = "Estimated survival risk", ylim = c(-0.02, 1.02))
#       }else{
#         points(df_plot$risk[sample_plot], df_plot$out[sample_plot]+0.02, col = ggplot2::alpha(2,0.025), cex = 0.2)
#       }
#       tt <- seq(0,1,0.01)
#       preds <- predict(gam1, newdata = list(risk=tt), type = "link", se.fit = TRUE)
#       
#       critval <- 1.96;
#       upperCI <- preds$fit + (critval * preds$se.fit); 
#       lowerCI <- preds$fit - (critval * preds$se.fit)
#       
#       fit <- preds$fit
#       fitPlotF <- gam1$family$linkinv(fit); 
#       CI1plotF <- gam1$family$linkinv(upperCI);  
#       CI2plotF <- gam1$family$linkinv(lowerCI)
#       
#       ## Plot GAM fits
#       if(test_train == "train"){
#         polygon(c(tt,rev(tt)),c(CI1plotF,rev(CI2plotF)),col=ggplot2::alpha(4,0.2),lty=0)
#         lines(tt, fitPlotF ,col=4,lwd=1)
#         legend("right", legend = c("Train", "Test"), col = c(4,2), lty = 1, bty = "n")
#       }else{
#         polygon(c(tt,rev(tt)),c(CI1plotF,rev(CI2plotF)),col=ggplot2::alpha(2,0.2),lty=0)
#         lines(tt, fitPlotF ,col=2,lwd=1)
#       }
#       
#     }
#     plot_prediction("train")
#       mtext(side = 3, paste0(exposure, " ~ ", outcome))
#     plot_prediction("test")  
#     
#     plot_roc <- function(test_train){
#       # plot ROC as well --------------------------------------------------------
#       tt = 365*5
#       base_surv_t <- surv_probs_baseline %>% filter(t == tt) %>%  pull(surv) 
#       
#       predict_surv <- get(paste0("predict_surv_", test_train)) %>% 
#         mutate(basesurv = base_surv_t,
#                lp = get(paste0("lin_predictor_", test_train)),
#                risk = 1-(basesurv^(exp(lp))))  
#       
#       ## who had the outcome at time t 
#       predict_roc <- predict_surv %>% 
#         filter(t<=tt)
#       
#       # report AUC  -------------------------------------------------------------
#       roc_calc <- pROC::roc(predict_roc$out, predict_roc$risk)
#       roc_calc$auc
#       
#       if(test_train=="train"){
#         plot(roc_calc, xlim = c(1,0), col = 4, main = "5 year survival")
#         text(1, 1, round(roc_calc$auc,2), font = 2, pos = 4, col = 4)
#         legend("bottomright", legend = c("Train", "Test"), col = c(4,2), lty = 1, bty = "n")
#       }else{
#         ii <- 1
#         lines(roc_calc$specificities, roc_calc$sensitivities, col = ii+1)
#         text(1, 1-(0.2*ii), round(roc_calc$auc,2), font = 2, pos = 4, col = ii+1)
#       }
#     }
#     plot_roc("train")
#     plot_roc("test")
#     
#     ## plot_boxplot of risk scores 
#       df <- risk_predictions_train
#       
#       # Create 2 vectors
#       risk_nonoutcome <- df[df$out == 0, "risk"] %>% pull()
#       risk_withoutcome <- df[df$out == 1, "risk"] %>% pull()
#       
#       # Make a list of these 2 vectors
#       risk_list <- list(risk_nonoutcome, risk_withoutcome)
#       
#       # Change the names of the elements of the list :
#       names(risk_list) <- c(paste("Control \n n=" , length(risk_nonoutcome) , sep=""), 
#                     paste("Case \n n=" , length(risk_withoutcome),  sep=""))
#       
#       # Change the mgp argument: avoid text overlaps axis
#       # Final Boxplot
#       mu1 <- signif(mean(risk_nonoutcome, na.rm = T), digits = 3)
#       text1 <- bquote(mu ~ "=" ~ .(mu1))
#       mu2 <- signif(mean(risk_withoutcome, na.rm = T), digits = 2)
#       text2 <- bquote(mu ~ "=" ~ .(mu2))
#       
#       col1 <- 4
#       boxplot(risk_list , col= ggplot2::alpha(col1, 0.2), ylab="Survival risk", outline = FALSE, ylim = c(0,1))
#       text(0.75, 0.6, text1, pos = 4, cex = 0.7)
#       text(1.75, 0.6, text2, pos = 4, cex = 0.7)
#   }
# }
#     dev.off()
