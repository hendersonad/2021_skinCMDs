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
load_data_fn <- function(X, Y, fupmax = Inf){
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
                  carstairs, cci, age, cal_period, out,
                  tstart, tstop, t)

  ## recode variables with `max` value during follow up (bmi, comorbidity, alc, sleep)
  df_exp_tuc <- df_exp_select %>%
    ungroup() %>%
    group_by(patid) %>%
    mutate(rownum = 1:n(), 
           sumfup = cumsum(t)) %>% 
    ungroup() %>% 
    select(rownum,sumfup, patid, comorbid, cci, severity, alc, sleep, gc90days, out) %>%
    filter(rownum == 1 | sumfup <= fupmax) %>% # filter to events only up to fupmax (argument to function)
    mutate_if(is.factor, ~as.integer(ordered(.))) %>% 
    group_by(patid) %>%
    summarise(across(everything(), max)) 
  
  df_exp_tuc$out[df_exp_tuc$rownum == 1 & df_exp_tuc$sumfup < fupmax] <- 0 # suppress out variable = 0 if t > fupmax
  
  ## special for smoking because of weird categories
  df_exp_smok <- df_exp_select %>%
    group_by(patid) %>%
    summarise(smoker = ifelse(any(smokstatus %in% c("Current Smoker", "Ex-Smoker", "Current Or Ex-Smoker")), 1, 0))
  
  ## variables we just want the value at indexdate
  df_exp_index <- df_exp_select %>%
    group_by(patid) %>%
    select(indexdate, enddate, exposed, gender, dob, age, bmi, bmi2, eth_edited, country, ruc, carstairs, cal_period) %>%
    slice(1)
  
  ## need to add duration of disease (in 1-year increments)
  df_exp_fup <- df_exp %>%
    select(setid, patid, tstart, tstop, t) %>%
    group_by(setid, patid) %>%
    mutate(t = tstop[n()] - tstart[1]) %>% 
    mutate(years = t/365.25) %>% 
    slice(1) %>% 
    ungroup()
  
  df_exp_static <- df_exp_index %>%
    left_join(df_exp_tuc, by = c("patid")) %>%
    left_join(df_exp_smok, by = c("patid")) %>% 
    left_join(df_exp_fup, by = c("patid"))
  
  df_exp_static$gender <- factor(df_exp_static$gender, levels = c(NA, "Male", "Female", "Indeterminate", NA))
  df_exp_static$age <- (df_exp_static$tstart)/365.25
  mean_age <- mean(df_exp_static$age, na.rm = T)
  df_exp_static$age <- df_exp_static$age - mean_age
  
  df_exp_static$smoker <- factor(df_exp_static$smoker)
  
  ### remove ordering of factor variables (this was used to select the max(var) per patid but will mess up the regression presentation)
  df_exp_static$cci <- factor(df_exp_static$cci, levels = 1:3, labels = c("Low", "Moderate", "Severe"))
  df_exp_static$comorbid <- factor(df_exp_static$comorbid, levels = 1:2, labels = c("No", "Yes"))
  df_exp_static$alc <- factor(df_exp_static$alc, levels = 1:2, labels = c("No", "Yes"))
  df_exp_static$sleep <- factor(df_exp_static$sleep, levels = 1:2, labels = c("No", "Yes"))
  df_exp_static$gc90days <- factor(df_exp_static$gc90days, levels = 1:2, labels = c("No", "Yes"))
  
  if(X == "eczema"){
    df_exp_static$severity <- factor(df_exp_static$severity, levels = 1:3, labels = c("Mild", "Moderate", "Severe"))
  }else{
    df_exp_static$severity <- factor(df_exp_static$severity, levels = 1:2, labels = c("Mild", "Moderate/severe"))
  }
  
  df_exp_static
}

# 1a - Logistic regression method -----------------------------------------------
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



# 2a - up to 1 year -------------------------------------------------------
pdf(paste0(here("out/predictions/"), "03_multimodel_logisticpredict_1year.pdf"), 10, 10)
par(mfrow = c(4,4), mgp=c(3,1,0))
ii <- 0
for(exposure in XX) {
  ABBRVexp <- substr(exposure, 1, 3)
  for(outcome in YY) {
    ii <- ii+1
    
    df_exp_static <- load_data_fn(X = exposure, Y = outcome, fupmax = 365.25)
    
    # restrict to 1 year follow up 
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
    df_exp_train$out %>% table()
    if(ABBRVexp == "pso"){
      glm(out ~ age + gender + carstairs + cci, data = df_exp_train, family = "binomial")
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
      filename =  paste0("tab1_", ABBRVexp, "_", substr(outcome, 1, 3), "_predictmodel_logistic_1yr.rtf"),
      path = here::here("out//predictions//")
    )
    gt::gtsave(
      predict_gt,
      filename =  paste0("tab1_", ABBRVexp, "_", substr(outcome, 1, 3), "_predictmodel_logistic_1yr.html"),
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


# 2b - up to 3 years -------------------------------------------------------
pdf(paste0(here("out/predictions/"), "03_multimodel_logisticpredict_3year.pdf"), 10, 10)
par(mfrow = c(4,4), mgp=c(3,1,0))
ii <- 0
for(exposure in XX) {
  ABBRVexp <- substr(exposure, 1, 3)
  for(outcome in YY) {
    ii <- ii+1
    
    df_exp_static <- load_data_fn(X = exposure, Y = outcome, fupmax = 365.25*3)
    
    # restrict to 1 year follow up 
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
      filename =  paste0("tab1_", ABBRVexp, "_", substr(outcome, 1, 3), "_predictmodel_logistic_1yr.rtf"),
      path = here::here("out//predictions//")
    )
    gt::gtsave(
      predict_gt,
      filename =  paste0("tab1_", ABBRVexp, "_", substr(outcome, 1, 3), "_predictmodel_logistic_1yr.html"),
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

