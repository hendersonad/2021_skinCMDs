library(tidyverse)
library(here)
library(labelled)
library(stringr)
library(janitor)
library(timetk)
library(gtsummary)
library(magrittr)
library(survival)
library(flexsurv)
library(survminer)
library(htmltools)
library(lubridate)
library(biostat3)
library(gt)
library(lmtest)
library(flexsurv)

if (Sys.info()["user"] == "lsh1510922") {
  if (Sys.info()["sysname"] == "Darwin") {
    datapath <- "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/"
    datapath <- "/Users/lsh1510922/Documents/Postdoc/2021_extract/"
  }
  if (Sys.info()["sysname"] == "Windows") {
    datapath <- "Z:/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
}
dir.create(file.path(here("out", "analysis")), showWarnings = FALSE)
dir.create(file.path(datapath, "out/models_data"), showWarnings = FALSE)

XX <- c("psoriasis", "eczema")[2]
YY <- c("anxiety", "depression")[2]

exposure <- XX[1]

# regression interaction with time since indexdate ------------------------
tibble_out <- NULL
for(exposure in XX) {
  ABBRVexp <- str_sub(exposure, 1 , 3)
  if (exposure == "eczema") {
    df_anx_split <-
      readRDS(paste0(datapath, "out/ecz-anxiety_split.rds"))
    df_dep_split <-
      readRDS(paste0(datapath, "out/ecz-depression_split.rds"))
  } else if (exposure == "psoriasis") {
    df_anx_split <-
      readRDS(paste0(datapath, "out/pso-anxiety_split.rds"))
    df_dep_split <-
      readRDS(paste0(datapath, "out/pso-depression_split.rds"))
  }
  
  .dib(exposure)

  for(outcome in YY) { 
    outcome = "anxiety"
    
    if (outcome == "anxiety") {
      df_outcome <- df_anx_split
    } else if (outcome == "depression") {
      df_outcome <- df_dep_split
    }
    
    df_outcome$gc90days <-
      factor(df_outcome$gc90days, levels = c("not active", "active"))
    
    df_outcome$bmi2 <-
      df_outcome$bmi - mean(df_outcome$bmi, na.rm = T)
    
    build_split <- function(df_in = anxiety_static){
      tosplit <- df_in %>% 
        arrange(setid, patid)
      
      ## re-index on indexdate not DOB
      df_split_pre <- tosplit %>%
        mutate(dobNum = as.numeric(dob, origin = ("1970-01-01")),
               tstart= dobNum + tstart - as.numeric(indexdate, origin = "1970-01-01"),
               tstop = dobNum + tstop - as.numeric(indexdate, origin = "1970-01-01")
              )
      ## split on age
      df_split_post <- survSplit(Surv(time = as.numeric(tstart), 
                                     time2 = as.numeric(tstop),
                                     event = out) ~ ., 
                                data=df_split_pre, 
                                cut=c(1, 3, 5, 100)*365.25,
                                episode="t_fup")
      #Make factors
      df_split_post$t_fup <- factor(df_split_post$t_fup, levels = 1:4,
                                        labels = c("<1 yr", "1-3 yrs", "3-5 yrs", "5+ yrs"))
      
      df_split_post
    }
    
    df_split_index <- build_split(df_outcome)
    
    df_split_index$t <-
      as.numeric(df_split_index$tstop - df_split_index$tstart)
    
    
    # Bit of variable formatting for output in regression tables --------------
    df_model <- df_split_index %>% 
      mutate(exposed = case_when(
        exposed == 0 ~ "Unexposed",
        exposed == 1 ~ stringr::str_to_title(paste0(exposure))
      ))
    #rename(!!paste0(exposure) := exposed) %>% 
    df_model$exposed <- factor(df_model$exposed, levels = c("Unexposed", str_to_title(exposure)))
    var_label(df_model$exposed) <- "Exposure"
    var_label(df_model$carstairs) <- "Carstairs index of deprivation"
    levels(df_model$carstairs)[1] <- "1 (least deprived)"
    levels(df_model$carstairs)[5] <- "5 (most deprived)"
    var_label(df_model$cal_period) <- "Calendar Period"
    var_label(df_model$cci) <- "Charlson's comorbidity index"
    levels(df_model$cci) <- c("Low", "Moderate", "Severe")
    var_label(df_model$agegroup) <- "Age group"
    df_model$agegroup <- relevel(df_model$agegroup, ref = "50-59")
    
    # Mediators
    var_label(df_model$bmi2) <- "BMI (centred)"
    var_label(df_model$alc) <- "Alcohol misuse"
    levels(df_model$alc) <- c("No", "Yes")
    var_label(df_model$smokstatus) <- "Smoking status"
    levels(df_model$smokstatus) <- stringr::str_to_title(levels(df_model$smokstatus))
    var_label(df_model$gc90days) <- "Recent high dose glucocorticoid steroid use (<30days)"
    levels(df_model$gc90days) <- c("No", "Yes")
    df_model$sleep <- factor(df_model$sleep, levels = 0:1, labels = c("No", "Yes"))
    var_label(df_model$sleep) <- "Sleep problems"
    var_label(df_model$severity) <- paste(str_to_title(exposure), "severity", sep = " ")
    levels(df_model$severity) <- str_to_title(levels(df_model$severity))
    if (ABBRVexp == "ecz") {
      df_model$comorbid <- df_model$asthma
      var_label(df_model$comorbid) <- "Asthma"
      levels(df_model$comorbid) <- c("No", "Yes")
    } else{
      df_model$comorbid <- df_model$arthritis
      var_label(df_model$comorbid) <- "Arthritis"
      levels(df_model$comorbid) <- c("No", "Yes")
    }
    
    .dib("loading confounder adjusted model")
    simple_model <-
      readRDS(
        paste0(
          datapath,
          "out/models_data/",
          ABBRVexp,
          "_",
          outcome,
          "_mod2_modeldata.rds"
        )
      )
    
    .dib("Running model 5 (interaction with time since index)")
    mod5 <-
      coxph(
        Surv(t, out) ~ exposed + t_fup + exposed * t_fup + agegroup + carstairs + cal_period + comorbid + cci + strata(setid),
        data = df_model
      ) 
    
    summary(mod5)
    summary(simple_model)
    
    int_test <- lrtest(simple_model, mod5)
    
    # calculating the effect modification of the hazard ratios ----------------
    interaction_model <- mod5
    
    coeffs <- interaction_model$coefficients %>% names()
    int_var <- "t_fup"
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
    lincom_out <- lincom(interaction_model,
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
      mutate(y = outcome, x = exposure, z = "t_fup", pval = int_test$`Pr(>Chisq)`[2])
    
    tibble_out <- tibble_out %>% 
      bind_rows(tibble_int)
  }
}

tibble_out

tibble_plot2 <- tibble_out %>% 
  mutate(lincom2 = lincom) %>% 
  mutate_at("lincom2", ~str_replace_all(lincom, "\\+", ",")) %>% 
  mutate_at("lincom2", ~str_replace_all(lincom2, "5,", "5plus"))  %>% 
  separate(lincom2, into = c("exposed", "z_level"), sep = ",") %>% 
  mutate_at("z_level", ~str_remove(., "exposedPsoriasis:")) %>% 
  mutate_at("z_level", ~str_remove(., "exposedEczema:")) %>% 
  dplyr::select(-exposed) %>%
  mutate(z_level = as.factor(z_level)) %>% 
  mutate_at(c("x", "y", "z"), ~stringr::str_to_title(.)) %>% 
  mutate(nice_z = "Time since index date") %>% 
  mutate(nice_z_level= str_remove(z_level, "t_fup")) 

pd <- position_dodge(width = 0.2)

ggplot(tibble_plot2, aes(x = nice_z_level, y = estimate, ymin = conf.low, ymax = conf.high, group = z, colour = z)) + 
  geom_point(position = pd, size = 3) +
  geom_errorbar(position = pd, width = 0.25) +
  geom_hline(yintercept = 1, lty=2) +  
  #ylim(c(0,NA)) +
  scale_y_log10(breaks=seq(0,4,0.1), limits = c(0.9, NA)) +
  scale_x_discrete(limits=rev) +
  coord_flip() +
  facet_grid(y~x, drop = TRUE, space = "free", scales = "free") +
  guides(colour = "none") +
  labs(y = "Hazard ratio", x = "Level") +
  scale_alpha_identity() +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

dev.copy(pdf, here::here("out/analysis/forest_plot5_fup.pdf"), width = 8, height = 6); dev.off()


# Royston-Parmer flexible parametric model --------------------------------
test <- df_model

test$age <- (as.numeric(test$indexdate, origin = "1970-01-01") + test$tstart - test$dobNum)/365.25
meanage <- mean(test$age, na.rm = T)
meanage
test$age <- test$age - meanage

hgengammaPH <- function(x, dummy, mu=0, sigma=1, Q){
dummy * hgengamma(x=x, mu=mu, sigma=sigma, Q=Q)
}
HgengammaPH <- function(x, dummy, mu=0, sigma=1, Q){
dummy * Hgengamma(x=x, mu=mu, sigma=sigma, Q=Q)
}
custom.gengammaPH <- list(name="gengammaPH",
                      pars=c("dummy","mu","sigma","Q"), location="dummy",
                      transforms=c(log, identity, log, identity),
                      inv.transforms=c(exp, identity, exp, identity),
                      inits=function(t){
                      lt <- log(t[t>0])
                      c(1, mean(lt), sd(lt), 0)
                      })

modcox <- coxph(
  Surv(t/365, out) ~ 
    exposed + 
    gender + 
    age , data=test
)
modph <- flexsurvreg(
  Surv(t/365, out) ~ 
    exposed + 
    gender + 
    age , data=test, dist=custom.gengammaPH, fixedpars=1)
mod1 <- flexsurvspline(
  Surv(t/365, out) ~ 
    exposed + 
    gamma1(exposed) +
    gamma2(exposed) +
    gamma3(exposed) +
    gender +
    age , data = test, k = 3, scale = "hazard")
summary(test$t/365)
mod1$knots

modcox$coefficients["exposedEczema"]
modph$res["exposedEczema",]
mod1$res["exposedEczema",]

summarydf <- data.frame(exposed=c("Unexposed","Eczema"),
                        gender = rep("Male", 2),
                        age = rep(meanage, 2))
B <- 1000
t <- seq(0.1, 6, by=0.1)
hrspl.est <-
  summary(mod1, t=t, type="hazard",
          newdata=summarydf[2,],ci=FALSE)[[1]][,"est"] /
  summary(mod1, t=t, type="hazard",
          newdata=summarydf[1,],ci=FALSE)[[1]][,"est"]
pars <- normboot.flexsurvreg(mod1, B=B, newdata=summarydf)
hrSPL <- matrix(nrow=B, ncol=length(t))
for (i in seq_along(t)){
  haz.rep.exp <- hsurvspline(t[i], gamma=pars[[2]], knots=mod1$knots)
  haz.rep.un <- hsurvspline(t[i], gamma=pars[[1]], knots=mod1$knots)
  hrSPL[,i] <- haz.rep.exp / haz.rep.un
}
hrSPL <- apply(hrSPL, 2, quantile, c(0.025, 0.975))
hrPH.est <-
  summary(modph, t=t, type="hazard",
          newdata=summarydf[2,],ci=FALSE)[[1]][,"est"] /
  summary(modph, t=t, type="hazard",
          newdata=summarydf[1,],ci=FALSE)[[1]][,"est"]
pars <- normboot.flexsurvreg(modph, B=B, newdata=summarydf)
hrPH <- matrix(nrow=B, ncol=length(t))
for (i in seq_along(t)){
  haz.rep.exp <- do.call(hgengammaPH, c(list(t[i]), as.data.frame(pars[[2]])))
  haz.rep.un <- do.call(hgengammaPH, c(list(t[i]), as.data.frame(pars[[1]])))
  hrPH[,i] <- haz.rep.exp / haz.rep.un
}
hrPH <- apply(hrPH, 2, quantile, c(0.025, 0.975))
plot(t, hrSPL[1,], type="l", ylim = c(0,2), col="red", xlab="years",
     ylab="Hazard ratio (Medium / Good)", lwd=1, lty=2)
lines(t, hrSPL[2,], col="red", lwd=1, lty=2)
lines(t, hrPH[1,], col="darkgray", lwd=1, lty=2)
lines(t, hrPH[2,], col="darkgray", lwd=1, lty=2)
lines(t, hrspl.est, col="red", lwd=2)
lines(t, hrPH.est, col="darkgray", lwd=2)
abline(v = exp(mod1$knots[1]), lty = 2, col = "blue")
abline(v = exp(mod1$knots[2]), lty = 2, col = "blue")
abline(v = exp(mod1$knots[3]), lty = 2, col = "blue")
legend("topright", lwd=c(2,2), col=c("red","darkgray"), bty="n",
       c("Generalized gamma: spline", "Generalized gamma: proportional hazards"))
