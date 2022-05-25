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
library(biostat3)
library(gt)
library(lmtest)

if (Sys.info()["user"] == "lsh1510922") {
  if (Sys.info()["sysname"] == "Darwin") {
    datapath <- "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
  if (Sys.info()["sysname"] == "Windows") {
    datapath <- "Z:/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
}
dir.create(file.path(here("out", "analysis")), showWarnings = FALSE)
dir.create(file.path(datapath, "out/models_data"), showWarnings = FALSE)

XX <- c("psoriasis", "eczema")
YY <- c("anxiety", "depression")

full_ethnicity <- NULL
st <- Sys.time()
for (exposure in XX) {
  #exposure <- XX[1]
  ABBRVexp <- str_sub(exposure, 1 , 3)
  if (exposure == "eczema") {
    df_anx_split <-
      readRDS(paste0(datapath, "out/df_model", ABBRVexp, "_anxiety_2006on.rds"))
    df_dep_split <-
      readRDS(paste0(datapath, "out/df_model", ABBRVexp, "_depression_2006on.rds"))
  } else if (exposure == "psoriasis") {
    df_anx_split <-
      readRDS(paste0(datapath, "out/df_model", ABBRVexp, "_anxiety_2006on.rds"))
    df_dep_split <-
      readRDS(paste0(datapath, "out/df_model", ABBRVexp, "_depression_2006on.rds"))
  }
  
  .dib(exposure)
  
  for (outcome in YY) {
    .dib(outcome)
    
    if (outcome == "anxiety") {
      df_model_eth_nonmiss <- df_anx_split
    } else if (outcome == "depression") {
      df_model_eth_nonmiss <- df_dep_split
    }
    
    # Run a simple Cox regression ----------------------------------------
    .dib("Running model 2 (confounder)")
    modEth <-
      coxph(
        Surv(t, out) ~ exposed + eth_edited + carstairs + cal_period + comorbid + cci + strata(setid),
        data = df_model_eth_nonmiss)
    .dib("Running model 3 (mediation)")
    if (ABBRVexp == "ecz") {
      modEth2 <-
        coxph(
          Surv(t, out) ~ exposed + eth_edited + carstairs + cal_period + comorbid + cci + bmi_cat + sleep + alc + smokstatus + gc90days + strata(setid),
          data = df_model_eth_nonmiss
        ) 
    } else if (ABBRVexp == "pso") {
      modEth2 <-
        coxph(
          Surv(t, out) ~ exposed + eth_edited + carstairs + cal_period + comorbid + cci + bmi_cat + alc + smokstatus + strata(setid),
          data = df_model_eth_nonmiss
        ) 
    }
    
    # interaction model ------------------------------------
    .dib("Running model 4 (interaction)")
    modEth3 <-
      coxph(
        Surv(t, out) ~ exposed + agegroup + exposed * eth_edited + eth_edited + carstairs + cal_period + comorbid + cci + strata(setid),
        data = df_model_eth_nonmiss
      ) 
    
    saveRDS(
      modEth,
      file = paste0(
        datapath,
        "out/models_data/",
        ABBRVexp,
        "_",
        outcome,
        "_mod_ethnicity_confound.rds"
      )
    )
    saveRDS(
      modEth2,
      file = paste0(
        datapath,
        "out/models_data/",
        ABBRVexp,
        "_",
        outcome,
        "_mod_ethnicity_mod.rds"
      )
    )
    saveRDS(
      modEth3,
      file = paste0(
        datapath,
        "out/models_data/",
        ABBRVexp,
        "_",
        outcome,
        "_mod_ethnicity_int.rds"
      )
    )
    
    # rerun original models with smaller dataset for lrtest ---------------------------
    df_model_small <- df_model_eth_nonmiss %>% 
      filter(!is.na(eth_edited))
    
    mod1_small <-
      coxph(Surv(t, out) ~ exposed + carstairs + cal_period + comorbid + cci + strata(setid),
            data = df_model_small)
    
    if (ABBRVexp == "ecz") {
      mod2_small <-
        coxph(
          Surv(t, out) ~ exposed + carstairs + cal_period + comorbid + cci + bmi_cat + sleep + alc + smokstatus + gc90days + strata(setid),
          data = df_model_small
        ) 
    } else if (ABBRVexp == "pso") {
      mod2_small <-
        coxph(
          Surv(t, out) ~ exposed + carstairs + cal_period + comorbid + cci + bmi_cat + alc + smokstatus + strata(setid),
          data = df_model_small
        ) 
    }
    lr1 <- lrtest(mod1_small, modEth)
    lr2 <- lrtest(mod2_small, modEth2)
    lr3 <- lrtest(modEth, modEth3)
    
    .dib(paste0(outcome, "~", exposure))
    
    # load models -------------------------------------------------------------
    mod_confound_old <-
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
    mod_mediator_old <-
      readRDS(
        paste0(
          datapath,
          "out/models_data/",
          ABBRVexp,
          "_",
          outcome,
          "_mod3_modeldata.rds"
        )
      )
    mod_confound_new <- modEth
    mod_mediator_new <- modEth2
    
    mod1_old <- broom::tidy(mod_confound_old, conf.int = T, exponentiate = T, conf.level = 0.95)
    mod2_old <- broom::tidy(mod_mediator_old, conf.int = T, exponentiate = T, conf.level = 0.95)
    mod1_new <- broom::tidy(mod_confound_new, conf.int = T, exponentiate = T, conf.level = 0.95)
    mod2_new <- broom::tidy(mod_mediator_new, conf.int = T, exponentiate = T, conf.level = 0.95)
    mod3_new <- broom::tidy(modEth3, conf.int = T, exponentiate = T, conf.level = 0.95)
    
    fmt_tab <- function(modelframe) {
      modelframe %>% 
        filter(str_detect(string = term, pattern = "^exposed.|^eth.")) %>% 
        mutate(Y = paste0(outcome), X = paste0(exposure)) %>% 
        dplyr::select(X, Y, term, estimate, conf.low, conf.high, p.value)
    }
    fmt_mod1_old <- fmt_tab(mod1_old) %>% mutate(model = "mod1_old", n = mod_confound_old$n, n_event = mod_confound_old$nevent, pval = NA)
    fmt_mod2_old <- fmt_tab(mod2_old) %>% mutate(model = "mod2_old", n = mod_mediator_old$n, n_event = mod_mediator_old$nevent, pval = NA)
    fmt_mod1_new <- fmt_tab(mod1_new) %>% mutate(model = "mod1_new", n = mod_confound_new$n, n_event = mod_confound_new$nevent, pval = lr1$`Pr(>Chisq)`[2])
    fmt_mod2_new <- fmt_tab(mod2_new) %>% mutate(model = "mod2_new", n = mod_mediator_new$n, n_event = mod_mediator_new$nevent, pval = lr2$`Pr(>Chisq)`[2])
    fmt_mod3_new <- fmt_tab(mod3_new) %>% mutate(model = "mod3_new", n = modEth3$n, n_event = modEth3$nevent, pval = lr3$`Pr(>Chisq)`[2])
    
    ## add lincom estimates for interaction
    coeffs <- modEth3$coefficients %>% names()
    int_var <- "eth_edited"
    int_levels <- df_model_eth_nonmiss[, int_var] %>% levels()
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
    lincom_out <- lincom(modEth3,
                         interactions,
                         eform = TRUE,
                         level = 0.95)
    rownames(lincom_out)[1] <-
      paste0(rownames(lincom_out)[1], "+", int_var, int_levels[1])
    
    lincom_out <- as_tibble(lincom_out, rownames = "lincom")
    tibble_int <- lincom_out %>%
      janitor::clean_names() %>%
      dplyr::select(lincom,
                    estimate,
                    conf.low = x2_5_percent,
                    conf.high = x97_5_percent,
                    p.value = pr_chisq) %>%
      mutate_at(c("estimate", "conf.low", "conf.high", "p.value"), ~ unlist(.))  %>% 
      dplyr::select(term = lincom, estimate, conf.low, conf.high, p.value) %>% 
      mutate(X = exposure, Y = outcome, model = "interaction",
             n = modEth3$n, n_event = modEth3$nevent, pval = lr3$`Pr(>Chisq)`[2])
    
    full_ethnicity <- full_ethnicity %>% 
      bind_rows(fmt_mod1_old,
                fmt_mod2_old,
                fmt_mod1_new,
                fmt_mod2_new,
                fmt_mod3_new,
                tibble_int)
  }
}
end <- Sys.time()
end-st

saveRDS(full_ethnicity, file = paste0(datapath,  "out/models_data/ethnicity_models.rds"))

# make the plot of interaction --------------------------------------------
tibble_out <- full_ethnicity %>% 
  filter(model == "interaction")
tibble_plot2 <- tibble_out %>% 
  mutate(term2 = term)  %>% 
  mutate_at("term2", ~str_replace_all(., "\\+", ",")) %>% 
  separate(term2, into = c("exposed", "z_level"), sep = ",") %>% 
  mutate_at("z_level", ~str_remove(., "exposedPsoriasis:")) %>% 
  mutate_at("z_level", ~str_remove(., "exposedEczema:")) %>% 
  dplyr::select(-exposed) %>%
  mutate(z_level = as.factor(z_level)) %>% 
  mutate_at(c("X", "Y"), ~stringr::str_to_title(.)) %>% 
  mutate(nice_z = "Ethnicity") %>% 
  mutate(nice_z_level = str_remove(z_level, "eth_edited"))

tibble_plot2$nice_z_level <- fct_reorder(tibble_plot2$nice_z_level, tibble_plot2$estimate)

pd <- position_dodge(width = 0.2)
ggplot(tibble_plot2, aes(x = nice_z_level, y = estimate, ymin = conf.low, ymax = conf.high, group = nice_z_level, colour = nice_z_level)) + 
  geom_point(position = pd, size = 3) +
  geom_errorbar(position = pd, width = 0.25) +
  geom_hline(yintercept = 1, lty=2) +  
  #ylim(c(0,NA)) +
  scale_y_log10(breaks=seq(0,4,0.5), limits = c(0.5, NA)) +
  scale_x_discrete(limits=rev) +
  coord_flip() +
  facet_grid(Y~X, drop = TRUE, space = "free", scales = "free") +
  guides(colour = "none") +
  labs(y = "Hazard ratio", x = "Ethnicity") +
  scale_alpha_identity() +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")
dev.copy(pdf, here::here("out/analysis/forest_plot4_ethnicity.pdf"), width = 8, height = 6); dev.off()


# table of skin HRs with and without ethnicity -----------------------------
ethnicity_table <- full_ethnicity %>% 
  filter(!model %in% c("interaction", "mod3_new"), 
         term %in% c("exposedPsoriasis", "exposedEczema")) %>% 
  mutate(waldP = ifelse(p.value<0.0001, "*", prettyNum(p.value)),
         lrP = ifelse(pval<0.0001, "*", prettyNum(pval, digits = 2))) %>% 
  mutate(estimate = paste0(signif(estimate, 3), " (", signif(conf.low, 3), "-", signif(conf.high, 3), ")")) %>% 
  mutate_if(is.numeric, ~prettyNum(., big.mark = ",")) %>% 
  mutate_at(c("X", "Y"), ~str_to_title(.)) %>% 
  dplyr::select(-p.value, -pval, -term, -conf.low, -conf.high, -waldP) 


ethnicity_table_long <- ethnicity_table %>% 
  pivot_longer(cols = model) %>% 
  mutate(rowname = ifelse(str_detect(value, "^mod1."), "Confounder adjusted", "Mediator adjusted")) %>% 
  mutate(runnum = str_extract(value, "old|new")) %>% 
  mutate(runnum = ifelse(runnum == "old", "Main", "With Ethnicity")) %>% 
  dplyr::select(rowname, n, n_event, everything(), -name, - value) %>% 
  arrange(X, rowname) 


ethnicity_table_long %>%
  mutate(rowname = paste0(rowname, ": ", runnum)) %>% 
  dplyr::select(-runnum) %>% 
  mutate_at("lrP", ~ifelse(is.na(lrP), "-", .)) %>% 
  gt(rowname_col = "rowname",
     groupname_col = c("X", "Y")) %>% 
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups(groups = everything())
  ) %>%
  cols_align(align = "right", 
             columns = c("n", "n_event", "estimate")) %>% 
  cols_label(
    n = md("**N**"),
    n_event = md("**No. events**"),
    estimate = md("**HR (95% CI)**"),
    lrP = md("**p**")
  ) %>% 
  tab_footnote(footnote = md("*p* value from a likelihood ratio test comparing models with and without ethnicity as a covariate"),
               cells_column_labels(columns = lrP)) %>% 
  gt::gtsave(filename = paste0("tab12_ethnicity_mainestimates.html"), 
             path = here::here("out/tables/")) 
