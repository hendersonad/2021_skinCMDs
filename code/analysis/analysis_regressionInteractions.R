library(gtsummary)
library(survival)
library(lmtest)

exposure <- "psoriasis"
ABBRVexp <- substr(exposure, 1, 3)
outcome <- "depression"

YY <- c("depression", "anxiety")
XX <- c("psoriasis", "eczema")

pval_interactions <- function(exposure, outcome) {
  ABBRVexp <- substr(exposure, 1, 3)
  
  df_model <-
    readRDS(paste0(
      datapath,
      "out/models_data/df_model",
      ABBRVexp,
      "_anxiety.rds"
    ))
  
  mod5 <-
    readRDS(
      paste0(
        datapath,
        "out/models_data/",
        ABBRVexp,
        "_",
        outcome,
        "_mod5_interaction_age_modeldata.rds"
      )
    )
  
  mod6 <-
    readRDS(
      paste0(
        datapath,
        "out/models_data/",
        ABBRVexp,
        "_",
        outcome,
        "_mod6_interaction_comorbid_modeldata.rds"
      )
    )
  mod7 <-
    readRDS(
      paste0(
        datapath,
        "out/models_data/",
        ABBRVexp,
        "_",
        outcome,
        "_mod7_interaction_calendar_modeldata.rds"
      )
    )
  # tbl_regression(mod5, exp = T)
  #
  # mod5_tab <- mod5 %>%
  #   tbl_regression(
  #     exp = T,
  #     include = c('exposed', 'exposed:agegroup', 'agegroup'),
  #     conf.level = 0.99
  #   )
  
  # lrtest ------------------------------------------------------------------
  mod2 <-
    coxph(Surv(t, out) ~ exposed + carstairs + cal_period + comorbid + cci + strata(setid),
          data = df_model)
  
  mod_simple <-
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
  summary(mod_simple, conf.int = F)
  summary(mod5, conf.int = F)
  summary(mod6, conf.int = F)
  summary(mod7, conf.int = F)
  browser()
  lr1 <- lrtest(mod_simple, mod5)
  lr2 <- lrtest(mod_simple, mod6)
  lr3 <- lrtest(mod_simple, mod7)
  
  interactions <- c("Age group", "Comorbidity", "Calendar period")
  inter_pval <- cbind.data.frame(interactions = interactions,
                                 pval = rbind(
                                   signif(lr1$`Pr(>Chisq)`[2], digits = 2),
                                   signif(lr2$`Pr(>Chisq)`[2], digits = 2),
                                   signif(lr3$`Pr(>Chisq)`[2], digits = 2)
                                 )) %>%
    as.tibble()
  names(inter_pval)[2] <-
    paste0(ABBRVexp, "_", substr(outcome, 1, 3))
  
  inter_pval %>%
    mutate_at(2, ~ prettyNum(.))
}

int1 <- pval_interactions("eczema", "depression")
int2 <- pval_interactions("eczema", "anxiety")
int3 <- pval_interactions("psoriasis", "depression")
int4 <- pval_interactions("psoriasis", "anxiety")

int2 %>% 
  left_join(int1, by = "interactions") %>% 
  left_join(int4, by = "interactions") %>% 
  left_join(int3, by = "interactions") %>% 
  gt() %>% 
  cols_label(
    interactions = "",
    ecz_anx = md("**Anxiety**"),
    ecz_dep = md("**Depression**"),
    pso_anx = md("**Anxiety**"),
    pso_dep = md("**Depression**"),
  ) %>% 
  tab_spanner(label = md("**Eczema**"),
              columns = c(ecz_anx, ecz_dep)) %>% 
  tab_spanner(label = md("**Psoriasis**"),
              columns = c(pso_anx, pso_dep)) 
mod5 %>% summary()
library(broom)

ggforest2(mod_simple, data = df_model)
ggforest2(mod5, data = df_model)
