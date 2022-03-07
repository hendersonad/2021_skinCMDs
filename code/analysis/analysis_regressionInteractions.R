library(here)
library(gtsummary)
library(gt)
library(survival)
library(lmtest)
library(broom)
library(tidyverse)

if (Sys.info()["user"] == "lsh1510922") {
  if (Sys.info()["sysname"] == "Darwin") {
    #datapath <- "/Users/lsh1510922/Documents/Postdoc/2021_extract/"
    datapath <- "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
  if (Sys.info()["sysname"] == "Windows") {
    datapath <- "Z:/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
}
dir.create(paste0(datapath, "out/supplementary/"))
           
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
  
  # lrtest ------------------------------------------------------------------
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
  summary(simple_model, conf.int = F)
  summary(mod5, conf.int = F)
  summary(mod6, conf.int = F)
  summary(mod7, conf.int = F)
  
  lr1 <- lrtest(simple_model, mod5)
  lr2 <- lrtest(simple_model, mod6)
  lr3 <- lrtest(simple_model, mod7)
  
  interactions <- c("Age group", "Comorbidity", "Calendar period")
  inter_pval <- cbind.data.frame(interactions = interactions,
                                 pval = rbind(
                                   signif(lr1$`Pr(>Chisq)`[2], digits = 2),
                                   signif(lr2$`Pr(>Chisq)`[2], digits = 2),
                                   signif(lr3$`Pr(>Chisq)`[2], digits = 2)
                                 )) %>%
    as_tibble()
  names(inter_pval)[2] <-
    paste0(ABBRVexp, "_", substr(outcome, 1, 3))
  
  inter_pval %>%
    mutate_at(2, ~ prettyNum(.))
}

int1 <- pval_interactions("eczema", "depression")
int2 <- pval_interactions("eczema", "anxiety")
int3 <- pval_interactions("psoriasis", "depression")
int4 <- pval_interactions("psoriasis", "anxiety")

int_table <- int2 %>% 
  left_join(int1, by = "interactions") %>% 
  left_join(int4, by = "interactions") %>% 
  left_join(int3, by = "interactions")

saveRDS(
  int_table,
  file = paste0(
    datapath,
    "out/supplementary/",
    "interaction_pvals.rds"
  )
)

int_table %>% 
  mutate(ecz_anx = as.numeric(ecz_anx),
         ecz_dep = as.numeric(ecz_dep),
         pso_anx = as.numeric(pso_anx),
         pso_dep = as.numeric(pso_dep),
         ) %>% 
  mutate_if(is.numeric, ~ifelse(.<0.0001, paste0("p<0.0001"), paste0(signif(., digits = 1)))) %>% 
  gt() %>% 
  cols_label(
    interactions = "",
    ecz_anx = md("**Anxiety**"),
    ecz_dep = md("**Depression**"),
    pso_anx = md("**Anxiety**"),
    pso_dep = md("**Depression**"),
  ) %>% 
  cols_align(columns = contains("ecz"),
      align = "right") %>% 
  cols_align(columns = contains("pso"),
      align = "right") %>% 
  tab_spanner(label = md("**Eczema**"),
              columns = c(ecz_anx, ecz_dep)) %>% 
  tab_spanner(label = md("**Psoriasis**"),
              columns = c(pso_anx, pso_dep)) 
