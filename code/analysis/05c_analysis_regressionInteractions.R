library(here)
library(tidyverse)
library(gtsummary)
library(gt)
library(survival)
library(lmtest)
library(broom)
library(biostat3)

if (Sys.info()["user"] == "lsh1510922") {
  if (Sys.info()["sysname"] == "Darwin") {
    datapath <- "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
  if (Sys.info()["sysname"] == "Windows") {
    datapath <- "Z:/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
}
dir.create(paste0(datapath, "out/tables/"), showWarnings = FALSE)
dir.create(paste0(datapath, "out/analysis/"), showWarnings = FALSE)
           
YY <- c("depression", "anxiety")
XX <- c("psoriasis", "eczema")

# exposure <- XX[1]
# outcome <- YY[1]
pval_interactions <- function(exposure, outcome) {
  ABBRVexp <- substr(exposure, 1, 3)
  
  df_model <-
    readRDS(paste0(
      datapath,
      "out/df_model",
      ABBRVexp,
      "_",
      outcome, 
      ".rds"
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
  mod8 <-
    readRDS(
      paste0(
        datapath,
        "out/models_data/",
        ABBRVexp,
        "_",
        outcome,
        "_mod8_interaction_carstairs_modeldata.rds"
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
  summary(mod8, conf.int = F)
  
  lr1 <- lrtest(simple_model, mod5)
  lr2 <- lrtest(simple_model, mod6)
  lr3 <- lrtest(simple_model, mod7)
  lr4 <- lrtest(simple_model, mod8)
  
  interactions <- c("Age group", "Comorbidity", "Calendar period", "Carstairs")
  inter_pval <- cbind.data.frame(interactions = interactions,
                                 pval = rbind(
                                   signif(lr1$`Pr(>Chisq)`[2], digits = 2),
                                   signif(lr2$`Pr(>Chisq)`[2], digits = 2),
                                   signif(lr3$`Pr(>Chisq)`[2], digits = 2),
                                   signif(lr4$`Pr(>Chisq)`[2], digits = 2)
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

gt_int_table <- int_table %>% 
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
gt_int_table
gt::gtsave(
  gt_int_table,
  filename =  paste0("tab8_interaction_pvals.html"),
  path = here::here("out/tables/")
)

# calculating the effect modification of the hazard ratios ----------------
tibble_out <- NULL
for(exposure in XX) {
  ABBRVexp <- substr(exposure, 1, 3)
  for (outcome in YY) {
    df_model <-
      readRDS(paste0(
        datapath,
        "out/df_model",
        ABBRVexp, 
        "_",
        outcome, 
        ".rds"
      ))
    
    for (ZZ in c("agegroup", "comorbid", "cal_period", "carstairs")) {
      
      if(ZZ == "agegroup"){
        data_name <- "_mod5_interaction_age_modeldata"
      }
      if(ZZ == "comorbid"){
        data_name <- "_mod6_interaction_comorbid_modeldata"
      }
      if(ZZ == "cal_period"){
        data_name <- "_mod7_interaction_calendar_modeldata"
      }
      if(ZZ == "carstairs"){
        data_name <- "_mod8_interaction_carstairs_modeldata"
      }
      # take a model: pso - dep modified by age group
      interaction_model <-
        readRDS(
          paste0(
            datapath,
            "out/models_data/",
            ABBRVexp,
            "_",
            outcome, 
            data_name, 
            ".rds"
          )
        )

      coeffs <- interaction_model$coefficients %>% names()
      int_var <- ZZ
      int_levels <- df_model[, int_var] %>% pull() %>% levels()
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
        mutate(y = outcome, x = exposure, z = ZZ)
      
      tibble_out <- tibble_out %>% 
        bind_rows(tibble_int)
    }
  }
}

tibble_plot <- tibble_out %>% 
  mutate_at("lincom", ~str_remove(., "Psoriasis")) %>% 
  mutate_at("lincom", ~str_remove(., "Eczema")) %>% 
  mutate_at("lincom", ~str_remove(., "exposed")) %>% 
  mutate_at("lincom", ~str_remove(., "exposedPsoriasis")) %>% 
  mutate_at("lincom", ~str_remove(., "exposedEczema")) %>% 
  mutate_at("lincom", ~str_remove(., "^+.")) 
pd <- position_dodge(width = 0.3)

ybase <- -0.1 + tibble_plot$conf.low %>% min() %>% round(digits = 2) 
yheight <- 0.1 + tibble_plot$conf.high %>% max() %>% round(digits = 2) 

tibble_plot2 <- tibble_out %>% 
  mutate(lincom2 = lincom) %>% 
  mutate_at("lincom2", ~str_replace_all(lincom, "\\+", ",")) %>% 
  mutate_at("lincom2", ~str_replace_all(lincom2, "80,", "80+"))  %>% 
  separate(lincom2, into = c("exposed", "z_level"), sep = ",") %>% 
  mutate_at("z_level", ~str_remove(., "exposedPsoriasis:")) %>% 
  mutate_at("z_level", ~str_remove(., "exposedEczema:")) %>% 
  dplyr::select(-exposed) %>%
  mutate(z_level = as.factor(z_level)) %>% 
  mutate_at(c("x", "y", "z"), ~stringr::str_to_title(.)) %>% 
  mutate(nice_z = ifelse(z == "Carstairs", "Carstairs (deprivation)", 
                         ifelse(z == "Cal_period", "Calendar period", 
                                ifelse(z == "Comorbid" & x == "Psoriasis", "Arthritis", 
                                       ifelse(z == "Comorbid" & x == "Eczema", "Asthma", 
                                              "Age group"))))) %>% 
  mutate(nice_z_level= paste0(nice_z, ": ", str_remove(z_level, str_to_lower(z)))) 

pdf(here::here("out/analysis/forest_plot3_interactions_v2.pdf"), width = 8, height = 6)
p1 <- ggplot(tibble_plot2, aes(x = nice_z_level, y = estimate, ymin = conf.low, ymax = conf.high, group = z, colour = z)) + 
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
print(p1)
dev.off()

