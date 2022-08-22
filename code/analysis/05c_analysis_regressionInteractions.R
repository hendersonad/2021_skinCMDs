library(here)
library(tidyverse)
library(ggprism)
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
        "_mod6_interaction_sex_modeldata.rds"
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
  lr1 <- lrtest(simple_model, mod5)
  lr2 <- lrtest(simple_model, mod6)
  
  interactions <- c("Age group", "Sex")
  inter_pval <- cbind.data.frame(interactions = interactions,
                                 pval = rbind(
                                   signif(lr1$`Pr(>Chisq)`[2], digits = 2),
                                   signif(lr2$`Pr(>Chisq)`[2], digits = 2)
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
    # bodge of gender levels 
    levels(df_model$gender) <- c(NA, "Male", "Female", NA, NA)
    df_model$gender <- factor(df_model$gender, levels = c("Female", "Male")) ## set Female as baseline
    
    for (ZZ in c("agegroup", "gender")) {
      
      if(ZZ == "agegroup"){
        data_name <- "_mod5_interaction_age_modeldata"
      }
      if(ZZ == "gender"){
        data_name <- "_mod6_interaction_sex_modeldata"
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
      
      coeffs <- tidy(interaction_model) %>% drop_na() %>% pull(term)
      int_var <- ZZ
      int_levels_tbl <- table(df_model[, int_var]) %>% as.data.frame()
      int_levels <- int_levels_tbl %>% filter(Freq>0) %>% pull(ZZ) %>% as.character()
      n_int_levels <- length(int_levels) - 1 # -1 because of reference category
      coeffs_interaction <- coeffs[str_detect(coeffs, paste0(":", int_var))]
      
      interactions <- c(coeffs[1])
      for (ii in 1:n_int_levels) {
        X <-
          paste0(c(coeffs[1], coeffs_interaction[ii]),
                 collapse = "+")
        interactions <- c(interactions, X)
      }
      interactions
      lincom_out <- lincom(interaction_model,
                           interactions,
                           eform = TRUE,
                           level = 0.95,
                           singular.ok = TRUE) ## may be collinearity in the model for sex because the "gender" coefficient is NA (because of the matched sets)
      rownames(lincom_out)[1] <-
        paste0(rownames(lincom_out)[1], "+", int_var, int_levels[1])
      
      lincom_out <- as_tibble(lincom_out, rownames = "lincom")
      tibble_int <- lincom_out %>%
        janitor::clean_names() %>%
        dplyr::select(lincom,
                      estimate,
                      conf.low = x2_5_percent,
                      conf.high = x97_5_percent) %>%
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
  mutate(nice_z_level = paste0(ifelse(z == "Gender", "Sex", z), ": ", str_remove(z_level, str_to_lower(z)))) 

## merge on p-values
int_table_merge <- int_table %>% 
  pivot_longer(cols = ecz_anx:pso_dep) %>% 
  separate(name, into = c("exposure", "outcome"), sep = "_") %>% 
  mutate(x = ifelse(exposure == "ecz", "Eczema", "Psoriasis")) %>% 
  mutate(y = ifelse(outcome == "anx", "Anxiety", "Depression")) %>% 
  mutate(z = ifelse(interactions == "Age group", "Agegroup", "Gender")) %>% 
  mutate(p = as.numeric(value)) %>% 
  mutate(nice_p = ifelse(p<0.001, "*", paste0("p=", signif(p, digits = 1))))
tibble_plot3 <- tibble_plot2 %>% 
  left_join(int_table_merge, by = c("z", "x", "y"))
write.csv(tibble_plot3, here::here("out/supplementary/interaction_HRs.csv"))

tibble_plot3$text_hr <- str_pad(round(tibble_plot3$estimate,2), 4, pad = "0", side = "right")
tibble_plot3$text_ciL <- str_pad(round(tibble_plot3$conf.low,2), 4, pad = "0", side = "right")
tibble_plot3$text_ciU <- str_pad(round(tibble_plot3$conf.high,2), 4, pad = "0", side = "right")

tibble_plot3$text_to_plot <- paste0(tibble_plot3$text_hr, 
                                    " [", 
                                    tibble_plot3$text_ciL, 
                                    ",", 
                                    tibble_plot3$text_ciU,
                                    "]")
# add pvalue
tibble_plot3$text_to_plot <- str_pad(paste(tibble_plot3$text_to_plot, tibble_plot3$nice_p), 23, pad = " ", side = "right")

pdf(here::here("out/analysis/forest_plot3_interactions_v2.pdf"), width = 8, height = 6)
pd <- position_dodge(width = 0.3)
p1 <- ggplot(tibble_plot3, aes(x = nice_z_level, group = z, colour = z, label = text_to_plot)) + 
  geom_point(aes(y = estimate), position = pd, size = 3, shape = 1) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = pd, width = 0.25) +
  geom_hline(yintercept = 1, lty=2, col = alpha(1,0.4)) +  
  geom_text(aes(y = 0.999), 
            size = 3.6,
            colour = 1,
            show.legend = FALSE, 
            hjust = 1) +
  scale_y_log10(breaks=seq(0,4,0.1), limits = c(0.9, NA)) +
  scale_x_discrete(limits=rev) +
  coord_flip() +
  facet_grid(y~x, drop = TRUE, space = "free", scales = "free") +
  guides(colour = "none") +
  labs(y = "Hazard ratio", x = "Level") +
  labs(caption = "* p<0.001") + 
  scale_alpha_identity() +
  theme_ali() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")
print(p1)
dev.off()
